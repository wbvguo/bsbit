use std::collections::BTreeMap;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Condvar, Mutex};
use std::thread;

use crate::CliError;
use crate::cpu_placement::CpuPlacement;

const NO_ACTIVE_BATCH: u64 = u64::MAX;

pub(crate) enum WorkerOutcome<T> {
    Completed(T),
    Failed(CliError),
    Cancelled,
}

pub(crate) enum ProducerOutcome {
    Completed,
    Failed(CliError),
    Cancelled,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum DispatchError {
    Cancelled,
    Disconnected { ordinal: u64 },
    OrdinalOverflow,
}

pub(crate) struct WorkDispatcher<T> {
    sender: SyncSender<Indexed<T>>,
    cancellation: Arc<AtomicBool>,
    sent: Arc<AtomicU64>,
    progress: Arc<DispatchProgress>,
    maximum_in_flight: u64,
    next_ordinal: u64,
}

struct DispatchProgress {
    completed: Mutex<u64>,
    changed: Condvar,
}

impl DispatchProgress {
    fn new() -> Self {
        Self {
            completed: Mutex::new(0),
            changed: Condvar::new(),
        }
    }

    fn publish(&self, completed: u64) {
        *self
            .completed
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = completed;
        self.changed.notify_all();
    }

    fn wake(&self) {
        self.changed.notify_all();
    }
}

impl<T> WorkDispatcher<T> {
    fn new(
        sender: SyncSender<Indexed<T>>,
        cancellation: Arc<AtomicBool>,
        sent: Arc<AtomicU64>,
        progress: Arc<DispatchProgress>,
        maximum_in_flight: u64,
    ) -> Self {
        Self {
            sender,
            cancellation,
            sent,
            progress,
            maximum_in_flight,
            next_ordinal: 0,
        }
    }

    pub(crate) fn send(&mut self, work: T) -> Result<u64, DispatchError> {
        if self.cancellation.load(Ordering::Relaxed) {
            return Err(DispatchError::Cancelled);
        }
        let ordinal = self.next_ordinal;
        let mut completed = self
            .progress
            .completed
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        while ordinal.saturating_sub(*completed) >= self.maximum_in_flight {
            if self.cancellation.load(Ordering::Relaxed) {
                return Err(DispatchError::Cancelled);
            }
            completed = self
                .progress
                .changed
                .wait(completed)
                .unwrap_or_else(std::sync::PoisonError::into_inner);
        }
        drop(completed);
        self.sender
            .send(Indexed {
                ordinal,
                value: work,
            })
            .map_err(|_| DispatchError::Disconnected { ordinal })?;
        self.next_ordinal = ordinal
            .checked_add(1)
            .ok_or(DispatchError::OrdinalOverflow)?;
        self.sent.store(self.next_ordinal, Ordering::Relaxed);
        Ok(ordinal)
    }

    #[must_use]
    pub(crate) fn is_cancelled(&self) -> bool {
        self.cancellation.load(Ordering::Relaxed)
    }
}

struct Indexed<T> {
    ordinal: u64,
    value: T,
}

enum WorkerMessage<T> {
    Completed(T),
    Failed(CliError),
    Cancelled,
    Panicked { worker: usize },
}

enum ProducerReady {
    Prepared,
    Failed(CliError),
    Panicked,
}

enum ProducerStart {
    Start,
}

enum ProducerExit {
    Completed { sent: u64 },
    Failed { error: CliError },
    Cancelled { sent: u64 },
    Panicked { sent: u64 },
    ReportedBeforeStart,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn run_ordered_parallel<
    Work,
    Output,
    Prepared,
    Sink,
    ResultValue,
    Prepare,
    Produce,
    Map,
    SetupSink,
    Consume,
    Finish,
>(
    workers: usize,
    prepare: Prepare,
    produce: Produce,
    map: Map,
    setup_sink: SetupSink,
    mut consume: Consume,
    finish: Finish,
) -> Result<ResultValue, CliError>
where
    Work: Send,
    Output: Send,
    Prepare: FnOnce() -> Result<Prepared, CliError> + Send,
    Produce: FnOnce(Prepared, &mut WorkDispatcher<Work>, &AtomicBool) -> ProducerOutcome + Send,
    Map: Fn(usize, Work, &AtomicBool) -> WorkerOutcome<Output> + Sync,
    SetupSink: FnOnce() -> Result<Sink, CliError>,
    Consume: FnMut(&mut Sink, Output) -> Result<(), CliError>,
    Finish: FnOnce(Sink) -> Result<ResultValue, CliError>,
{
    if workers < 2 {
        return Err(CliError::operation(
            "align: internal parallel processing requires at least two workers",
        ));
    }
    let cpu_placement = CpuPlacement::detect(workers);
    thread::scope(|scope| {
        let cancellation = Arc::new(AtomicBool::new(false));
        let sent = Arc::new(AtomicU64::new(0));
        let progress = Arc::new(DispatchProgress::new());
        // Batches have uneven repeat-search cost, so idle workers pull from
        // one bounded queue instead of owning a static ordinal lane.
        let (work_sender, work_receiver) = sync_channel(workers);
        let work_receiver = Arc::new(Mutex::new(work_receiver));
        let (result_sender, result_receiver) = sync_channel(workers);
        let mut worker_handles = Vec::new();
        let mut active_ordinals = Vec::new();
        worker_handles
            .try_reserve_exact(workers)
            .map_err(|_| CliError::operation("align: allocate parallel worker handles: failed"))?;
        active_ordinals
            .try_reserve_exact(workers)
            .map_err(|_| CliError::operation("align: allocate parallel worker state: failed"))?;

        for worker in 0..workers {
            let work_receiver = Arc::clone(&work_receiver);
            let result_sender = result_sender.clone();
            let worker_cancellation = Arc::clone(&cancellation);
            let active = Arc::new(AtomicU64::new(NO_ACTIVE_BATCH));
            active_ordinals.push(Arc::clone(&active));
            let map = &map;
            let worker_cpu_placement = &cpu_placement;
            worker_handles.push((
                worker,
                scope.spawn(move || {
                    worker_cpu_placement.pin_mapping_worker(worker);
                    worker_loop(
                        worker,
                        &work_receiver,
                        &result_sender,
                        &worker_cancellation,
                        &active,
                        map,
                    );
                }),
            ));
        }
        drop(result_sender);

        let (ready_sender, ready_receiver) = sync_channel(0);
        let (start_sender, start_receiver) = sync_channel(0);
        let producer_cancellation = Arc::clone(&cancellation);
        let producer_sent = Arc::clone(&sent);
        let producer_progress = Arc::clone(&progress);
        let producer_cpu_placement = &cpu_placement;
        let producer_handle = scope.spawn(move || {
            producer_cpu_placement.pin_auxiliary_worker();
            let prepared = match catch_unwind(AssertUnwindSafe(prepare)) {
                Ok(Ok(prepared)) => prepared,
                Ok(Err(error)) => {
                    let _ = ready_sender.send(ProducerReady::Failed(error));
                    return ProducerExit::ReportedBeforeStart;
                }
                Err(_) => {
                    let _ = ready_sender.send(ProducerReady::Panicked);
                    return ProducerExit::ReportedBeforeStart;
                }
            };
            if ready_sender.send(ProducerReady::Prepared).is_err() || start_receiver.recv().is_err()
            {
                return ProducerExit::Cancelled {
                    sent: producer_sent.load(Ordering::Relaxed),
                };
            }
            let mut dispatcher = WorkDispatcher::new(
                work_sender,
                Arc::clone(&producer_cancellation),
                Arc::clone(&producer_sent),
                producer_progress,
                u64::try_from(workers.saturating_mul(2)).expect("worker window fits u64"),
            );
            match catch_unwind(AssertUnwindSafe(|| {
                produce(prepared, &mut dispatcher, &producer_cancellation)
            })) {
                Ok(ProducerOutcome::Completed) => ProducerExit::Completed {
                    sent: dispatcher.next_ordinal,
                },
                Ok(ProducerOutcome::Failed(error)) => ProducerExit::Failed { error },
                Ok(ProducerOutcome::Cancelled) => ProducerExit::Cancelled {
                    sent: dispatcher.next_ordinal,
                },
                Err(_) => ProducerExit::Panicked {
                    sent: dispatcher.next_ordinal,
                },
            }
        });

        let ready = match ready_receiver.recv() {
            Ok(ready) => ready,
            Err(_) => ProducerReady::Panicked,
        };
        let ready_error = match ready {
            ProducerReady::Prepared => None,
            ProducerReady::Failed(error) => Some(error),
            ProducerReady::Panicked => Some(CliError::operation(
                "align: parallel input producer panicked before opening input",
            )),
        };
        if let Some(error) = ready_error {
            cancellation.store(true, Ordering::Relaxed);
            progress.wake();
            drop(start_sender);
            drop(result_receiver);
            let _ = producer_handle.join();
            for (_, handle) in worker_handles {
                let _ = handle.join();
            }
            return Err(error);
        }

        let _coordinator_affinity = cpu_placement.pin_auxiliary_scoped();
        let mut sink = match setup_sink() {
            Ok(sink) => sink,
            Err(error) => {
                cancellation.store(true, Ordering::Relaxed);
                progress.wake();
                drop(start_sender);
                drop(result_receiver);
                let _ = producer_handle.join();
                for (_, handle) in worker_handles {
                    let _ = handle.join();
                }
                return Err(error);
            }
        };
        start_sender
            .send(ProducerStart::Start)
            .map_err(|_| CliError::operation("align: parallel input producer disconnected"))?;
        drop(start_sender);

        let mut next_ordinal = 0_u64;
        let mut primary_error = None;
        let mut saw_cancellation = false;
        // Completion may be out of order.  The dispatcher limits the distance
        // ahead to twice the worker count, which bounds this reorder buffer.
        let mut pending = BTreeMap::new();
        while let Ok(message) = result_receiver.recv() {
            if message.ordinal < next_ordinal
                || pending.insert(message.ordinal, message.value).is_some()
            {
                cancellation.store(true, Ordering::Relaxed);
                progress.wake();
                primary_error.get_or_insert_with(|| {
                    CliError::operation(format!(
                        "align: parallel result ordinal {} was completed more than once",
                        message.ordinal
                    ))
                });
                continue;
            }
            while let Some(message) = pending.remove(&next_ordinal) {
                match message {
                    WorkerMessage::Completed(output) if primary_error.is_none() => {
                        if let Err(error) = consume(&mut sink, output) {
                            cancellation.store(true, Ordering::Relaxed);
                            progress.wake();
                            primary_error = Some(error);
                        }
                    }
                    WorkerMessage::Completed(_) => {}
                    WorkerMessage::Failed(error) => {
                        cancellation.store(true, Ordering::Relaxed);
                        progress.wake();
                        primary_error.get_or_insert(error);
                    }
                    WorkerMessage::Cancelled => saw_cancellation = true,
                    WorkerMessage::Panicked { worker } => {
                        cancellation.store(true, Ordering::Relaxed);
                        progress.wake();
                        primary_error.get_or_insert_with(|| {
                            CliError::operation(format!(
                                "align: parallel worker {worker} panicked in batch {next_ordinal}"
                            ))
                        });
                    }
                }
                let Some(completed) = next_ordinal.checked_add(1) else {
                    cancellation.store(true, Ordering::Relaxed);
                    progress.wake();
                    primary_error.get_or_insert_with(|| {
                        CliError::operation("align: parallel completed-batch ordinal overflow")
                    });
                    break;
                };
                next_ordinal = completed;
                progress.publish(next_ordinal);
            }
        }
        if primary_error.is_some() {
            cancellation.store(true, Ordering::Relaxed);
            progress.wake();
        }
        drop(result_receiver);

        let producer_exit = match producer_handle.join() {
            Ok(exit) => exit,
            Err(_) => ProducerExit::Panicked {
                sent: sent.load(Ordering::Relaxed),
            },
        };
        let mut outer_worker_panic = None;
        for (worker, handle) in worker_handles {
            if handle.join().is_err() && outer_worker_panic.is_none() {
                let active = active_ordinals[worker].load(Ordering::Relaxed);
                outer_worker_panic = Some((worker, active));
            }
        }

        if let Some(error) = primary_error {
            return Err(error);
        }
        if let Some((worker, active)) = outer_worker_panic {
            let context = if active == NO_ACTIVE_BATCH {
                format!("align: parallel worker {worker} panicked outside a batch")
            } else {
                format!("align: parallel worker {worker} panicked in batch {active}")
            };
            return Err(CliError::operation(context));
        }
        let sent_batches = match producer_exit {
            ProducerExit::Completed { sent } => sent,
            ProducerExit::Failed { error } => return Err(error),
            ProducerExit::Panicked { sent } => {
                return Err(CliError::operation(format!(
                    "align: parallel input producer panicked before batch {sent}"
                )));
            }
            ProducerExit::Cancelled { sent } => {
                return Err(CliError::operation(format!(
                    "align: parallel processing cancelled without a primary error after {sent} batches"
                )));
            }
            ProducerExit::ReportedBeforeStart => {
                return Err(CliError::operation(
                    "align: parallel input producer stopped before start",
                ));
            }
        };
        if saw_cancellation {
            return Err(CliError::operation(
                "align: parallel worker cancelled without a primary error",
            ));
        }
        if next_ordinal != sent_batches {
            return Err(CliError::operation(format!(
                "align: parallel result count expected {sent_batches} batches, observed {next_ordinal}"
            )));
        }
        finish(sink)
    })
}

fn worker_loop<Work, Output, Map>(
    worker: usize,
    work_receiver: &Arc<Mutex<Receiver<Indexed<Work>>>>,
    result_sender: &SyncSender<Indexed<WorkerMessage<Output>>>,
    cancellation: &AtomicBool,
    active: &AtomicU64,
    map: &Map,
) where
    Map: Fn(usize, Work, &AtomicBool) -> WorkerOutcome<Output>,
{
    loop {
        let work = work_receiver
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .recv();
        let Ok(work) = work else {
            break;
        };
        active.store(work.ordinal, Ordering::Relaxed);
        let message = if cancellation.load(Ordering::Relaxed) {
            WorkerMessage::Cancelled
        } else {
            match catch_unwind(AssertUnwindSafe(|| map(worker, work.value, cancellation))) {
                Ok(WorkerOutcome::Completed(output)) => {
                    if cancellation.load(Ordering::Relaxed) {
                        WorkerMessage::Cancelled
                    } else {
                        WorkerMessage::Completed(output)
                    }
                }
                Ok(WorkerOutcome::Failed(error)) => WorkerMessage::Failed(error),
                Ok(WorkerOutcome::Cancelled) => WorkerMessage::Cancelled,
                Err(_) => WorkerMessage::Panicked { worker },
            }
        };
        let ordinal = work.ordinal;
        if result_sender
            .send(Indexed {
                ordinal,
                value: message,
            })
            .is_err()
        {
            active.store(NO_ACTIVE_BATCH, Ordering::Relaxed);
            return;
        }
        active.store(NO_ACTIVE_BATCH, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Barrier;
    use std::time::Duration;

    #[test]
    fn slow_first_lane_backpressures_and_order_is_exact() {
        let release = Arc::new(AtomicBool::new(false));
        let sent = Arc::new(AtomicU64::new(0));
        let entered = Arc::new(Barrier::new(2));
        let thread_release = Arc::clone(&release);
        let thread_sent = Arc::clone(&sent);
        let thread_entered = Arc::clone(&entered);
        let handle = thread::spawn(move || {
            run_ordered_parallel(
                2,
                || Ok(()),
                move |(), dispatcher, _| {
                    for value in 0..100_u64 {
                        match dispatcher.send(value) {
                            Ok(_) => {
                                thread_sent.fetch_add(1, Ordering::Relaxed);
                            }
                            Err(DispatchError::Cancelled) => return ProducerOutcome::Cancelled,
                            Err(error) => {
                                return ProducerOutcome::Failed(CliError::operation(format!(
                                    "dispatch failed: {error:?}"
                                )));
                            }
                        }
                    }
                    ProducerOutcome::Completed
                },
                move |_, value, _| {
                    if value == 0 {
                        thread_entered.wait();
                        while !thread_release.load(Ordering::Relaxed) {
                            thread::yield_now();
                        }
                    }
                    WorkerOutcome::Completed(value)
                },
                || Ok(Vec::new()),
                |output, value| {
                    output.push(value);
                    Ok(())
                },
                Ok,
            )
        });
        entered.wait();
        while sent.load(Ordering::Relaxed) < 4 {
            thread::yield_now();
        }
        thread::sleep(Duration::from_millis(20));
        assert_eq!(sent.load(Ordering::Relaxed), 4);
        release.store(true, Ordering::Relaxed);
        let output = handle
            .join()
            .expect("scheduler thread joins")
            .expect("parallel processing succeeds");
        assert_eq!(output, (0..100_u64).collect::<Vec<_>>());
    }

    #[test]
    fn worker_panic_is_contained_and_sink_is_not_finished() {
        let finished = Arc::new(AtomicBool::new(false));
        let observed_finished = Arc::clone(&finished);
        let error = run_ordered_parallel(
            2,
            || Ok(()),
            |(), dispatcher, _| {
                for value in 0..8_u64 {
                    if dispatcher.send(value).is_err() {
                        return ProducerOutcome::Cancelled;
                    }
                }
                ProducerOutcome::Completed
            },
            |_, value, _| {
                assert!(value != 2, "injected worker panic");
                WorkerOutcome::Completed(value)
            },
            || Ok(Vec::new()),
            |output, value| {
                output.push(value);
                Ok(())
            },
            move |output| {
                observed_finished.store(true, Ordering::Relaxed);
                Ok(output)
            },
        )
        .expect_err("panic fails parallel processing");
        let message = error.to_string();
        assert!(message.starts_with("align: parallel worker "));
        assert!(message.ends_with(" panicked in batch 2"));
        assert!(!finished.load(Ordering::Relaxed));
    }

    #[test]
    fn earlier_worker_error_outranks_later_producer_error() {
        let error = run_ordered_parallel(
            2,
            || Ok(()),
            |(), dispatcher, _| {
                dispatcher.send(0_u64).expect("first batch sends");
                ProducerOutcome::Failed(CliError::operation("injected producer failure"))
            },
            |_, _, _| {
                WorkerOutcome::<u64>::Failed(CliError::operation("injected batch-zero failure"))
            },
            || Ok(Vec::new()),
            |output, value| {
                output.push(value);
                Ok(())
            },
            Ok,
        )
        .expect_err("earlier worker failure rejects parallel processing");
        assert_eq!(error.to_string(), "injected batch-zero failure");
    }

    #[test]
    fn lowest_worker_ordinal_wins_independent_of_completion_order() {
        let release = Arc::new(AtomicBool::new(false));
        let entered = Arc::new(Barrier::new(2));
        let map_release = Arc::clone(&release);
        let map_entered = Arc::clone(&entered);
        let error = thread::scope(|scope| {
            let handle = scope.spawn(|| {
                run_ordered_parallel(
                    2,
                    || Ok(()),
                    |(), dispatcher, _| {
                        dispatcher.send(0_u64).expect("batch zero sends");
                        dispatcher.send(1_u64).expect("batch one sends");
                        ProducerOutcome::Completed
                    },
                    move |_, value, _| {
                        if value == 0 {
                            map_entered.wait();
                            while !map_release.load(Ordering::Relaxed) {
                                thread::yield_now();
                            }
                        } else {
                            map_release.store(true, Ordering::Relaxed);
                        }
                        WorkerOutcome::<u64>::Failed(CliError::operation(format!(
                            "injected batch-{value} failure"
                        )))
                    },
                    || Ok(Vec::new()),
                    |output, value| {
                        output.push(value);
                        Ok(())
                    },
                    Ok,
                )
            });
            entered.wait();
            handle.join().expect("scheduler thread joins")
        })
        .expect_err("worker failure rejects parallel processing");
        assert_eq!(error.to_string(), "injected batch-0 failure");
    }

    #[test]
    fn every_parallel_worker_count_preserves_order_under_yield_perturbation() {
        for workers in 2..=64_usize {
            let item_count = u64::try_from(workers * 4).expect("small item count");
            for rotation in 0..3_usize {
                let output = run_ordered_parallel(
                    workers,
                    || Ok(()),
                    |(), dispatcher, _| {
                        for value in 0..item_count {
                            if dispatcher.send(value).is_err() {
                                return ProducerOutcome::Cancelled;
                            }
                        }
                        ProducerOutcome::Completed
                    },
                    |worker, value, _| {
                        for _ in 0..((workers + rotation - worker) % workers) {
                            thread::yield_now();
                        }
                        WorkerOutcome::Completed(value)
                    },
                    || Ok(Vec::new()),
                    |output, value| {
                        output.push(value);
                        Ok(())
                    },
                    Ok,
                )
                .expect("perturbed schedule succeeds");
                assert_eq!(
                    output,
                    (0..item_count).collect::<Vec<_>>(),
                    "order differs with {workers} workers in rotation {rotation}"
                );
            }
        }
    }

    #[test]
    fn sink_failure_cancels_workers_and_never_finishes() {
        let finished = Arc::new(AtomicBool::new(false));
        let observed_finished = Arc::clone(&finished);
        let error = run_ordered_parallel(
            4,
            || Ok(()),
            |(), dispatcher, _| {
                for value in 0..100_u64 {
                    if dispatcher.send(value).is_err() {
                        return ProducerOutcome::Cancelled;
                    }
                }
                ProducerOutcome::Completed
            },
            |_, value, _| WorkerOutcome::Completed(value),
            || Ok(Vec::new()),
            |output, value| {
                if value == 2 {
                    Err(CliError::operation("injected sink failure"))
                } else {
                    output.push(value);
                    Ok(())
                }
            },
            move |output| {
                observed_finished.store(true, Ordering::Relaxed);
                Ok(output)
            },
        )
        .expect_err("sink failure rejects parallel processing");
        assert_eq!(error.to_string(), "injected sink failure");
        assert!(!finished.load(Ordering::Relaxed));
    }

    #[test]
    fn producer_panic_is_contained_and_never_finishes() {
        let finished = Arc::new(AtomicBool::new(false));
        let observed_finished = Arc::clone(&finished);
        let error = run_ordered_parallel(
            2,
            || Ok(()),
            |(), _, _| -> ProducerOutcome { panic!("injected producer panic") },
            |_, value: u64, _| WorkerOutcome::Completed(value),
            || Ok(Vec::new()),
            |output, value| {
                output.push(value);
                Ok(())
            },
            move |output| {
                observed_finished.store(true, Ordering::Relaxed);
                Ok(output)
            },
        )
        .expect_err("producer panic rejects parallel processing");
        assert_eq!(
            error.to_string(),
            "align: parallel input producer panicked before batch 0"
        );
        assert!(!finished.load(Ordering::Relaxed));
    }

    #[test]
    fn sink_setup_failure_never_starts_the_producer_or_finishes() {
        let produced = Arc::new(AtomicBool::new(false));
        let finished = Arc::new(AtomicBool::new(false));
        let observed_produced = Arc::clone(&produced);
        let observed_finished = Arc::clone(&finished);
        let error = run_ordered_parallel(
            2,
            || Ok(()),
            move |(), _, _| {
                observed_produced.store(true, Ordering::Relaxed);
                ProducerOutcome::Completed
            },
            |_, value: u64, _| WorkerOutcome::Completed(value),
            || Err::<Vec<u64>, _>(CliError::operation("injected sink setup failure")),
            |output, value| {
                output.push(value);
                Ok(())
            },
            move |output| {
                observed_finished.store(true, Ordering::Relaxed);
                Ok(output)
            },
        )
        .expect_err("sink setup failure rejects parallel processing");
        assert_eq!(error.to_string(), "injected sink setup failure");
        assert!(!produced.load(Ordering::Relaxed));
        assert!(!finished.load(Ordering::Relaxed));
    }
}
