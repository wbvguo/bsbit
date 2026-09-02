//! Public-contract tests for identity-safe generic file publication.

use std::fs;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Barrier};
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_io::{
    PublicationPhase, StagedFile, reopen_read_write, select_sibling_staging_path,
    validate_create_target, validate_distinct_paths, validate_regular_file_or_absent,
    validate_replace_target,
};

fn unique_path(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-io-publication-{label}-{}-{nonce}",
        std::process::id()
    ))
}

fn complete_bytes(target: &Path, label: &str, bytes: &[u8]) -> bsbit_io::CompletedFile {
    let mut staged = StagedFile::create_sibling(target, label).expect("stage");
    let mut file = staged.take_file().expect("descriptor");
    file.write_all(bytes).expect("write");
    staged.complete(file).expect("complete")
}

fn complete_replacement_bytes(target: &Path, label: &str, bytes: &[u8]) -> bsbit_io::CompletedFile {
    let mut staged = StagedFile::create_sibling_replace(target, label).expect("stage replacement");
    let mut file = staged.take_file().expect("descriptor");
    file.write_all(bytes).expect("write");
    staged.complete(file).expect("complete")
}

#[test]
fn generic_create_target_validation_owns_parent_and_absence_policy() {
    let directory = unique_path("create-target-directory");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("result.dat");
    validate_create_target(&target).expect("absent target with directory parent");

    fs::write(&target, b"owner").expect("occupied target");
    assert_eq!(
        validate_create_target(&target)
            .expect_err("occupied target fails")
            .kind(),
        std::io::ErrorKind::AlreadyExists
    );
    fs::remove_file(&target).expect("target cleanup");

    let parent_file = directory.join("not-a-directory");
    fs::write(&parent_file, b"file").expect("parent fixture");
    assert_eq!(
        validate_create_target(&parent_file.join("result.dat"))
            .expect_err("file parent fails")
            .kind(),
        std::io::ErrorKind::NotADirectory
    );
    fs::remove_file(parent_file).expect("parent fixture cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn generic_reader_path_validation_rejects_non_files_but_preserves_missing_paths() {
    let directory = unique_path("reader-path-directory");
    fs::create_dir(&directory).expect("directory");
    let file = directory.join("input.dat");
    let missing = directory.join("missing.dat");

    fs::write(&file, b"input").expect("file fixture");
    validate_regular_file_or_absent(&file).expect("regular file");
    validate_regular_file_or_absent(&missing).expect("missing path passes to format opener");
    assert_eq!(
        validate_regular_file_or_absent(&directory)
            .expect_err("directory is not a file")
            .kind(),
        std::io::ErrorKind::Unsupported
    );

    fs::remove_file(file).expect("file cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn generic_replace_target_validation_accepts_files_and_rejects_directories() {
    let directory = unique_path("replace-target-directory");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("result.dat");
    validate_replace_target(&target).expect("missing target");
    fs::write(&target, b"existing").expect("target");
    validate_replace_target(&target).expect("regular target");
    assert_eq!(
        validate_replace_target(&directory)
            .expect_err("directory target fails")
            .kind(),
        std::io::ErrorKind::Unsupported
    );
    fs::remove_file(target).expect("target cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn selected_staging_candidates_are_absolute_unused_unique_siblings() {
    let directory = unique_path("selected-staging-directory");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("result.dat");
    let first = select_sibling_staging_path(&target, "index/unsafe?").expect("first candidate");
    let second = select_sibling_staging_path(&target, "index/unsafe?").expect("second candidate");

    assert!(first.is_absolute());
    assert_eq!(first.parent(), target.parent());
    assert_eq!(second.parent(), target.parent());
    assert_ne!(first, second);
    assert!(!first.exists());
    assert!(!second.exists());
    assert!(!target.exists());
    let first_name = first
        .file_name()
        .and_then(|name| name.to_str())
        .expect("UTF-8 generated name");
    assert!(first_name.starts_with(".bsbit-indexunsafe-"));
    assert_eq!(first.extension(), Some(std::ffi::OsStr::new("tmp")));
    validate_distinct_paths(&first, &target).expect("candidate and target differ");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn completed_bytes_publish_create_only_and_can_roll_back() {
    let target = unique_path("target");
    let mut staged = StagedFile::create_sibling(&target, "contract").expect("stage");
    let mut file = staged.take_file().expect("descriptor");
    file.write_all(b"complete bytes").expect("write");
    let completed = staged.complete(file).expect("complete");
    let published = match completed.publish_create_new() {
        Ok(published) => published,
        Err(error) if error.kind() == std::io::ErrorKind::Unsupported => return,
        Err(error) => panic!("publish: {error}"),
    };
    assert_eq!(
        fs::read(&target).expect("published bytes"),
        b"complete bytes"
    );
    published.rollback().expect("rollback");
    assert!(!target.exists());
}

#[test]
fn completed_bytes_replace_atomically_and_rollback_restores_old_target() {
    let directory = unique_path("replace-rollback-directory");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("result.dat");
    fs::write(&target, b"old bytes").expect("old target");
    let completed = complete_replacement_bytes(&target, "replace", b"new bytes");
    let published = completed.publish_replace().expect("replace target");
    assert_eq!(fs::read(&target).expect("new target"), b"new bytes");
    published.rollback().expect("restore old target");
    assert_eq!(fs::read(&target).expect("restored target"), b"old bytes");
    assert_eq!(
        fs::read_dir(&directory).expect("directory entries").count(),
        1
    );
    fs::remove_file(target).expect("target cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn successful_replacement_drop_removes_private_backup() {
    let directory = unique_path("replace-commit-directory");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("result.dat");
    fs::write(&target, b"old bytes").expect("old target");
    let completed = complete_replacement_bytes(&target, "replace", b"new bytes");
    let published = completed.publish_replace().expect("replace target");
    drop(published);
    assert_eq!(fs::read(&target).expect("new target"), b"new bytes");
    assert_eq!(
        fs::read_dir(&directory).expect("directory entries").count(),
        1
    );
    fs::remove_file(target).expect("target cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn replacement_of_staging_is_never_removed_or_published() {
    let target = unique_path("identity-target");
    let mut staged = StagedFile::create_sibling(&target, "identity").expect("stage");
    let staging = staged.path().to_path_buf();
    let displaced = unique_path("displaced");
    let file = staged.take_file().expect("descriptor");
    fs::rename(&staging, &displaced).expect("displace owned path");
    fs::write(&staging, b"replacement").expect("replacement");
    let error = staged.complete(file).expect_err("identity mismatch");
    assert_eq!(error.phase(), PublicationPhase::ValidateStaging);
    assert_eq!(error.cleanup_warning(), None);
    assert_eq!(
        fs::read(&staging).expect("replacement survives"),
        b"replacement"
    );
    assert!(!target.exists());
    fs::remove_file(staging).expect("replacement cleanup");
    fs::remove_file(displaced).expect("owned file cleanup");
}

#[test]
fn independent_descriptor_cursor_stays_bound_after_path_replacement() {
    let target = unique_path("reopen-target");
    let mut staged = StagedFile::create_sibling(&target, "reopen").expect("stage");
    let staging = staged.path().to_path_buf();
    let displaced = unique_path("reopen-displaced");
    let mut owner = staged.take_file().expect("owned descriptor");
    owner.write_all(b"owner-000").expect("initial owner bytes");
    let mut independent = reopen_read_write(&owner).expect("independent cursor");
    independent
        .seek(SeekFrom::Start(6))
        .expect("independent seek");

    fs::rename(&staging, &displaced).expect("displace owned namespace entry");
    fs::write(&staging, b"replacement").expect("replacement entry");
    independent.write_all(b"123").expect("write held object");
    independent.sync_all().expect("sync held object");

    assert_eq!(fs::read(&displaced).expect("owned bytes"), b"owner-123");
    assert_eq!(
        fs::read(&staging).expect("replacement bytes"),
        b"replacement"
    );
    drop(independent);
    drop(owner);
    drop(staged);
    fs::remove_file(staging).expect("replacement cleanup");
    fs::remove_file(displaced).expect("owned cleanup");
}

#[test]
fn existing_target_wins_without_modification() {
    let target = unique_path("existing-target");
    fs::write(&target, b"owner bytes").expect("target");
    let error = StagedFile::create_sibling(&target, "existing").expect_err("create-only");
    assert_eq!(error.phase(), PublicationPhase::ValidatePaths);
    assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(fs::read(&target).expect("target survives"), b"owner bytes");
    fs::remove_file(target).expect("cleanup");
}

#[test]
fn target_created_after_staging_wins_and_private_bytes_are_removed() {
    let target = unique_path("late-target");
    let completed = complete_bytes(&target, "late-target", b"private bytes");
    let staging = completed.staging_path().to_path_buf();
    fs::write(&target, b"concurrent owner").expect("late target");

    let error = completed
        .publish_create_new()
        .expect_err("late target wins");
    assert_eq!(error.phase(), PublicationPhase::ValidatePaths);
    assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(
        fs::read(&target).expect("target survives"),
        b"concurrent owner"
    );
    assert!(!staging.exists());
    fs::remove_file(target).expect("target cleanup");
}

#[test]
fn concurrent_publications_have_exactly_one_winner() {
    const WORKERS: usize = 8;

    let target = unique_path("concurrent-target");
    let mut completed = Vec::new();
    let mut staging_paths = Vec::new();
    for worker in 0..WORKERS {
        let value = complete_bytes(&target, "concurrent", format!("worker-{worker}").as_bytes());
        staging_paths.push(value.staging_path().to_path_buf());
        completed.push(value);
    }
    let barrier = Arc::new(Barrier::new(WORKERS));
    let handles: Vec<_> = completed
        .into_iter()
        .map(|value| {
            let barrier = Arc::clone(&barrier);
            std::thread::spawn(move || {
                barrier.wait();
                value.publish_create_new()
            })
        })
        .collect();

    let mut successes = 0;
    let mut occupied = 0;
    let mut unsupported = 0;
    for handle in handles {
        match handle.join().expect("publisher thread") {
            Ok(publication) => {
                successes += 1;
                assert_eq!(publication.cleanup_warning(), None);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => occupied += 1,
            Err(error) if error.kind() == std::io::ErrorKind::Unsupported => unsupported += 1,
            Err(error) => panic!("unexpected publication error: {error}"),
        }
    }
    if unsupported == WORKERS {
        assert_eq!(successes, 0);
        assert_eq!(occupied, 0);
    } else {
        assert_eq!(unsupported, 0);
        assert_eq!(successes, 1);
        assert_eq!(occupied, WORKERS - 1);
        fs::remove_file(&target).expect("winner cleanup");
    }
    assert!(staging_paths.iter().all(|path| !path.exists()));
}

#[cfg(unix)]
#[test]
fn dangling_target_symlink_counts_as_an_existing_entry() {
    use std::os::unix::fs::symlink;

    let target = unique_path("dangling-target");
    let missing = unique_path("missing-target");
    symlink(&missing, &target).expect("dangling target symlink");
    let error =
        StagedFile::create_sibling(&target, "dangling").expect_err("symlink is an occupied target");
    assert_eq!(error.phase(), PublicationPhase::ValidatePaths);
    assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(fs::read_link(&target).expect("symlink survives"), missing);
    fs::remove_file(target).expect("symlink cleanup");
}
