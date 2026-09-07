# Alignment performance matrix

## Run identity

| Field | Value |
|---|---|
| bsbit benchmark date | 2026-09-06（America/Los_Angeles） |
| External-tool matrix date | 2026-09-05（同一主机、输入和资源口径） |
| bsbit software | `0.1.0` release candidate；最终来源由 Git tag `v0.1.0` 固定 |
| Alignment policy | `bounded-structural-alignment-v1` |
| MAPQ policy | `structural-origin-evidence-v1` |
| Metrics schemas | `bsbit-alignment-metrics-single-end-v1`、`bsbit-alignment-metrics-paired-end-v1` |
| bsbit benchmark binary SHA-256 | `b32acf7fd02a2a8373bbc27e9b8937c1d5393f130010391282e7dd8baf61a28b` |
| fast index SHA-256 | `f1d2d2a876b5721f7f86c16649cce6c9432593cc610e6244b359b99d9affb53a` |
| Performance input | 5,000,000 reads（SE）或 read pairs（PE） |
| Accuracy input | 200,000 个唯一 truth reads/pairs |
| Host | Intel Core i7-14700K，Linux/WSL，8 mapping workers |

四个实验分别成表，directional 与 non-directional 不混表。本页保存在
`docs/development/` 作为可审计开发证据，但不发布到 MkDocs 网站。

粗体表示该实验内原始数值最佳：Wall、CPU time、RSS 越低越好，其余越高越好。Q 列依次为 precision / recall / F1，并以 `+/-5 bp` 正确性计算；F1 后的 `*` 表示该结果没有通过相应 Q 阈值的单侧 95% Clopper–Pearson 错误率上界。加粗依据未四舍五入值，因此两个显示为相同四位小数的值不一定同时加粗。

本次重新编译并重跑当前工作树中的 bsbit default 与 sensitive。bsbit 性能值来自 5,000,000 reads（SE）或 read pairs（PE）；准确率来自对应的 200,000 个唯一 truth reads/pairs。其他工具保留自 2026-09-05 的同数据、同硬件科学矩阵，不在本次重复运行。

## SE directional

Throughput 单位：reads/s。

| Tool / mode | Wall | CPU time | RSS | Throughput | Output reads | Mapped | Exact rate | +/-5 bp rate | Q10 precision/recall/F1 | Q20 precision/recall/F1 | Q30 precision/recall/F1 | Q40 precision/recall/F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| bsbit default, current/fast-index | **25.80 s** | 207.83 s | 8.95 GiB | **193,833** | **100.000%** | 99.866% | 88.981% | 89.012% | 0.9985 / 0.8369 / 0.9106 | 0.9994 / 0.7825 / 0.8778 | 0.9998 / 0.6385 / 0.7793 | 1.0000 / 0.5768 / 0.7316 |
| bsbit sensitive, current/fast-index | 131.24 s | 1042.77 s | 8.95 GiB | 38,099 | **100.000%** | 99.928% | 89.001% | 89.040% | 0.9992 / **0.8454** / **0.9159** | 0.9997 / 0.8389 / **0.9123** | 0.9999 / 0.7964 / 0.8866 | 1.0000 / 0.7760 / 0.8739 |
| BitMapperBS | 48.90 s | **170.92 s** | **6.77 GiB** | 102,240 | **100.000%** | 98.978% | 87.902% | 87.934% | 0.9993 / 0.8021 / 0.8899 | 0.9994 / 0.8011 / 0.8893 | 0.9995 / 0.6821 / 0.8108 | 0.9996 / 0.6817 / 0.8106\* |
| HISAT-3N | 298.68 s | 1202.69 s | 7.90 GiB | 16,740 | 84.519% | 84.519% | 84.084% | 84.162% | 0.9958 / 0.8416 / 0.9122 | 0.9958 / **0.8416** / 0.9122 | 0.9958 / **0.8416** / **0.9122\*** | 0.9958 / **0.8416** / **0.9122\*** |
| Bismark | 585.34 s | 4181.91 s | 11.07 GiB | 8,542 | 84.183% | 84.183% | 84.011% | 84.042% | 0.9995 / 0.8070 / 0.8930 | 0.9996 / 0.7976 / 0.8872 | 0.9996 / 0.6927 / 0.8183 | 1.0000 / 0.6650 / 0.7988 |
| BISCUIT | 591.46 s | 4353.30 s | 11.34 GiB | 8,454 | **100.000%** | **100.000%** | **89.033%** | **89.072%** | **1.0000** / 0.8160 / 0.8986 | **1.0000** / 0.7990 / 0.8883 | **1.0000** / 0.7872 / 0.8809 | **1.0000** / 0.7749 / 0.8732 |
| BSBolt | 474.96 s | 3557.86 s | 13.12 GiB | 10,527 | **100.000%** | 91.213% | 85.047% | 85.077% | 1.0000 / 0.8190 / 0.9005 | 1.0000 / 0.8156 / 0.8984 | 1.0000 / 0.7993 / 0.8884 | 1.0000 / 0.7869 / 0.8808 |

通过校准后的 F1 最优：Q10 bsbit sensitive；Q20 bsbit sensitive；Q30 BSBolt；Q40 BSBolt。

## SE non-directional

Throughput 单位：reads/s。BitMapperBS 不支持 non-directional，因此不列入本组。

| Tool / mode | Wall | CPU time | RSS | Throughput | Output reads | Mapped | Exact rate | +/-5 bp rate | Q10 precision/recall/F1 | Q20 precision/recall/F1 | Q30 precision/recall/F1 | Q40 precision/recall/F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| bsbit default, current/fast-index | **58.55 s** | **466.88 s** | 8.94 GiB | **85,403** | **100.000%** | 99.876% | 88.976% | 89.015% | 0.9987 / 0.8406 / 0.9128 | 0.9995 / 0.7833 / 0.8783 | 0.9999 / 0.6389 / 0.7796 | 1.0000 / 0.5778 / 0.7324 |
| bsbit sensitive, current/fast-index | 336.08 s | 2650.64 s | 8.95 GiB | 14,878 | **100.000%** | 99.932% | **88.994%** | **89.039%** | 0.9992 / **0.8456** / **0.9160** | 0.9998 / 0.8389 / **0.9123** | 0.9999 / 0.7992 / 0.8883 | 1.0000 / 0.7786 / 0.8755 |
| HISAT-3N | 429.07 s | 2357.73 s | **7.90 GiB** | 11,653 | 84.576% | 84.576% | 84.111% | 84.189% | 0.9954 / 0.8419 / 0.9122 | 0.9954 / **0.8419** / 0.9122 | 0.9954 / **0.8419** / **0.9122\*** | 0.9954 / **0.8419** / **0.9122\*** |
| Bismark | 1083.03 s | 7734.20 s | 19.04 GiB | 4,617 | 84.239% | 84.239% | 84.079% | 84.105% | 0.9994 / 0.8080 / 0.8936 | 0.9995 / 0.7983 / 0.8877 | 0.9996 / 0.6928 / 0.8184 | 1.0000 / 0.6653 / 0.7990 |
| BISCUIT | 2251.54 s | 17571.24 s | 15.81 GiB | 2,221 | **100.000%** | **100.000%** | 88.981% | 89.019% | **1.0000** / 0.8163 / 0.8988 | **1.0000** / 0.7991 / 0.8883 | **1.0000** / 0.7875 / 0.8811 | **1.0000** / 0.7750 / 0.8732 |
| BSBolt | 543.83 s | 4092.64 s | 13.24 GiB | 9,194 | **100.000%** | 91.306% | 84.869% | 84.897% | 0.9992 / 0.8176 / 0.8993 | 0.9995 / 0.8141 / 0.8973 | 0.9996 / 0.7981 / 0.8876 | 0.9997 / 0.7860 / 0.8800\* |

通过校准后的 F1 最优：Q10、Q20、Q30、Q40 均为 bsbit sensitive。

## PE directional

Throughput 单位：pairs/s。

| Tool / mode | Wall | CPU time | RSS | Throughput | Output reads | Mapped | Exact rate | +/-5 bp rate | Q10 precision/recall/F1 | Q20 precision/recall/F1 | Q30 precision/recall/F1 | Q40 precision/recall/F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| bsbit default, current/fast-index | **21.90 s** | **181.98 s** | 8.90 GiB | **228,266** | **100.000%** | 99.057% | 90.257% | 90.352% | 0.9996 / 0.8525 / 0.9202 | 0.9996 / 0.8525 / 0.9202 | 0.9996 / 0.8525 / 0.9202 | 0.9996 / 0.8461 / 0.9165\* |
| bsbit sensitive, current/fast-index | 203.18 s | 1589.33 s | 8.88 GiB | 24,608 | **100.000%** | **99.910%** | **90.723%** | **90.787%** | 0.9994 / 0.8629 / **0.9262** | 0.9994 / **0.8626** / **0.9260** | 0.9998 / 0.8519 / 0.9199 | **1.0000** / 0.7795 / 0.8761 |
| BitMapperBS | 68.74 s | 282.01 s | **6.77 GiB** | 72,740 | 84.901% | 84.901% | 84.826% | 84.882% | **0.9998** / 0.8319 / 0.9082 | 0.9998 / 0.8310 / 0.9076 | 0.9999 / 0.6794 / 0.8091 | 0.9999 / 0.6789 / 0.8087\* |
| HISAT-3N | 398.10 s | 2184.28 s | 7.93 GiB | 12,560 | 86.847% | 86.847% | 85.947% | 86.120% | 0.9916 / 0.8612 / 0.9218 | 0.9916 / 0.8612 / 0.9218 | 0.9916 / **0.8612** / **0.9218\*** | 0.9916 / **0.8612** / **0.9218\*** |
| Bismark | 954.63 s | 6899.17 s | 11.16 GiB | 5,238 | 86.741% | 86.741% | 86.478% | 86.541% | 0.9990 / 0.8538 / 0.9207 | 0.9997 / 0.8324 / 0.9084 | 0.9998 / 0.7766 / 0.8742 | 1.0000 / 0.7620 / 0.8649 |
| BISCUIT | 1212.27 s | 9295.18 s | 11.33 GiB | 4,125 | **100.000%** | 99.527% | 90.530% | 90.611% | 0.9992 / **0.8631** / 0.9261 | 0.9995 / 0.8477 / 0.9174 | 0.9996 / 0.8455 / 0.9161 | 0.9997 / 0.8363 / 0.9107\* |
| BSBolt | 1079.59 s | 8396.53 s | 13.12 GiB | 4,631 | **100.000%** | 93.146% | 87.581% | 87.643% | 0.9998 / 0.8199 / 0.9010 | **0.9999** / 0.8197 / 0.9009 | **0.9999** / 0.8143 / 0.8976 | 1.0000 / 0.8130 / 0.8968 |

通过校准后的 F1 最优：Q10 bsbit sensitive；Q20 bsbit sensitive；Q30 bsbit default；Q40 BSBolt。

## PE non-directional

Throughput 单位：pairs/s。BitMapperBS 不支持 non-directional，因此不列入本组。

| Tool / mode | Wall | CPU time | RSS | Throughput | Output reads | Mapped | Exact rate | +/-5 bp rate | Q10 precision/recall/F1 | Q20 precision/recall/F1 | Q30 precision/recall/F1 | Q40 precision/recall/F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| bsbit default, current/fast-index | **44.08 s** | **358.00 s** | 8.85 GiB | **113,426** | **100.000%** | 99.096% | 90.267% | 90.363% | **0.9996** / 0.8527 / 0.9204 | 0.9996 / 0.8527 / 0.9204 | 0.9996 / 0.8527 / 0.9204 | 0.9997 / 0.8462 / 0.9166\* |
| bsbit sensitive, current/fast-index | 321.43 s | 2534.63 s | 8.88 GiB | 15,555 | **100.000%** | **99.920%** | **90.723%** | **90.787%** | 0.9994 / 0.8628 / 0.9261 | 0.9994 / **0.8625** / **0.9260** | **0.9998** / 0.8463 / 0.9167 | **1.0000** / 0.7748 / 0.8731 |
| HISAT-3N | 656.76 s | 4263.74 s | **7.93 GiB** | 7,613 | 86.856% | 86.856% | 85.952% | 86.124% | 0.9916 / 0.8612 / 0.9218 | 0.9916 / 0.8612 / 0.9218 | 0.9916 / **0.8612** / **0.9218\*** | 0.9916 / **0.8612** / **0.9218\*** |
| Bismark | 1412.83 s | 10258.00 s | 19.04 GiB | 3,539 | 86.743% | 86.743% | 86.476% | 86.540% | 0.9989 / 0.8541 / 0.9209 | **0.9998** / 0.8323 / 0.9084 | 0.9998 / 0.7760 / 0.8738 | 1.0000 / 0.7618 / 0.8648 |
| BISCUIT | 4615.43 s | 36530.53 s | 15.79 GiB | 1,083 | **100.000%** | 99.526% | 90.496% | 90.577% | 0.9992 / **0.8631** / **0.9262** | 0.9994 / 0.8477 / 0.9173 | 0.9995 / 0.8455 / 0.9161 | 0.9997 / 0.8364 / 0.9108\* |
| BSBolt | 1107.03 s | 8624.70 s | 13.19 GiB | 4,517 | **100.000%** | 93.072% | 87.477% | 87.540% | 0.9993 / 0.8188 / 0.9001 | 0.9993 / 0.8185 / 0.8999 | 0.9994 / 0.8132 / 0.8967 | 0.9996 / 0.8119 / 0.8960\* |

通过校准后的 F1 最优：Q10 BISCUIT；Q20 bsbit sensitive；Q30 bsbit default；Q40 bsbit sensitive。

## 口径与复现信息

- 当前 bsbit：`bounded-structural-alignment-v1`、`structural-origin-evidence-v1`，最大 edit distance 5，默认 adapter/soft-clip policy，tie-break seed 0，minimal output contract。
- 当前 bsbit 使用 fast index（SA stride 8）、8 个 mapping workers、2 个 BAM compression workers、batch size 16,384、queue batches 2、compression level 1；CPU 固定到 `0,2,4,6,8,10,12,14,16,18`。
- 性能运行二进制和 fast index 的 SHA-256 见页首 Run identity；当前 release-candidate 源码保留相同的 alignment/MAPQ 决策，版本和历史整理不改变本表的 alignment 决策。
- CPU：Intel Core i7-14700K，Linux/WSL；所有正式 performance cells 顺序执行。Wall 是端到端 elapsed time；CPU time 是整个子进程树累计 user + system time；RSS 是 20 ms 采样下同一时刻存活进程树 `VmRSS` 之和的峰值。
- 性能输入将同一 200,000-fragment corpus 重复 25 次得到 5M，只用于吞吐与大差异回归。准确率使用原始 200,000 个唯一 QNAME truth records，避免重复 QNAME 污染配对评估。
- Exact、`+/-5 bp` 和 Q recall 的分母始终是全部 truth。SE 要求 contig、strand 与 unclipped 5′ origin 正确；PE 使用 proper-pair view，两个 mate 都必须满足。PE MAPQ 为两端 MAPQ 的最小值。
- Output reads 是 accuracy corpus 中实际出现的 primary records / 期望 records；Mapped 在 SE 为 mapped primary reads，在 PE 为 proper pairs。bsbit 默认输出完整主记录，但不会强行为没有已验证 placement 的 read 分配坐标。
- 外部工具行来自同机冻结矩阵，采用各工具 native command 与索引；HISAT-3N 使用 sparse output，Bismark 使用原生 sparse output，BISCUIT 与 BSBolt 输出 unmapped records，BitMapperBS 仅支持 directional。BSBolt 版本为 1.6.0。
- 冻结外部版本：BISCUIT 1.10.2；Bismark Rust suite 3.1.0；BitMapperBS 1.0.2.3（commit `699a931`）；HISAT-3N 2.2.1-3n-0.0.3；BSBolt 1.6.0（source `ea4870e`）；samtools/HTSlib 1.22.1。
- 外部工具本轮未重跑，因此非常小的 wall-time 差异不应被解释为统计显著。类似 PE directional Q10（bsbit sensitive 与 BISCUIT 的 F1 仅差约 `1.8e-6`）以及 PE non-directional Q10（约 `2.2e-5`）的差距，也应在独立 references/seeds 上复核，而不应作为确定的普遍优势。

旧外部矩阵与机器可读数据保存在 `workspace/baselines/bsbit-0.1.0/RESULTS.md` 和 `workspace/baselines/bsbit-0.1.0/results/complete-matrix.tsv`。本次 bsbit 原始 runtime、BAM 与 evaluator 输出位于 `/tmp/bsbit-alignment-performance-20260906-UdYnCK/`。
