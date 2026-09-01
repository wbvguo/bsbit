#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
binary_dir=${BSBIT_BIN_DIR:-"$repository_root/target/release"}
output_dir=${1:-$(mktemp -d)}

mkdir -p "$output_dir"
cd "$repository_root"

cp docs/examples/quickstart-reference.fa "$output_dir/reference.fa"
samtools faidx "$output_dir/reference.fa"

"$binary_dir/bsbit" index \
  --reference "$output_dir/reference.fa" \
  --output "$output_dir/reference.bsbit" \
  --threads 2

"$binary_dir/bsbit" align \
  --index "$output_dir/reference.bsbit" \
  --read1 docs/examples/quickstart_R1.fastq \
  --read2 docs/examples/quickstart_R2.fastq \
  --output-bam "$output_dir/alignment.bam" \
  --threads 2 \
  --bam-threads 1 \
  --min-template-span 100 \
  --max-template-span 250 \
  --metrics \
  > "$output_dir/alignment.summary.tsv"

samtools quickcheck -v "$output_dir/alignment.bam"
awk -F '\t' 'NR == 2 && $1 == "bsbit-alignment-metrics-v2" && $2 == 4 && $3 == 4 && $6 == 8 { passed = 1 } END { exit !passed }' \
  "$output_dir/alignment.summary.tsv"
test "$(samtools view "$output_dir/alignment.bam" | awk '$5 == 60 { count++ } END { print count + 0 }')" -eq 8

echo "quick-start smoke test passed: $output_dir"
