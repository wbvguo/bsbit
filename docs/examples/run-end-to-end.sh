#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
binary_dir=${BSBIT_BIN_DIR:-"$repository_root/target/release"}
output_dir=${1:-$(mktemp -d)}

mkdir -p "$output_dir"
cd "$repository_root"

cp docs/examples/test.fa "$output_dir/reference.fa"
samtools faidx "$output_dir/reference.fa"

"$binary_dir/bsbit" index \
  --reference "$output_dir/reference.fa" \
  --output "$output_dir/reference.bsbit" \
  --threads 2

"$binary_dir/bsbit" align \
  --index "$output_dir/reference.bsbit" \
  --read1 docs/examples/test_R1.fastq \
  --read2 docs/examples/test_R2.fastq \
  --output "$output_dir/sample.bam" \
  --threads 2 \
  --min-template-span 100 \
  --max-template-span 250 \
  --metrics \
  > "$output_dir/alignment.summary.tsv"

samtools sort -n -o "$output_dir/sample.qname.bam" "$output_dir/sample.bam"
samtools fixmate -m "$output_dir/sample.qname.bam" "$output_dir/sample.fixmate.bam"
samtools sort -o "$output_dir/sample.sorted.bam" "$output_dir/sample.fixmate.bam"
samtools markdup "$output_dir/sample.sorted.bam" "$output_dir/sample.prep.bam"
samtools index "$output_dir/sample.prep.bam"
samtools quickcheck -v "$output_dir/sample.prep.bam"

"$binary_dir/bsbit" call joint \
  --input "$output_dir/sample.prep.bam" \
  --reference "$output_dir/reference.fa" \
  --prefix "$output_dir/sample" \
  --meth-format bed \
  --sample-name demo \
  --min-depth 1 \
  --threads 2

"$binary_dir/bsbit" combine \
  --input "$output_dir/sample.bed" \
  --sample-name demo \
  --output "$output_dir/cohort.bed" \
  --matrix both \
  --min-count 1 \
  --min-prop 1 \
  --threads 2

awk -F '\t' 'NR == 2 && $1 == "bsbit-alignment-metrics-paired-end-v1" && $2 == 4 && $3 == 4 && $6 == 8 { passed = 1 } END { exit !passed }' \
  "$output_dir/alignment.summary.tsv"
test "$(samtools view "$output_dir/sample.bam" | awk '$5 == 60 { count++ } END { print count + 0 }')" -eq 8
grep -q $'demo\t40\t.\tA\tG\t' "$output_dir/sample.vcf"
test -s "$output_dir/sample.bed"
test -s "$output_dir/cohort.level.bed"
test -s "$output_dir/cohort.count.bed"

echo "end-to-end smoke test passed: $output_dir"
