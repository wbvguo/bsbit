#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
binary_dir=${BSBIT_BIN_DIR:-"$repository_root/target/release"}
output_dir=${1:-$(mktemp -d)}

mkdir -p "$output_dir"
cd "$repository_root"

cp docs/examples/quickstart-reference.fa "$output_dir/reference.fa"

"$binary_dir/bsbit" index \
  -r "$output_dir/reference.fa" \
  -o "$output_dir/reference.bsbit" \
  -t 2

"$binary_dir/bsbit" align \
  --index "$output_dir/reference.bsbit" \
  --read1 docs/examples/quickstart_R1.fastq \
  --read2 docs/examples/quickstart_R2.fastq \
  --output "$output_dir/alignment.bam" \
  --threads 2 \
  --compression-threads 1 \
  --min-template-span 100 \
  --max-template-span 250 \
  --metrics \
  > "$output_dir/alignment.summary.tsv"

samtools sort -n -o "$output_dir/alignment.name.bam" "$output_dir/alignment.bam"
samtools fixmate -m "$output_dir/alignment.name.bam" "$output_dir/alignment.fixmate.bam"
samtools sort -o "$output_dir/alignment.position.bam" "$output_dir/alignment.fixmate.bam"
samtools markdup "$output_dir/alignment.position.bam" "$output_dir/alignment.analysis.bam"
samtools index "$output_dir/alignment.analysis.bam"
samtools quickcheck -v "$output_dir/alignment.analysis.bam"

"$binary_dir/bsbit" call joint \
  --input "$output_dir/alignment.analysis.bam" \
  -r "$output_dir/reference.fa" \
  --meth "$output_dir/methylation.bed" \
  --meth-format bed \
  --vcf "$output_dir/variants.vcf" \
  --sample-name demo \
  --threads 2

"$binary_dir/bsbit" combine \
  --input "$output_dir/methylation.bed" \
  --sample-name demo \
  --output "$output_dir/cohort.bed" \
  --matrix both \
  --min-count 1 \
  --min-prop 1 \
  --threads 2

awk -F '\t' 'NR == 2 && $1 == "bsbit-alignment-metrics-v1" && $2 == 4 && $3 == 4 && $6 == 8 { passed = 1 } END { exit !passed }' \
  "$output_dir/alignment.summary.tsv"
test "$(samtools view "$output_dir/alignment.bam" | awk '$5 == 60 { count++ } END { print count + 0 }')" -eq 8
grep -q $'demo\t40\t.\tA\tG\t' "$output_dir/variants.vcf"
test -s "$output_dir/methylation.bed"
test -s "$output_dir/cohort.level.bed"
test -s "$output_dir/cohort.count.bed"
test ! -e "$output_dir/reference.fa.fai"

echo "end-to-end smoke test passed: $output_dir"
