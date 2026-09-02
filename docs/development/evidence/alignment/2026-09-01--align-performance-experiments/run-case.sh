#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

if (( $# != 6 )); then
  die "usage: ${0##*/} LABEL BINARY single|paired MAPPING_THREADS BGZF_THREADS CPU_LIST"
fi

readonly label="$1"
readonly binary="$(realpath -e -- "$2")"
readonly layout="$3"
readonly mapping_threads="$4"
readonly bgzf_threads="$5"
readonly cpu_list="$6"
readonly run_root="${BSBIT_PERF_RUN_ROOT:-/tmp/bsbit-align-perf-experiments-20260901/runs}"
readonly fixture_root="${BSBIT_PERF_FIXTURE_ROOT:-/tmp/bsbit-current-benchmark-20260831}"
readonly index_path="${BSBIT_PERF_INDEX:-/tmp/bsbit-flattened-20260831/indices/bsbit/current.bsbit}"
readonly read1="${BSBIT_PERF_READ1:-${fixture_root}/inputs/R1.fastq.gz}"
readonly read2="${BSBIT_PERF_READ2:-${fixture_root}/inputs/R2.fastq.gz}"
readonly case_dir="${run_root}/${label}"
readonly output_bam="${case_dir}/output.bam"

[[ ${layout} == single || ${layout} == paired ]] || die "layout must be single or paired"
[[ ${mapping_threads} =~ ^[1-9][0-9]*$ ]] || die "invalid mapping thread count"
[[ ${bgzf_threads} =~ ^[1-9][0-9]*$ ]] || die "invalid BGZF thread count"
[[ ! -e ${case_dir} ]] || die "refusing to overwrite ${case_dir}"
mkdir -p -- "${case_dir}"

{
  printf 'captured_utc=%s\n' "$(date -u +%FT%TZ)"
  printf 'label=%s\nlayout=%s\n' "${label}" "${layout}"
  printf 'binary=%s\n' "${binary}"
  printf 'binary_sha256='; sha256sum -- "${binary}" | cut -d' ' -f1
  printf 'index=%s\nread1=%s\nread2=%s\n' "${index_path}" "${read1}" "${read2}"
  printf 'mapping_threads=%s\nbgzf_threads=%s\ncpu_list=%s\n' \
    "${mapping_threads}" "${bgzf_threads}" "${cpu_list}"
  uptime
} > "${case_dir}/environment.txt"

command=(
  "${binary}" align
  --index "${index_path}"
  --read1 "${read1}"
  --output-bam "${output_bam}"
  --threads "${mapping_threads}"
  --bam-threads "${bgzf_threads}"
)
if [[ ${layout} == paired ]]; then
  command+=(--read2 "${read2}" --metrics)
fi
printf '%q ' taskset -c "${cpu_list}" "${command[@]}" > "${case_dir}/command.txt"
printf '\n' >> "${case_dir}/command.txt"

/usr/bin/time -v -o "${case_dir}/time.txt" \
  taskset -c "${cpu_list}" "${command[@]}" \
  > "${case_dir}/stdout.txt" 2> "${case_dir}/stderr.txt"

sha256sum -- "${output_bam}" > "${case_dir}/output.sha256"
samtools quickcheck -v "${output_bam}" > "${case_dir}/quickcheck.txt"
samtools view -c "${output_bam}" > "${case_dir}/record-count.txt"
samtools view "${output_bam}" \
  | awk 'BEGIN { OFS="\t" }
         { total++; if (and($2,4)) unmapped++; else mapped++;
           if ($6 == "150M") exact_150m++; }
         END { print "total",total; print "mapped",mapped;
               print "unmapped",unmapped; print "cigar_150M",exact_150m }' \
  > "${case_dir}/bam-summary.tsv"

awk -F: '
  /User time \(seconds\)/ { gsub(/^ +/,"",$2); user=$2 }
  /System time \(seconds\)/ { gsub(/^ +/,"",$2); sys=$2 }
  /Elapsed \(wall clock\) time/ {
    value=$2; sub(/^.*: /,"",$0); value=$0;
    pieces=split(value,a,":");
    if (pieces == 3) wall=a[1]*3600+a[2]*60+a[3];
    else wall=a[1]*60+a[2]
  }
  /Maximum resident set size/ { gsub(/^ +/,"",$2); rss=$2 }
  END { print "wall_s\tuser_s\tsys_s\tpeak_rss_kib";
        printf "%.6f\t%s\t%s\t%s\n",wall,user,sys,rss }
' "${case_dir}/time.txt" > "${case_dir}/summary.tsv"

printf 'completed %s\n' "${case_dir}"
