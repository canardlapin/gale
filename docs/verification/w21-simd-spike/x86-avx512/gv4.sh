#!/bin/bash
# gale Breeze two-lane / W2.1 SIMD gate runner (no sbt): reproduces the build.sbt
# aliases breezeLaneA / breezeLaneB with `java -cp <exported Jmh/fullClasspath>`.
# HC-4 vector gemv; a copy of run.sh (same setup, pinning, SHA binding and sidecars), lane B.
# GALE_BENCH_ROOT selects the bundle.
# Usage: gv4.sh <hill|pairB> <outdir>
#   hill:  GemvHillClimbJmh, -f 3 -wi 5 -i 5 (new bundle only).
#   pairB: BlasL2 gemv/gemvT pairs and LbfgsBreezeJmh logistic, -p backend=vector -p budget=fixed.
set -euo pipefail
if (( $# < 2 )); then echo "usage: $0 <mode> <outdir> [extra JMH args...]" >&2; exit 2; fi
MODE=$1; OUT=$2; shift 2
case "$MODE" in hill|pairB) ;; *) echo "unknown mode $MODE" >&2; exit 2 ;; esac
ROOT=${GALE_BENCH_ROOT:-/scratch/brad/gale-bench}
BUNDLE=$ROOT/bundle
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
trap 'rc=$?; printf "rc=%s\n" "$rc" > "$OUT/rc.txt"' EXIT
module load java/25 >/dev/null 2>&1
# The site sets JAVA_TOOL_OPTIONS=-Xmx2g; JMH forks inherit the environment, so
# drop it: heap comes only from the lane's explicit -Xms4g -Xmx4g.
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
cd "$BUNDLE"
CP=$(cat jmh.cp)
SOURCE_SHA=$(cat SHA)
if [[ ! "$SOURCE_SHA" =~ ^[0-9a-f]{40}$ ]]; then echo "bundle SHA must be a full commit ID" >&2; exit 2; fi
if [[ -z "$CP" ]]; then echo "empty JMH classpath" >&2; exit 2; fi
JAVA=$(command -v java)

# Pin to the L3 domain (CCD) of the first CPU and its NUMA node.
CCD=$(cat /sys/devices/system/cpu/cpu0/cache/index3/shared_cpu_list)
NUMANODE=$(ls -d /sys/devices/system/cpu/cpu0/node* 2>/dev/null | head -1 | sed 's/.*node//' || true)
PIN=(numactl --physcpubind="$CCD" --membind="${NUMANODE:-0}")

{
  echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "job: ${SLURM_JOB_ID:-none}"
  echo "sha: $SOURCE_SHA"; echo "uptime: $(uptime)"
  echo "java: $JAVA"; "$JAVA" -version 2>&1
  echo "pin: ${PIN[*]}"
  echo "env JAVA_TOOL_OPTIONS=${JAVA_TOOL_OPTIONS:-<unset>} LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<unset>}"
  lscpu; numactl -H
  echo "--- blas/lapack visible to the loader:"
  ldconfig -p 2>/dev/null | grep -iE 'libblas|liblapack|openblas|flexiblas' || echo "(none in ldconfig cache)"
  LIB_DIRS=${LD_LIBRARY_PATH:-}
  for d in ${LIB_DIRS//:/ }; do ls "$d" 2>/dev/null | grep -iE '^lib(blas|lapack|openblas|flexiblas)' | sed "s|^|$d/|" || true; done
} > "$OUT/env.txt" 2>&1

HOST=("${PIN[@]}" "$JAVA" --add-modules=jdk.incubator.vector "-Dgale.bench.sourceSha=$SOURCE_SHA" -cp "$CP" org.openjdk.jmh.Main -foe true)

reset_sidecars() { rm -f "$1/netlib.jsonl" "$1/netlib.jsonl.validation.jsonl"; }

case "$MODE" in
  hill|pairB) reset_sidecars "$OUT"; rm -f "$OUT/result.json" ;;
esac
B_ARGS=(-jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.sourceSha=$SOURCE_SHA -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=$OUT/netlib.jsonl" -rf json)
case "$MODE" in
  hill) "${HOST[@]}" "${B_ARGS[@]}" -f 3 -wi 5 -i 5 -rff "$OUT/result.json" "$@" 'gale.bench.GemvHillClimbJmh.*' ;;
  pairB) "${HOST[@]}" "${B_ARGS[@]}" -p backend=vector -p budget=fixed -rff "$OUT/result.json" "$@" \
           'gale.bench.BlasL2BreezeJmh.(galeGemvT|breezeGemvT|galeGemv|breezeGemv)$' 'gale.bench.LbfgsBreezeJmh.*[Ll]ogistic' ;;
esac
