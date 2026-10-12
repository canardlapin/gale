#!/bin/bash
# gale Breeze two-lane / W2.1 SIMD gate runner (no sbt): reproduces the build.sbt
# aliases breezeLaneA / breezeLaneB with `java -cp <exported Jmh/fullClasspath>`.
# hc-max2 max hill-climb; a copy of run.sh (same setup, pinning and SHA binding), lane A
# (-jvmArgs replaces the host's --add-modules, so forks have no Vector module).
# Usage: mx2.sh <maxA|profA> <outdir>
#   maxA:  MaxHillClimbJmh candidates, -f 8 -wi 5 -i 15 -w 500ms -r 500ms, n in {1024, 65536, 1048576}.
#   profA: -f 2 at n=1048576 on scanMax|scan4Max|keyMax with -prof perfasm; if that fails,
#          -prof perfnorm; if both fail, nothing is profiled (see prof.txt).
set -euo pipefail
if (( $# < 2 )); then echo "usage: $0 <mode> <outdir> [extra JMH args...]" >&2; exit 2; fi
MODE=$1; OUT=$2; shift 2
case "$MODE" in maxA|profA) ;; *) echo "unknown mode $MODE" >&2; exit 2 ;; esac
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

LANE_A=(-jvmArgs "-Xms4g -Xmx4g -Dgale.bench.sourceSha=$SOURCE_SHA -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=$OUT/netlib.jsonl" -p backend=pure -rf json)
LANE_B=(-jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.sourceSha=$SOURCE_SHA -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=$OUT/netlib.jsonl" -p backend=pure,vector -rf json)
HOST=("${PIN[@]}" "$JAVA" --add-modules=jdk.incubator.vector "-Dgale.bench.sourceSha=$SOURCE_SHA" -cp "$CP" org.openjdk.jmh.Main -foe true)

reset_sidecars() { rm -f "$1/netlib.jsonl" "$1/netlib.jsonl.validation.jsonl"; }

case "$MODE" in
  maxA|profA) reset_sidecars "$OUT"; rm -f "$OUT/result.json" ;;
esac
MX='gale.bench.MaxHillClimbJmh.(curMax|scanMax|scan4Max|keyMax|reduceMax|breezeShapeMax)$'
PR='gale.bench.MaxHillClimbJmh.(scanMax|scan4Max|keyMax)$'
case "$MODE" in
  maxA) "${HOST[@]}" "${LANE_A[@]}" -f 8 -wi 5 -i 15 -w 500ms -r 500ms -p n=1024,65536,1048576 -rff "$OUT/result.json" "$MX" ;;
  profA)
    set +e
    { echo "perf: $(command -v perf || echo missing)"; ls /usr/lib*/hsdis* "$(dirname "$(dirname "$JAVA")")"/lib/hsdis* 2>/dev/null || echo "hsdis: not found"; cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null | sed 's/^/perf_event_paranoid: /'; } > "$OUT/prof.txt" 2>&1
    "${HOST[@]}" "${LANE_A[@]}" -f 2 -p n=1048576 -prof perfasm -rff "$OUT/perfasm.json" "$PR" > "$OUT/perfasm.txt" 2>&1
    rc=$?; echo "perfasm rc=$rc" >> "$OUT/prof.txt"
    if [ $rc -ne 0 ] || [ -n "${FORCE_PERFNORM:-}" ]; then
      "${HOST[@]}" "${LANE_A[@]}" -f 2 -p n=1048576 -prof perfnorm -rff "$OUT/perfnorm.json" "$PR" > "$OUT/perfnorm.txt" 2>&1
      echo "perfnorm rc=$?" >> "$OUT/prof.txt"
    fi
    true ;;
esac
