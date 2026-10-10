#!/bin/bash
# gale Breeze two-lane / W2.1 SIMD gate runner (no sbt): reproduces the build.sbt
# aliases breezeLaneA / breezeLaneB with `java -cp <exported Jmh/fullClasspath>`.
# Usage: run.sh <mode> <outdir> [extra JMH args...]
#   mode: env | smokeA | smokeB | smokeGate | laneA | laneB | gate | vtest
set -euo pipefail
if (( $# < 2 )); then echo "usage: $0 <mode> <outdir> [extra JMH args...]" >&2; exit 2; fi
MODE=$1; OUT=$2; shift 2
case "$MODE" in env|smokeA|smokeB|smokeGate|laneA|laneB|gate|vtest) ;; *) echo "unknown mode $MODE" >&2; exit 2 ;; esac
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
  smokeA|smokeB|smokeGate|laneA|laneB) reset_sidecars "$OUT"; rm -f "$OUT/result.json" ;;
esac
case "$MODE" in
  env) ;;
  smokeA) "${HOST[@]}" "${LANE_A[@]}" -f 1 -wi 1 -i 1 -w 300ms -r 300ms -p n=65536 -rff "$OUT/result.json" "$@" 'gale.bench.BlasL1BreezeJmh.*[Dd]ot$' ;;
  smokeB) "${HOST[@]}" "${LANE_B[@]}" -f 1 -wi 1 -i 1 -w 300ms -r 300ms -p n=256 -rff "$OUT/result.json" "$@" 'gale.bench.BlasL2BreezeJmh.*Gemv$' 'gale.bench.BlasL1BreezeJmh.*[Dd]ot$' ;;
  smokeGate) "${HOST[@]}" -jvmArgsAppend "-Dgale.bench.lane=B -Dgale.bench.netlibSidecar=$OUT/netlib.jsonl" -f 1 -wi 1 -i 1 -w 300ms -r 300ms -p n=4096 -prof gc -rf json -rff "$OUT/result.json" "$@" 'gale.bench.SimdL1SpikeJmh.*[Dd]ot' ;;
  laneA) "${HOST[@]}" "${LANE_A[@]}" -rff "$OUT/result.json" "$@" '.*BreezeJmh.*' ;;
  laneB) "${HOST[@]}" "${LANE_B[@]}" -rff "$OUT/result.json" "$@" '.*BreezeJmh.*' ;;
  gate)
    # docs/verification/w21-simd-spike/README.md recipe, AVX-512 (no UseAVX cap).
    mkdir -p "$OUT/ops" "$OUT/dispatch"
    reset_sidecars "$OUT/ops"; reset_sidecars "$OUT/dispatch"
    rm -f "$OUT/ops/ops-laneB.json" "$OUT/dispatch/dispatch.json"
    "${HOST[@]}" -jvmArgsAppend "-Dgale.bench.lane=B -Dgale.bench.netlibSidecar=$OUT/ops/netlib.jsonl" -f 3 -wi 5 -i 10 -w 1s -r 1s -rf json -rff "$OUT/ops/ops-laneB.json" gale.bench.SimdL1SpikeJmh
    "${HOST[@]}" -jvmArgsAppend "-Dgale.bench.lane=B -Dgale.bench.netlibSidecar=$OUT/dispatch/netlib.jsonl" -f 3 -wi 5 -i 10 -w 1s -r 1s -prof gc -rf json -rff "$OUT/dispatch/dispatch.json" gale.bench.BackendDispatchJmh
    ;;
  vtest)
    # Equivalent of `sbt vectorBackendTest` (forked, Test/javaOptions --add-modules).
    "${PIN[@]}" "$JAVA" --add-modules=jdk.incubator.vector -cp "$(cat vtest.cp)" org.junit.runner.JUnitCore \
      gale.backend.jvm.vector.VectorL1KernelsSuite gale.backend.jvm.vector.VectorGemmSuite gale.backend.jvm.vector.VectorBackendConformanceSuite
    ;;
  *) echo "unknown mode $MODE"; exit 2 ;;
esac
