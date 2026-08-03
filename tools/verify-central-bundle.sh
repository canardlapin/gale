#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "usage: $0 [--require-signatures] BUNDLE.zip" >&2
  exit 2
}

require_signatures=0
if [[ "${1:-}" == "--require-signatures" ]]; then
  require_signatures=1
  shift
fi
[[ $# -eq 1 ]] || usage
bundle=$1
[[ -f "$bundle" ]] || { echo "bundle does not exist: $bundle" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
unzip -q "$bundle" -d "$tmp"

artifacts=(
  gale-core_3
  gale-core_sjs1_3
  gale-laws_3
  gale-laws_sjs1_3
  gale-interop-breeze_3
  gale-backend-jvm-vector_3
  gale-backend-jvm-native_3
  gale-backend-jvm-blas-ffm_3
)

for artifact in "${artifacts[@]}"; do
  poms=$(find "$tmp" -type f -path "*/$artifact/*/$artifact-*.pom" | sort)
  pom_count=$(printf '%s\n' "$poms" | sed '/^$/d' | wc -l | tr -d ' ')
  [[ "$pom_count" -eq 1 ]] || {
    echo "expected exactly one POM for $artifact, found $pom_count" >&2
    exit 1
  }

  pom=$(printf '%s\n' "$poms" | sed -n '1p')
  version=$(basename "$(dirname "$pom")")
  [[ "$version" != *SNAPSHOT* ]] || {
    echo "snapshot version in Central bundle: $artifact:$version" >&2
    exit 1
  }

  directory=$(dirname "$pom")
  for suffix in pom jar sources.jar javadoc.jar; do
    file="$directory/$artifact-$version"
    case "$suffix" in
      pom) file+=".pom" ;;
      jar) file+=".jar" ;;
      sources.jar) file+="-sources.jar" ;;
      javadoc.jar) file+="-javadoc.jar" ;;
    esac
    [[ -f "$file" ]] || { echo "missing $suffix for $artifact:$version" >&2; exit 1; }
    if (( require_signatures )); then
      [[ -f "$file.asc" ]] || { echo "missing signature for $file" >&2; exit 1; }
    fi
  done
done

echo "Central bundle verified: ${#artifacts[@]} admitted artifacts, version $version"
