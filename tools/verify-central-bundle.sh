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
)

group_root="$tmp/io/github/canardlapin"
[[ -d "$group_root" ]] || {
  echo "missing io.github.canardlapin Maven group in Central bundle" >&2
  exit 1
}

expected_artifacts=$(printf '%s\n' "${artifacts[@]}" | sort)
actual_artifacts=$(
  find "$group_root" -type f -name '*.pom' -exec sh -c '
    for pom do
      basename "$(dirname "$(dirname "$pom")")"
    done
  ' sh {} + | sort -u
)
[[ "$actual_artifacts" == "$expected_artifacts" ]] || {
  echo "Central bundle artifact set differs from the admitted 0.1 manifest" >&2
  diff -u <(printf '%s\n' "$expected_artifacts") <(printf '%s\n' "$actual_artifacts") >&2 || true
  exit 1
}

bundle_version=""
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
  if [[ -z "$bundle_version" ]]; then
    bundle_version=$version
  elif [[ "$version" != "$bundle_version" ]]; then
    echo "mixed versions in Central bundle: expected $bundle_version, found $artifact:$version" >&2
    exit 1
  fi

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

echo "Central bundle verified: ${#artifacts[@]} admitted artifacts, version $bundle_version"
