#!/usr/bin/env bash
# Usage: bash prepare_qe75_source.sh /absolute/path/to/q-e-qe-7.5.tar.gz
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
TARBALL=${1:?Provide the QE 7.5 tar.gz path}
SOURCE_PARENT=${2:-${SCRATCH:?SCRATCH is not defined}/src}
TARGET="$SOURCE_PARENT/q-e-qe-7.5-atomic-fixed"
[[ -s "$TARBALL" ]] || { echo "Missing tarball: $TARBALL" >&2; exit 2; }
[[ ! -e "$TARGET" ]] || { echo "$TARGET already exists; refusing to mix source attempts." >&2; exit 3; }
mkdir -p "$SOURCE_PARENT"
work=$(mktemp -d "$SOURCE_PARENT/qe75_unpack.XXXXXX")
trap 'rm -rf "$work"' EXIT
tar -xzf "$TARBALL" -C "$work"
source_dir=$(find "$work" -mindepth 1 -maxdepth 1 -type d | head -1)
[[ -n "$source_dir" && -f "$source_dir/PW/src/input.f90" ]] || { echo "Unexpected QE archive layout." >&2; exit 4; }
mv "$source_dir" "$TARGET"
patch -d "$TARGET" -p1 < "$ROOT/qe75_atomic_constraint_zv.patch"
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$TARGET/PW/src/input.f90" || {
    echo "Patch verification failed." >&2; exit 5;
}
printf 'QE 7.5 atomic-constraint valence scaling patch\nPrepared: %s\nPatch: %s\n' "$(date -Is)" "$ROOT/qe75_atomic_constraint_zv.patch" > "$TARGET/ATOMIC_CONSTRAINT_PATCH_INFO.txt"
echo "Prepared patched source: $TARGET"
