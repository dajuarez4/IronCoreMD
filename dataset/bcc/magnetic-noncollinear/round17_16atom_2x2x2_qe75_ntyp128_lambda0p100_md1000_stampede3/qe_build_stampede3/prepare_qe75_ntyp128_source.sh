#!/usr/bin/env bash
# Usage: bash prepare_qe75_ntyp128_source.sh /absolute/path/q-e-qe-7.5.tar.gz
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
TARBALL=${1:?Provide the official QE 7.5 tar.gz path}
SOURCE_PARENT=${2:-${SCRATCH:?SCRATCH is not defined}/src}
TARGET="$SOURCE_PARENT/q-e-qe-7.5-atomic-fixed-ntyp128"
PATCH_FILE="$ROOT/qe75_atomic_constraint_ntyp128.patch"

[[ -s "$TARBALL" ]] || { echo "Missing tarball: $TARBALL" >&2; exit 2; }
[[ -s "$PATCH_FILE" ]] || { echo "Missing patch: $PATCH_FILE" >&2; exit 2; }
[[ ! -e "$TARGET" ]] || {
    echo "$TARGET already exists; refusing to mix source/build attempts." >&2
    exit 3
}

mkdir -p "$SOURCE_PARENT"
work=$(mktemp -d "$SOURCE_PARENT/qe75_ntyp128_unpack.XXXXXX")
trap 'rm -rf "$work"' EXIT
tar -xzf "$TARBALL" -C "$work"
source_dir=$(find "$work" -mindepth 1 -maxdepth 1 -type d | head -1)
[[ -n "$source_dir" && -f "$source_dir/Modules/parameters.f90" && -f "$source_dir/PW/src/input.f90" ]] || {
    echo "Unexpected QE archive layout." >&2
    exit 4
}
patch -d "$source_dir" -p1 < "$PATCH_FILE"

grep -Eq 'INTEGER,[[:space:]]*PARAMETER[[:space:]]*::[[:space:]]*ntypx[[:space:]]*=[[:space:]]*128' \
    "$source_dir/Modules/parameters.f90" || { echo "ntypx patch verification failed." >&2; exit 5; }
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$source_dir/PW/src/input.f90" || {
    echo "atomic-constraint patch verification failed." >&2
    exit 5
}
grep -q 'CHARACTER(LEN=6).*atm' "$source_dir/Modules/ions_base.f90" || {
    echo "QE 7.5 six-character species-label check failed." >&2
    exit 5
}

{
    echo "QE 7.5 private Stampede3 build"
    echo "Patch 1: ntypx increased from 10 to 128"
    echo "Patch 2: atomic mcons target scaled by pseudopotential valence zv"
    echo "Prepared UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Source archive: $TARBALL"
    sha256sum "$TARBALL"
    echo "Patch: $PATCH_FILE"
    sha256sum "$PATCH_FILE"
} > "$source_dir/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt"

mv "$source_dir" "$TARGET"
echo "Prepared patched source: $TARGET"
