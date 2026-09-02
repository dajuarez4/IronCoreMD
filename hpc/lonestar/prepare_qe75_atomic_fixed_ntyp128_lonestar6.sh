#!/usr/bin/env bash
# Prepare a clean QE 7.5 tree with the IronCoreMD atomic-constraint and ntypx patches.
# Usage:
#   bash prepare_qe75_atomic_fixed_ntyp128_lonestar6.sh /path/q-e-qe-7.5.tar.gz [source_parent]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
TARBALL=${1:?Provide the official QE 7.5 tar.gz path}
SOURCE_PARENT=${2:-${SCRATCH:?SCRATCH is not defined}/src}
TARGET=${TARGET_DIR:-$SOURCE_PARENT/q-e-qe-7.5-atomic-fixed-ntyp128}
PATCH_FILE="$SCRIPT_DIR/qe75_atomic_constraint_ntyp128.patch"
INFO_FILE=ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt

[[ -s "$TARBALL" ]] || { echo "Missing archive: $TARBALL" >&2; exit 2; }
[[ -s "$PATCH_FILE" ]] || { echo "Missing patch: $PATCH_FILE" >&2; exit 2; }
[[ ! -e "$TARGET" ]] || {
    echo "$TARGET already exists; refusing to mix source/build attempts." >&2
    exit 3
}

mkdir -p "$SOURCE_PARENT"
WORK=$(mktemp -d "$SOURCE_PARENT/qe75_ntyp128_unpack.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$TARBALL" -C "$WORK"

mapfile -t SOURCE_DIRS < <(find "$WORK" -mindepth 1 -maxdepth 1 -type d -print)
[[ ${#SOURCE_DIRS[@]} -eq 1 ]] || {
    echo "Expected one top-level source directory in the QE archive." >&2
    exit 4
}
SOURCE_DIR=${SOURCE_DIRS[0]}
[[ -f "$SOURCE_DIR/Modules/parameters.f90" && -f "$SOURCE_DIR/PW/src/input.f90" ]] || {
    echo "Unexpected QE archive layout." >&2
    exit 4
}

patch -d "$SOURCE_DIR" -p1 < "$PATCH_FILE"
grep -Eq 'ntypx[[:space:]]*=[[:space:]]*128' "$SOURCE_DIR/Modules/parameters.f90" || {
    echo "ntypx=128 verification failed." >&2
    exit 5
}
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$SOURCE_DIR/PW/src/input.f90" || {
    echo "Atomic-constraint patch verification failed." >&2
    exit 5
}

{
    echo "Quantum ESPRESSO 7.5 private IronCoreMD build for Lonestar6"
    echo "Patch 1: Modules/parameters.f90 ntypx increased from 10 to 128"
    echo "Patch 2: PW/src/input.f90 atomic mcons target multiplied by zv(nt)"
    echo "Prepared UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Source archive: $TARBALL"
    sha256sum "$TARBALL"
    echo "Patch file: $PATCH_FILE"
    sha256sum "$PATCH_FILE"
} > "$SOURCE_DIR/$INFO_FILE"

mv "$SOURCE_DIR" "$TARGET"
echo "Prepared patched source: $TARGET"
echo "Next: sbatch build_qe75_atomic_fixed_ntyp128_lonestar6.sbatch"
