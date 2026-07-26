#!/usr/bin/env bash

set -euo pipefail

: "${QE_ARCHIVE:?Set QE_ARCHIVE to the QE source archive path}"

QE_INSTALL_ROOT=${QE_INSTALL_ROOT:-"${SCRATCH:-$HOME}/software"}
BUILD_JOBS=${BUILD_JOBS:-8}

QE_ARCHIVE=$(cd "$(dirname "$QE_ARCHIVE")" && pwd)/$(basename "$QE_ARCHIVE")

if [[ ! -f "$QE_ARCHIVE" ]]; then
    echo "QE archive not found: $QE_ARCHIVE" >&2
    exit 1
fi

mkdir -p "$QE_INSTALL_ROOT"
cd "$QE_INSTALL_ROOT"

archive_listing=$(tar -tf "$QE_ARCHIVE")
first_entry=${archive_listing%%$'\n'*}
source_dir=${first_entry%%/*}

if [[ -z "$source_dir" ]]; then
    echo "Could not determine the source directory in $QE_ARCHIVE" >&2
    exit 1
fi

if [[ ! -d "$source_dir" ]]; then
    tar -xf "$QE_ARCHIVE"
fi

cd "$source_dir"
mkdir -p build
cd build

../configure
make -j "$BUILD_JOBS" all

QE_BIN="$PWD/bin/pw.x"
if [[ ! -x "$QE_BIN" ]]; then
    echo "Build completed, but pw.x was not found at $QE_BIN" >&2
    exit 1
fi

echo "Quantum ESPRESSO executable: $QE_BIN"
"$QE_BIN" -help >/dev/null
echo "pw.x startup check passed. Run the QE test suite before production use."
