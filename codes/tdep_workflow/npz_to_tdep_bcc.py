#!/usr/bin/env python3
"""Backward-compatible BCC wrapper for the generic NPZ -> TDEP converter."""

from __future__ import annotations

from npz_to_tdep import main, sanitize_stem, write_tdep_folder


if __name__ == "__main__":
    main(default_phase="bcc")
