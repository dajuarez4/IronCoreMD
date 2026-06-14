#!/usr/bin/env python3
"""Convenience FCC wrapper for the generic NPZ -> TDEP converter."""

from __future__ import annotations

from npz_to_tdep import main


if __name__ == "__main__":
    main(default_phase="fcc")
