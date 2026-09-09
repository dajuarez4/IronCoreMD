#!/usr/bin/env python3
"""Generate Round-21 SCF-stabilization pilots from Round-20's difficult seed."""

from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROUND20 = ROOT.parent / "round20_16atom_dlm8_conv1e3_nraise5_md100_jakar"
SOURCE = ROUND20 / "cases" / "dlm01_spin200001_vel210001"
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
TESTS = (
    ("conv3e-4_beta1e-2", 3.0e-4, 0.010, 800),
    ("conv1e-4_beta1e-2", 1.0e-4, 0.010, 1000),
    ("conv3e-4_beta5e-3", 3.0e-4, 0.005, 1000),
)


def replace_setting(text: str, name: str, value: str) -> str:
    pattern = rf"(?m)^(\s*{re.escape(name)}\s*=\s*)[^,\n]+"
    updated, count = re.subn(pattern, rf"\g<1>{value}", text)
    if count != 1:
        raise RuntimeError(f"Expected one {name}, found {count}")
    return updated


def main() -> None:
    if not SOURCE.is_dir():
        raise FileNotFoundError(SOURCE)
    (ROOT / "cases").mkdir(parents=True, exist_ok=True)
    (ROOT / "logs").mkdir(exist_ok=True)
    (ROOT / "pseudo").mkdir(exist_ok=True)
    shutil.copy2(ROUND20 / "pseudo" / PSEUDO, ROOT / "pseudo" / PSEUDO)
    rows = []
    for index, (case, threshold, beta, maxstep) in enumerate(TESTS):
        folder = ROOT / "cases" / case
        folder.mkdir(exist_ok=True)
        if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
            raise RuntimeError(f"Refusing to overwrite results in {folder}")
        prefix = f"Fe_bcc16_r21_{index:02d}"
        for stem in ("stage1_4e-3", "stage2_1e-3", "stage3_3e-4", "md_100steps"):
            text = (SOURCE / f"{stem}.in").read_text()
            text = replace_setting(text, "prefix", f"'{prefix}'")
            if stem == "stage3_3e-4":
                text = replace_setting(text, "conv_thr", f"{threshold:.1e}")
                text = replace_setting(text, "electron_maxstep", str(maxstep))
                text = replace_setting(text, "mixing_beta", f"{beta:.3f}d0")
            elif stem == "md_100steps":
                text = replace_setting(text, "conv_thr", f"{threshold:.1e}")
                text = replace_setting(text, "electron_maxstep", str(maxstep))
                text = replace_setting(text, "mixing_beta", f"{beta:.3f}d0")
            (folder / f"{stem}.in").write_text(text)
        metadata = {
            "serial_index": index, "case": case, "source_round": 20,
            "source_case": SOURCE.name, "spin_seed": 200001, "velocity_seed": 210001,
            "md_steps": 100, "md_conv_thr_Ry": threshold,
            "mixing_mode": "local-TF", "mixing_beta": beta, "mixing_ndim": 20,
            "electron_maxstep": maxstep, "nraise": 5, "mpi_tasks": 16,
        }
        (folder / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        rows.append(metadata)
    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (ROOT / "cases.txt").write_text("\n".join(row["case"] for row in rows) + "\n")
    print(f"Generated {len(rows)} Round-21 stabilization pilots")


if __name__ == "__main__":
    main()
