#!/usr/bin/env python3
"""Extract converged exploratory BCC MD cycles to CSV."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re


OUTPUT = Path("md_exploratory_4000K_10steps.out")
CSV = Path("md_exploratory_4000K_10steps.csv")


def last(pattern: str, text: str) -> str | None:
    values = re.findall(pattern, text, flags=re.MULTILINE)
    return values[-1] if values else None


def value(text: str | None) -> float | None:
    return float(text) if text is not None else None


def main() -> None:
    if not OUTPUT.exists():
        raise SystemExit(f"Missing {OUTPUT}")

    text = OUTPUT.read_text(encoding="utf-8", errors="replace")
    convergence = list(re.finditer(r"convergence has been achieved in\s+(\d+) iterations", text))
    rows: list[dict[str, object]] = []

    for cycle, match in enumerate(convergence):
        start = convergence[cycle - 1].end() if cycle else 0
        block = text[start:match.end()]
        post_end = convergence[cycle + 1].start() if cycle + 1 < len(convergence) else len(text)
        post_block = text[match.end():post_end]
        vectors = re.findall(
            r"total magnetization\s*=\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            block,
        )
        mx, my, mz = map(float, vectors[-1]) if vectors else (None, None, None)
        force = value(last(r"Total force\s*=\s*([-+0-9.Ee]+)", post_block))
        correction = value(last(r"Total SCF correction\s*=\s*([-+0-9.Ee]+)", post_block))
        rows.append(
            {
                "electronic_cycle": cycle,
                "scf_iterations": int(match.group(1)),
                "accuracy_Ry": value(last(r"estimated scf accuracy\s*<\s*([-+0-9.Ee]+)", block)),
                "Mx_muB": mx,
                "My_muB": my,
                "Mz_muB": mz,
                "Mnorm_muB": math.sqrt(mx * mx + my * my + mz * mz) if mx is not None else None,
                "Mabs_muB": value(last(r"absolute magnetization\s*=\s*([-+0-9.Ee]+)", block)),
                "total_energy_Ry": value(last(r"^!\s+total energy\s*=\s*([-+0-9.Ee]+)", block)),
                "temperature_K": value(last(r"temperature\s*=\s*([-+0-9.Ee]+)\s+K", block)),
                "pressure_kbar": value(last(r"\bP=\s*([-+0-9.Ee]+)", post_block)),
                "total_force_Ry_per_bohr": force,
                "scf_force_correction_Ry_per_bohr": correction,
                "force_correction_ratio": correction / force if force not in (None, 0.0) and correction is not None else None,
            }
        )

    fields = list(rows[0]) if rows else []
    with CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"converged_electronic_cycles={len(rows)}")
    print(f"job_done={'JOB DONE.' in text}")
    print(f"csv={CSV}")


if __name__ == "__main__":
    main()
