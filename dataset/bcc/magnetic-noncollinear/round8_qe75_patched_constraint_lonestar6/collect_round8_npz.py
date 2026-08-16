#!/usr/bin/env python3
"""Create a safe lightweight NPZ for patched-QE Round-8 SCF validation."""

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
STAGES = ("stage1_4e-3", "stage2_1e-3", "stage3_3e-4")
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
RE_ACC = re.compile(r"estimated scf accuracy\s*<\s*(" + FLOAT + r")\s*Ry", re.I)
RE_ENERGY = re.compile(r"!\s*total energy\s*=\s*(" + FLOAT + r")\s*Ry", re.I)
RE_LOCAL = re.compile(r"^\s*magnetization\s*:\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")", re.I | re.M)
RE_CONS = re.compile(r"^\s*constrained moment\s*:\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")", re.I | re.M)
RE_TOTAL = re.compile(r"total magnetization\s*=\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")", re.I)
RE_ABS = re.compile(r"absolute magnetization\s*=\s*(" + FLOAT + r")", re.I)
RE_FORCE = re.compile(r"Total force\s*=\s*(" + FLOAT + r")", re.I)
RE_PRESS = re.compile(r"P=\s*(" + FLOAT + r")", re.I)
RE_WALL = re.compile(r"PWSCF\s*:.*?CPU\s+(.+?)\s+WALL", re.I)
SUPPORT = ("README.md", "qe75_atomic_constraint_zv.patch", "prepare_qe75_source.sh",
           "build_qe75_patched_lonestar6.sbatch", "generate_round8.py",
           "run_round8_serial_lonestar6.sbatch", "check_round8.sh", "plot_round8.sh",
           "collect_round8_npz.py", "case_manifest.csv", "cases.txt", "round8_serial_status.tsv")


def qfloat(value):
    return float(value.replace("D", "E").replace("d", "e"))


def uarray(values):
    values = [str(value) for value in values]
    return np.asarray(values, dtype=f"U{max([len(x) for x in values] or [1])}")


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def vals(regex, text, group=1):
    return [qfloat(match.group(group)) for match in regex.finditer(text)]


def ragged(payload, name, rows, dtype):
    lengths = np.asarray([len(row) for row in rows], dtype=np.int32)
    offsets = np.empty(len(rows) + 1, dtype=np.int64); offsets[0] = 0
    np.cumsum(lengths, out=offsets[1:])
    payload[name + "_flat"] = np.asarray([x for row in rows for x in row], dtype=dtype)
    payload[name + "_offsets"] = offsets


def wall_seconds(text):
    found = RE_WALL.findall(text)
    if not found:
        return math.nan
    factors = {"d": 86400, "h": 3600, "m": 60, "s": 1}
    return sum(float(n) * factors[u] for n, u in re.findall(r"([0-9.]+)\s*([dhms])", found[-1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=ROOT / "round8_results_light.npz")
    parser.add_argument("--require-terminal", action="store_true")
    args = parser.parse_args()
    with (ROOT / "case_manifest.csv").open(newline="") as handle:
        cases = list(csv.DictReader(handle))

    source_case=[]; source_stem=[]; status=[]; iterations=[]; last_accuracy=[]
    input_path=[]; input_text=[]; input_hash=[]; output_path=[]; output_hash=[]; output_bytes=[]
    accuracy=[]; energy=[]; local_x=[]; local_y=[]; local_z=[]; cons_x=[]; cons_y=[]; cons_z=[]
    total_x=[]; total_y=[]; total_z=[]; absolute=[]; pressure=[]; force=[]
    final_local=[]; final_constrained=[]; walls=[]; cbands=[]; cholesky=[]; errors=[]; excerpts=[]
    for case_index, case in enumerate(cases):
        for stem in STAGES:
            folder = ROOT / "cases" / case["case"]
            inp, out = folder / f"{stem}.in", folder / f"{stem}.out"
            text = out.read_text(errors="replace") if out.is_file() else ""
            if not out.is_file(): state = "NOT_STARTED"
            elif "convergence has been achieved" in text: state = "CONVERGED"
            elif "Error in routine" in text or "convergence NOT achieved" in text: state = "FAILED"
            else: state = "INCOMPLETE"
            mags = [[qfloat(m.group(i)) for i in range(1,4)] for m in RE_LOCAL.finditer(text)]
            constrained = [[qfloat(m.group(i)) for i in range(1,4)] for m in RE_CONS.finditer(text)]
            final = mags[-8:] if len(mags) >= 8 else []
            final_cons = constrained[-8:] if len(constrained) >= 8 else []
            source_case.append(case_index); source_stem.append(stem); status.append(state)
            iterations.append(len(re.findall(r"iteration\s*#", text, re.I)))
            av = vals(RE_ACC, text); last_accuracy.append(av[-1] if av else math.nan)
            input_path.append(str(inp.relative_to(ROOT))); input_text.append(inp.read_text()); input_hash.append(sha(inp))
            output_path.append(str(out.relative_to(ROOT))); output_hash.append(sha(out) if out.is_file() else "")
            output_bytes.append(out.stat().st_size if out.is_file() else 0)
            accuracy.append(av); energy.append(vals(RE_ENERGY, text))
            local_x.append([v[0] for v in mags]); local_y.append([v[1] for v in mags]); local_z.append([v[2] for v in mags])
            cons_x.append([v[0] for v in constrained]); cons_y.append([v[1] for v in constrained]); cons_z.append([v[2] for v in constrained])
            total_x.append(vals(RE_TOTAL, text, 1)); total_y.append(vals(RE_TOTAL, text, 2)); total_z.append(vals(RE_TOTAL, text, 3))
            absolute.append(vals(RE_ABS, text)); pressure.append(vals(RE_PRESS, text)); force.append(vals(RE_FORCE, text))
            final_local.append(sum(math.sqrt(sum(x*x for x in v)) for v in final)/8 if len(final)==8 else math.nan)
            final_constrained.append(sum(math.sqrt(sum(x*x for x in v)) for v in final_cons)/8 if len(final_cons)==8 else math.nan)
            walls.append(wall_seconds(text)); cbands.append(text.lower().count("c_bands:")); cholesky.append(text.lower().count("cholesky")); errors.append(text.count("Error in routine"))
            selected=[line for line in text.splitlines() if any(x in line for x in ("convergence has been achieved","convergence NOT achieved","Error in routine","JOB DONE"))]
            excerpts.append("\n".join(selected + (["--- tail ---"] + text.splitlines()[-30:] if text else [])))

    payload={
        "schema_version":uarray(["ironcoremd.qe.bcc8_patched_qe75_atomic_scf_round8.v1"]),
        "created_utc":uarray([datetime.now(timezone.utc).isoformat()]),
        "case_count":np.asarray([len(cases)],dtype=np.int32), "case_name":uarray([x["case"] for x in cases]),
        "lambda_stage1_Ry":np.asarray([float(x["lambda_stage1"]) for x in cases]),
        "lambda_stage2_Ry":np.asarray([float(x["lambda_stage2"]) for x in cases]),
        "lambda_stage3_Ry":np.asarray([float(x["lambda_stage3"]) for x in cases]),
        "target_local_moment_Bohr":np.asarray([float(x["intended_target_Bohr"]) for x in cases]),
        "source_case_index":np.asarray(source_case,dtype=np.int16), "source_stem":uarray(source_stem),
        "status":uarray(status), "iteration_count":np.asarray(iterations,dtype=np.int32),
        "last_accuracy_Ry":np.asarray(last_accuracy), "final_mean_local_moment_Bohr":np.asarray(final_local),
        "printed_constrained_moment_Bohr":np.asarray(final_constrained),
        "wall_seconds":np.asarray(walls), "c_bands_warning_count":np.asarray(cbands,dtype=np.int32),
        "cholesky_count":np.asarray(cholesky,dtype=np.int32), "qe_error_count":np.asarray(errors,dtype=np.int32),
        "input_relative_path":uarray(input_path), "input_text":uarray(input_text), "input_sha256":uarray(input_hash),
        "output_relative_path":uarray(output_path), "output_sha256":uarray(output_hash),
        "output_bytes":np.asarray(output_bytes,dtype=np.int64), "diagnostic_excerpt":uarray(excerpts),
    }
    for name, rows in (("accuracy_Ry",accuracy),("energy_Ry",energy),("local_mag_x_Bohr",local_x),
                       ("local_mag_y_Bohr",local_y),("local_mag_z_Bohr",local_z),
                       ("constrained_x_Bohr",cons_x),("constrained_y_Bohr",cons_y),("constrained_z_Bohr",cons_z),
                       ("total_mag_x_Bohr",total_x),("total_mag_y_Bohr",total_y),("total_mag_z_Bohr",total_z),
                       ("absolute_magnetization_Bohr",absolute),("pressure_kbar",pressure),("total_force_Ry_Bohr",force)):
        ragged(payload,name,rows,np.float64)
    support=[ROOT/name for name in SUPPORT if (ROOT/name).is_file()]
    payload["support_relative_path"]=uarray([str(x.relative_to(ROOT)) for x in support])
    payload["support_text"]=uarray([x.read_text(errors="replace") for x in support])
    payload["support_sha256"]=uarray([sha(x) for x in support])
    logs=sorted((ROOT/"logs").glob("*")); logs=[x for x in logs if x.is_file() and x.stat().st_size <= 512*1024]
    payload["slurm_log_relative_path"]=uarray([str(x.relative_to(ROOT)) for x in logs])
    payload["slurm_log_text"]=uarray([x.read_text(errors="replace") for x in logs])
    terminal=sum(1 for i in range(len(cases)) if any(status[3*i+j]=="FAILED" for j in range(3)) or all(status[3*i+j]=="CONVERGED" for j in range(3)))
    if args.require_terminal and terminal != len(cases):
        raise SystemExit(f"Refusing incomplete archive: terminal={terminal}/{len(cases)}")
    args.output.parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(args.output,**payload)
    print(args.output.resolve()); print(f"cases={len(cases)} terminal={terminal} arrays={len(payload)}")


if __name__ == "__main__":
    main()
