#!/usr/bin/env python3
"""Build a Round-22 HELD trajectory, phonons, still image, and animated dashboard."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-round22-held-dashboard")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import numpy as np
from scipy.optimize import linear_sum_assignment

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize


HERE = Path(__file__).resolve().parent
IRONCOREMD = HERE.parents[3]
WORKSPACE = IRONCOREMD.parent
HELD_ROOT = WORKSPACE / "HELD"
if str(HELD_ROOT) not in sys.path:
    sys.path.insert(0, str(HELD_ROOT))

from HELD import fit_case, plot_heatmap, plot_mean_dispersion  # noqa: E402


CASE = "conv1e-4_beta1e-2_32atom"
DEFAULT_INPUT = HERE / "cases" / CASE / "md_100steps.in"
DEFAULT_OUTPUT = HERE / "md_100steps.out"
DEFAULT_NPZ = HERE / "round22_held_trajectory.npz"
DEFAULT_COEFFICIENTS = HERE / "round22_held_coefficients.csv"
DEFAULT_STEPS_CSV = HERE / "round22_held_steps.csv"
DEFAULT_DISPERSION_DATA = HERE / "round22_held_mean_dispersion.dat"
DEFAULT_DISPERSION_PNG = HERE / "round22_held_mean_dispersion.png"
DEFAULT_HEATMAP_CACHE = HERE / "round22_held_phonon_frames.npz"
DEFAULT_HEATMAP_PNG = HERE / "round22_held_heatmap.png"
DEFAULT_SUMMARY = HERE / "round22_held_summary.json"
DEFAULT_PNG = HERE / "round22_held_md_dashboard_final.png"
DEFAULT_GIF = HERE / "round22_held_md_dashboard.gif"
NATOMS = 32
REPETITIONS = (4, 2, 2)
REQUESTED_STEPS = 100
REQUESTED_MOMENT = 2.0
DT_AU = 20.67
AU_TIME_SECONDS = 2.4188843265864e-17
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

COLORS = {
    "navy": "#0b2d4d", "orange": "#f47c20", "purple": "#7b4ab5",
    "green": "#2a9d68", "blue": "#2878b5", "red": "#c44e52",
    "gray": "#687386", "light": "#d8dee7",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--previous-output", type=Path, default=HERE / "md_100steps-0.out")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep every Nth frame; the final frame is always included")
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--points-per-segment", type=int, default=60)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--gif", type=Path, default=DEFAULT_GIF)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def scalars(pattern: str, text: str, flags: int = re.I) -> np.ndarray:
    return np.asarray([number(value) for value in re.findall(pattern, text, flags)], dtype=float)


def input_positions_cell_symbols(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    lines = path.read_text().splitlines()
    p0 = lines.index("ATOMIC_POSITIONS crystal") + 1
    position_rows = [row.split() for row in lines[p0:p0 + NATOMS]]
    positions = np.asarray([[number(value) for value in row[1:4]] for row in position_rows])
    symbols = ["Fe"] * NATOMS
    c0 = lines.index("CELL_PARAMETERS angstrom") + 1
    cell = np.asarray([[number(value) for value in row.split()[:3]] for row in lines[c0:c0 + 3]])
    return positions, cell, symbols


def output_position_blocks(text: str) -> np.ndarray:
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    for index, line in enumerate(lines):
        if not line.startswith("ATOMIC_POSITIONS (crystal)"):
            continue
        rows: list[list[float]] = []
        for candidate in lines[index + 1:index + 1 + NATOMS]:
            fields = candidate.split()
            if len(fields) < 4:
                break
            rows.append([number(value) for value in fields[1:4]])
        if len(rows) == NATOMS:
            blocks.append(rows)
    return np.asarray(blocks, dtype=float).reshape(-1, NATOMS, 3)


def force_blocks(text: str) -> np.ndarray:
    """Read only the total force block following each QE force header."""
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    force_pattern = re.compile(
        rf"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        re.I,
    )
    for index, line in enumerate(lines):
        if "Forces acting on atoms (cartesian axes, Ry/au)" not in line:
            continue
        rows: list[list[float]] = []
        cursor = index + 1
        while cursor < len(lines) and len(rows) < NATOMS:
            match = force_pattern.search(lines[cursor])
            if match:
                rows.append([number(value) for value in match.groups()])
            elif rows and lines[cursor].strip():
                break
            cursor += 1
        if len(rows) == NATOMS:
            blocks.append(rows)
    return np.asarray(blocks, dtype=float).reshape(-1, NATOMS, 3)


def vector_frames(text: str, label: str) -> np.ndarray:
    rows = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    values = np.asarray([[number(value) for value in row] for row in rows], dtype=float)
    complete = len(values) // NATOMS
    return values[:complete * NATOMS].reshape(complete, NATOMS, 3)


def converged_cell_history(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vectors: list[list[float]] = []
    absolute: list[float] = []
    iterations: list[int] = []
    current_vector = [np.nan, np.nan, np.nan]
    current_absolute = np.nan
    count = 0
    for line in text.splitlines():
        match = re.search(
            rf"total magnetization\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})", line
        )
        if match:
            current_vector = [number(value) for value in match.groups()]
        match = re.search(rf"absolute magnetization\s*=\s*({FLOAT})", line)
        if match:
            current_absolute = number(match.group(1))
        if "iteration #" in line:
            count += 1
        if "convergence has been achieved" in line:
            vectors.append(current_vector)
            absolute.append(current_absolute)
            iterations.append(count)
            count = 0
    return np.asarray(vectors), np.asarray(absolute), np.asarray(iterations)


def canonical_bcc_fractional(repetitions: tuple[int, int, int]) -> np.ndarray:
    nx, ny, nz = repetitions
    sites = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                sites.append([ix / nx, iy / ny, iz / nz])
                sites.append([(ix + 0.5) / nx, (iy + 0.5) / ny, (iz + 0.5) / nz])
    return np.asarray(sites, dtype=float)


def ideal_sites_in_atom_order(initial: np.ndarray, cell: np.ndarray) -> np.ndarray:
    canonical = canonical_bcc_fractional(REPETITIONS)
    delta = initial[:, None, :] - canonical[None, :, :]
    delta -= np.rint(delta)
    squared_distance = np.sum((delta @ cell) ** 2, axis=2)
    rows, columns = linear_sum_assignment(squared_distance)
    if not np.array_equal(rows, np.arange(NATOMS)):
        raise ValueError("Unexpected atom rows from ideal-site assignment")
    distances = np.sqrt(squared_distance[rows, columns])
    if float(np.max(distances)) > 1.0:
        raise ValueError(f"Ideal-site assignment is unsafe; max distance={np.max(distances):.3f} Å")
    return canonical[columns]


def selected_indices(count: int, stride: int) -> np.ndarray:
    if stride < 1:
        raise ValueError("--stride must be positive")
    indices = np.arange(0, count, stride, dtype=int)
    if indices[-1] != count - 1:
        indices = np.append(indices, count - 1)
    return indices


def source_text(args: argparse.Namespace) -> str:
    previous = args.previous_output.read_text(errors="replace")
    # Drop the interrupted SCF attempt, keeping the last propagated configuration.
    last_dynamics = previous.rfind("Entering Dynamics:")
    interrupted = previous.find("Self-consistent Calculation", last_dynamics)
    if interrupted >= 0:
        previous = previous[:interrupted]
    return previous + "\n" + args.output.read_text(errors="replace")


def prepare_trajectory(args: argparse.Namespace) -> dict[str, np.ndarray]:
    text = source_text(args)
    initial, cell, symbols = input_positions_cell_symbols(args.input)
    propagated = output_position_blocks(text)
    forces = force_blocks(text)
    total_magnetization, absolute_magnetization, scf_iterations = converged_cell_history(text)
    count = min(len(propagated), len(forces), len(total_magnetization))
    if count < 2:
        raise ValueError(f"Need at least two complete QE frames, found {count}")

    # QE evaluates force record 1 at the input positions. Each printed position
    # block is the configuration propagated after that force evaluation.
    force_aligned_positions = np.concatenate((initial[None], propagated[:-1]), axis=0)[:count]
    local_magnetization = vector_frames(text, "magnetization")
    constrained_moments = vector_frames(text, "constrained moment")
    if len(local_magnetization) >= count + 1:
        local_magnetization = local_magnetization[1:count + 1]
    else:
        local_magnetization = local_magnetization[:count]
    if len(constrained_moments) >= count + 1:
        constrained_moments = constrained_moments[1:count + 1]
    else:
        constrained_moments = constrained_moments[:count]

    temperature = scalars(
        rf"^\s*temperature\s*=\s*({FLOAT})\s*K", text, re.I | re.M
    )[:count]
    pressure = scalars(rf"(?m)(?:^|\s)P\s*=\s*({FLOAT})", text)[:count] / 10.0
    if min(len(local_magnetization), len(constrained_moments), len(temperature), len(pressure)) < count:
        raise ValueError("One or more QE observables do not cover every complete force frame")

    ideal = ideal_sites_in_atom_order(initial, cell)
    dt_ps = DT_AU * 2.0 * AU_TIME_SECONDS * 1.0e12
    all_iterations = np.arange(1, count + 1, dtype=int)
    all_time = np.arange(count, dtype=float) * dt_ps
    keep = selected_indices(count, args.stride)
    payload: dict[str, np.ndarray] = {
        "input_cell_parameters": np.asarray(cell, dtype=np.float64),
        "input_cell_unit": np.asarray("angstrom", dtype="U16"),
        "initial_positions_alat": np.asarray(ideal, dtype=np.float64),
        "initial_cell_alat": np.eye(3, dtype=np.float64),
        "positions": np.asarray(force_aligned_positions[keep], dtype=np.float64),
        "positions_unit": np.full(len(keep), "crystal", dtype="U16"),
        "forces_ry_au": np.asarray(forces[:count][keep], dtype=np.float64),
        "symbols": np.asarray(symbols, dtype="U2"),
        "species": np.asarray(symbols, dtype="U2"),
        "iteration": all_iterations[keep],
        "time_ps": all_time[keep],
        "temperature_K": temperature[keep],
        "pressure_GPa": pressure[keep],
        "mag_total_vector_Bohr": total_magnetization[:count][keep],
        "mag_total_Bohr": np.linalg.norm(total_magnetization[:count][keep], axis=1),
        "abs_mag_total_Bohr": absolute_magnetization[:count][keep],
        "local_magnetization_Bohr": local_magnetization[keep],
        "printed_constrained_moments_Bohr": constrained_moments[keep],
        "scf_iterations": scf_iterations[:count][keep],
        "source_complete_frames": np.asarray([count], dtype=np.int32),
        "source_requested_frames": np.asarray([REQUESTED_STEPS], dtype=np.int32),
        "source_job_done": np.asarray(["JOB DONE" in text], dtype=bool),
    }
    args.npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.npz, **payload)
    return payload


def compute(args):
    from held.model import build_model_from_npz
    from held.phases import build_q_path
    from held.io import write_fc_csv, write_step_csv
    import hashlib

    data = prepare_trajectory(args)
    text = source_text(args)
    boundaries = list(re.finditer(r'Entering Dynamics:\s*iteration =\s*(\d+)', text))
    n = int(data['source_complete_frames'].item())
    assert len(boundaries) == n
    assert [int(m[1]) for m in boundaries] == list(range(1, n + 1))
    segments = [text[0 if i == 0 else boundaries[i-1].end():m.start()]
                for i, m in enumerate(boundaries)]
    keep = selected_indices(n, args.stride)
    def last_values(pattern):
        values = []
        for segment in segments:
            matches = re.findall(pattern, segment, re.M)
            if not matches:
                raise ValueError(f'Missing observable: {pattern}')
            values.append(number(matches[-1]))
        return np.asarray(values)[keep]
    data['temperature_K'] = last_values(rf'(?:Starting )?temperature\s*=\s*({FLOAT})\s*K')
    data['energy_Ry'] = last_values(rf'!\s+total energy\s*=\s*({FLOAT})')
    data['constraint_energy_Ry'] = last_values(rf'constraint energy \(Ry\)\s*=\s*({FLOAT})')
    data['scf_accuracy_Ry'] = last_values(rf'estimated scf accuracy\s*<\s*({FLOAT})')
    data['force_correction_Ry_Bohr'] = last_values(rf'Total SCF correction\s*=\s*({FLOAT})')
    data['total_force_Ry_Bohr'] = last_values(rf'Total force\s*=\s*({FLOAT})')
    data['force_warning'] = np.asarray(['SCF correction compared to forces is large' in s for s in segments])[keep]
    # Use the final atomic report of each completed SCF cycle, excluding initial guesses.
    for key, label in [('local_magnetization_Bohr', 'magnetization'),
                       ('printed_constrained_moments_Bohr', 'constrained moment')]:
        data[key] = np.asarray([vector_frames(s, label)[-1] for s in segments])[keep]
    assert all('convergence has been achieved' in s and 'convergence NOT achieved' not in s for s in segments)
    np.savez_compressed(args.npz, **data)
    model, dataset, metadata = build_model_from_npz('bcc', args.npz, num_shells=1)
    assert max(metadata.selected_shell_distances) < min(np.linalg.norm(model.ss_cell, axis=1)) / 2
    result = model.fit_series(dataset.positions_frac, dataset.forces_ev_ang, dataset.step_ids, verbose=True)
    write_fc_csv(DEFAULT_COEFFICIENTS, result)
    write_step_csv(DEFAULT_STEPS_CSV, result, frame_indices=np.arange(len(result.step_ids)), observables=dataset.observables)
    q, x, labels, ticks = build_q_path('bcc', model.uc_cell, path_labels=['GM','H','N','GM','P','H'], points_per_segment=args.points_per_segment)
    frequencies = np.asarray([model.dispersion_thz_from_reduced_path(c, q) for c in result.step_values])
    mean_freq = model.dispersion_thz_from_reduced_path(result.mean_values, q)
    ranks, errors, conditions = [], [], []
    for pos, force, coeff in zip(dataset.positions_frac, dataset.forces_ev_ang, result.step_values):
        design = model.build_design_matrix(model.frame_displacements_cart(pos))
        ranks.append(np.linalg.matrix_rank(design))
        conditions.append(np.linalg.cond(design))
        errors.append(np.sqrt(np.mean((design @ coeff - force.ravel())**2)))
    assert min(ranks) == len(result.labels), 'Underdetermined per-frame fit'
    assert np.isfinite(frequencies).all()
    gamma_max = float(np.max(np.abs(frequencies[:,np.linalg.norm(q,axis=1)<1e-10,:])))
    assert gamma_max < 1e-4, gamma_max
    data['fit_rmse_eV_A'] = np.asarray(errors)
    data['fit_condition'] = np.asarray(conditions)
    data['fit_rank'] = np.asarray(ranks)
    data['imaginary_path_percent'] = 100*np.mean(frequencies < -1e-4, axis=(1,2))
    np.savez_compressed(args.npz, **data)
    np.savez_compressed(DEFAULT_HEATMAP_CACHE, step_ids=dataset.step_ids, qpoints=q, x_values=x,
                        tick_labels=labels, tick_positions=ticks, step_frequencies_thz=frequencies,
                        mean_fc_frequencies_thz=mean_freq, coefficients=result.step_values,
                        coefficient_labels=result.labels, num_shells=1)
    # Export every scalar/vector observable to a readable table.
    columns = {}
    for key, value in data.items():
        if value.shape == (len(keep),) and value.dtype.kind in 'ifb': columns[key] = value
        elif key == 'mag_total_vector_Bohr':
            for j, axis in enumerate('xyz'): columns['mag_'+axis+'_Bohr'] = value[:,j]
    with (HERE/'round22_dashboard_parameters.csv').open('w') as f:
        writer=csv.writer(f); writer.writerow(columns); writer.writerows(zip(*columns.values()))
    plot_mean_dispersion('bcc', args.npz, DEFAULT_COEFFICIENTS,
                         output_data=DEFAULT_DISPERSION_DATA, output_plot=DEFAULT_DISPERSION_PNG,
                         num_shells=1, path_labels=['GM','H','N','GM','P','H'], points_per_segment=args.points_per_segment)
    summary = dict(natoms=32, repetitions=[4,2,2], frames=len(keep), source_complete_frames=n,
                   requested_frames=REQUESTED_STEPS, source_outputs=[str(args.previous_output), str(args.output)], source_job_done='JOB DONE' in text, num_shells=1,
                   cutoff_angstrom=metadata.selected_shell_distances[-1],
                   force_position_alignment='Force frame 1 uses input positions; frame k uses propagated positions k-1. Temperature is the pre-propagation value. Time starts at 0 fs.',
                   source_sha256=hashlib.sha256(text.encode()).hexdigest(),
                   model='Cubic-symmetry, first-neighbor HELD per-frame harmonic force fit; mean curve uses mean force constants.',
                   limitations=['Combined initial segment and completed 100-step restart; 162 force frames.', 'First-neighbor truncation: longer-range interactions are omitted to avoid minimum-image ambiguity in the 4x2x2 cell.',
                                'Instantaneous diagnostic phonons; not an equilibrated finite-temperature spectrum. Negative frequencies denote imaginary modes.'],
                   independent_coefficients=len(result.labels), minimum_fit_rank=int(min(ranks)),
                   maximum_condition_number=float(max(conditions)), gamma_max_abs_THz=gamma_max,
                   frequency_range_THz=[float(frequencies.min()),float(frequencies.max())],
                   fit_rmse_eV_A_range=[float(min(errors)),float(max(errors))],
                   source_output=str(args.output), source_input=str(args.input))
    args.summary.write_text(json.dumps(summary, indent=2)+'\n')
    return data, frequencies, mean_freq, x, labels, ticks


def render(args, data, frequencies, mean_freq, x, labels, ticks):
    from PIL import Image
    from io import BytesIO
    steps=data['iteration']; n=len(steps)
    cell=data['input_cell_parameters']; ideal=data['initial_positions_alat']
    delta=data['positions']-ideal; delta-=np.rint(delta); delta=delta@cell
    rms=np.sqrt(np.mean(np.sum(delta**2,axis=2),axis=1))
    maximum=np.max(np.linalg.norm(delta,axis=2),axis=1)
    colors=['#2563eb','#e87924','#059669']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold'})
    fig=plt.figure(figsize=(16,10), facecolor='#f5f7fb')
    grid=fig.add_gridspec(3,3,left=.055,right=.95,bottom=.12,top=.79,hspace=.62,wspace=.43,
                         height_ratios=[1.45,1,1])
    band=fig.add_subplot(grid[0,:2]); structure=fig.add_subplot(grid[0,2],projection='3d')
    axes=[fig.add_subplot(grid[r,c]) for r in (1,2) for c in range(3)]
    temp,mag,scf,displ,energy,fit=axes
    fig.text(.045,.958,'ROUND 22  |  HELD phonons + MD dashboard',fontsize=23,weight='bold',color='#12304d')
    fig.text(.045,.921,f'32 Fe atoms • 4×2×2 BCC • a = 2.39 Å • {n} complete force frames • steps {steps[0]}–{steps[-1]}',fontsize=11)
    fig.text(.045,.894,'1NN HELD (cutoff 2.070 Å) • Γ–H–N–Γ–P–H • dt ≈ 1 fs • SVR: 4000 K, nraise 5 • Γ-point DFT',fontsize=10)
    fig.text(.045,.869,'DFT: ecut 71/496 Ry • FD smearing 0.025334 Ry • conv_thr 10⁻⁴ Ry • β 0.01, local-TF • λ 0.1 Ry • no SOC',fontsize=10)
    status=fig.text(.045,.828,'',fontsize=11,weight='bold',color='#12304d')
    fig.text(.045,.035,'Diagnostic first-neighbor fits: longer-range forces omitted; finite trajectory. Negative ν = imaginary mode.\nForce frame k uses positions before propagation k; thermodynamic histories use that same configuration.',fontsize=9,color='#526278')
    for j in range(3): band.plot(x,mean_freq[:,j],color='#94a3b8',ls='--',lw=1.4,label='Mean force constants' if j==0 else None)
    bands=[band.plot(x,frequencies[0,:,j],color=colors[j],lw=1.7)[0] for j in range(3)]
    for v in ticks: band.axvline(v,color='#cbd5e1',lw=.7)
    band.axhline(0,color='#64748b',lw=.8); band.set_xticks(ticks,labels)
    band.set(xlim=(x[0],x[-1]),ylim=(min(-.5,frequencies.min()-1),frequencies.max()+1),ylabel='Frequency (THz)',title='Per-frame phonon dispersion')
    band.legend(fontsize=8,loc='lower right')
    cursors=[]
    def history(ax, ys, names, title, ylabel):
        for j,(y,name) in enumerate(zip(ys,names)): ax.plot(steps,y,color=colors[j%3],lw=1.4,label=name)
        ax.set(title=title,ylabel=ylabel,xlabel='Force frame / MD step',xlim=(steps[0]-.5,steps[-1]+.5))
        ax.grid(alpha=.18);ax.legend(fontsize=7,loc='best')
        cursors.append(ax.axvline(steps[0],color='#dc2626',lw=1))
    history(temp,[data['temperature_K']],['Temperature'],'Temperature / pressure','Temperature (K)')
    temp.axhline(4000,color='#64748b',ls='--',lw=.8)
    pressure=temp.twinx();pressure.plot(steps,data['pressure_GPa'],color='#e87924',alpha=.65,lw=1)
    pressure.set_ylabel('Pressure (GPa)',color='#e87924',fontsize=8)
    history(mag,list(data['mag_total_vector_Bohr'].T),['Mx','My','Mz'],'Cell magnetization','μB / cell')
    history(scf,[data['scf_iterations']],['SCF iterations'],'SCF effort / force correction','Iterations')
    corr=scf.twinx();corr.plot(steps,100*data['force_correction_Ry_Bohr']/data['total_force_Ry_Bohr'],color='#e87924',lw=1)
    corr.set_ylabel('SCF correction / force (%)',color='#e87924',fontsize=8)
    history(displ,[rms,maximum],['RMS','Maximum'],'Displacement from ideal BCC','Å')
    e=(data['energy_Ry']-data['energy_Ry'][0])*13.605693122994/32
    history(energy,[e],['Δ electronic energy'],'Energy / magnetic constraint','ΔE (eV / atom)')
    e2=energy.twinx();e2.plot(steps,data['constraint_energy_Ry'],color='#e87924',lw=1)
    e2.set_ylabel('Constraint energy (Ry)',color='#e87924',fontsize=8)
    history(fit,[data['fit_rmse_eV_A']],['HELD force RMSE'],'Fit error / imaginary branches','RMSE (eV / Å)')
    f2=fit.twinx();f2.plot(steps,data['imaginary_path_percent'],color='#e87924',lw=1)
    f2.set_ylabel('Imaginary path samples (%)',color='#e87924',fontsize=8)
    images=[]
    for i in range(n):
        for j,line in enumerate(bands): line.set_ydata(frequencies[i,:,j])
        for cursor in cursors: cursor.set_xdata([steps[i],steps[i]])
        structure.clear()
        pos=(data['positions'][i]%1)@cell
        spin=data['local_magnetization_Bohr'][i]; norm=np.linalg.norm(spin,axis=1)
        unit=spin/np.maximum(norm[:,None],1e-12)
        structure.scatter(*pos.T,c=unit[:,2],cmap='coolwarm',vmin=-1,vmax=1,s=30,depthshade=True)
        structure.quiver(*pos.T,*unit.T,length=.55,normalize=False,color='#334155',linewidth=.65)
        structure.set(xlim=(0,cell[0,0]),ylim=(0,cell[1,1]),zlim=(0,cell[2,2]),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)')
        structure.set_box_aspect(np.diag(cell).copy());structure.view_init(elev=22,azim=-58)
        structure.set_title(f'Atoms + local spins | ⟨|m|⟩ = {norm.mean():.2f} μB',fontsize=9)
        structure.tick_params(labelsize=6)
        status.set_text(f'Step {steps[i]:02d}/{steps[-1]}   •   t = {data["time_ps"][i]*1000:.1f} fs   •   T = {data["temperature_K"][i]:.0f} K   •   P = {data["pressure_GPa"][i]:.1f} GPa   •   SCF = {data["scf_iterations"][i]} iterations   •   |M| = {data["mag_total_Bohr"][i]:.2f} μB')
        buffer=BytesIO();fig.savefig(buffer,format='png',dpi=args.dpi);buffer.seek(0)
        images.append(Image.open(buffer).convert('RGB').quantize(colors=256))
        if i==0: fig.savefig(HERE/'round22_held_md_dashboard_first.png',dpi=args.dpi)
        if i==n-1: fig.savefig(args.png,dpi=args.dpi)
        if i%10==0 or i==n-1: print(f'[GIF] rendered {i+1}/{n}',flush=True)
    images[0].save(args.gif,save_all=True,append_images=images[1:],duration=round(1000/args.fps),loop=0,optimize=False,disposal=2)
    plt.close(fig)
    with Image.open(args.gif) as gif: assert gif.n_frames==n,(gif.n_frames,n)
    print(f'Saved {args.gif}',flush=True)


if __name__ == '__main__':
    args=arguments()
    data,freq,mean,x,labels,ticks=compute(args)
    render(args,data,freq,mean,x,labels,ticks)
