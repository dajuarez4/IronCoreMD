#!/usr/bin/env python3

import argparse
import csv
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


BOHR_TO_ANG = 0.529177210903
RY_TO_EV = 13.605693009
RY_BOHR_TO_EV_ANG = RY_TO_EV / BOHR_TO_ANG

PHASE_COLORS = {
    "BCC": "#0072B2",
    "FCC": "#D55E00",
    "HCP": "#009E73",
}

STATE_MARKERS = {
    "non-magnetic": "o",
    "collinear DLM": "s",
    "noncollinear": "*",
}

FEATURE_LABELS = [
    r"$V/N$",
    r"$E/N$",
    r"$T$",
    r"$P$",
    r"$F_{\mathrm{RMS}}$",
    r"$u_{\mathrm{RMS}}$",
    r"$|M|/N$",
]


@dataclass(frozen=True)
class DatasetSpec:
    phase: str
    magnetic_state: str
    path: Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the combined BCC/FCC/HCP Quantum ESPRESSO dataset space."
    )
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path("dataset"),
        help="Canonical dataset directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("dataset/figures"),
        help="Directory for figures and the source-level summary.",
    )
    parser.add_argument(
        "--poster-dir",
        type=Path,
        default=None,
        help="Optional directory receiving copies of the PNG and PDF.",
    )
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def discover_datasets(root):
    specs = []
    specs.extend(
        DatasetSpec("BCC", "non-magnetic", path)
        for path in sorted((root / "bcc" / "non-mag").glob("*.npz"))
    )
    specs.extend(
        DatasetSpec("FCC", "non-magnetic", path)
        for path in sorted((root / "fcc" / "non-mag").glob("*.npz"))
    )
    specs.extend(
        DatasetSpec("HCP", "non-magnetic", path)
        for path in sorted((root / "hcp").glob("*.npz"))
    )

    magnetic = [
        DatasetSpec(
            "BCC",
            "collinear DLM",
            root / "bcc" / "magnetic-collinear" / "simulation.npz",
        ),
        DatasetSpec(
            "BCC",
            "noncollinear",
            root / "bcc" / "magnetic-non_coll" / "simulation.npz",
        ),
    ]
    specs.extend(spec for spec in magnetic if spec.path.exists())
    return specs


def read_metadata(data):
    if "metadata_json" not in data:
        return {}
    return json.loads(str(data["metadata_json"]))


def infer_target_temperature(path, values):
    match = re.search(r"_(\d{4,5})K", path.stem)
    if match:
        return float(match.group(1))
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else np.nan


def cell_from_archive(data, metadata):
    alat_ang = float(metadata["alat_bohr"]) * BOHR_TO_ANG
    cell_alat = np.asarray(data["initial_cell_alat"], dtype=np.float64)
    return cell_alat * alat_ang, cell_alat, alat_ang


def periodic_displacement_rms(data, valid_mask, cell_ang, cell_alat):
    positions = np.asarray(data["positions"], dtype=np.float64)[valid_mask]
    units = np.asarray(data["positions_unit"])[valid_mask]
    initial_alat = np.asarray(data["initial_positions_alat"], dtype=np.float64)
    initial_fractional = initial_alat @ np.linalg.inv(cell_alat)
    output = np.full(positions.shape[0], np.nan, dtype=np.float64)

    for unit in np.unique(units):
        selector = units == unit
        unit_name = str(unit).strip().lower()
        selected = positions[selector]
        if unit_name == "crystal":
            fractional = selected
        elif unit_name == "alat":
            fractional = selected @ np.linalg.inv(cell_alat)
        elif unit_name == "angstrom":
            fractional = selected @ np.linalg.inv(cell_ang)
        elif unit_name == "bohr":
            fractional = (selected * BOHR_TO_ANG) @ np.linalg.inv(cell_ang)
        else:
            continue

        delta_fractional = fractional - initial_fractional
        delta_fractional -= np.rint(delta_fractional)
        delta_cartesian = delta_fractional @ cell_ang
        output[selector] = np.sqrt(np.mean(np.sum(delta_cartesian**2, axis=2), axis=1))
    return output


def load_archive(spec):
    with np.load(spec.path, allow_pickle=False) as data:
        metadata = read_metadata(data)
        natoms = int(metadata.get("natoms", len(data["symbols"])))
        energy = np.asarray(data["internal_energy_ry"], dtype=np.float64)
        if not np.isfinite(energy).any():
            energy = np.asarray(data["energy_ry"], dtype=np.float64)
        temperature = np.asarray(data["temperature_K"], dtype=np.float64)
        pressure = np.asarray(data["pressure_GPa"], dtype=np.float64)
        forces = np.asarray(data["forces_ry_au"], dtype=np.float64)

        valid = (
            np.isfinite(energy)
            & np.isfinite(temperature)
            & np.isfinite(pressure)
            & np.isfinite(forces).all(axis=(1, 2))
        )
        if "frame_valid" in data:
            valid &= np.asarray(data["frame_valid"], dtype=bool)

        cell_ang, cell_alat, _ = cell_from_archive(data, metadata)
        volume_per_atom = abs(np.linalg.det(cell_ang)) / natoms
        energy_per_atom = energy[valid] * RY_TO_EV / natoms
        temperature = temperature[valid]
        pressure = pressure[valid]
        forces = forces[valid] * RY_BOHR_TO_EV_ANG
        force_rms = np.sqrt(np.mean(forces**2, axis=(1, 2)))
        displacement_rms = periodic_displacement_rms(
            data, valid, cell_ang, cell_alat
        )

        magnetization = np.asarray(
            data["abs_mag_total_Bohr"]
            if "abs_mag_total_Bohr" in data
            else np.zeros(len(valid)),
            dtype=np.float64,
        )[valid]
        magnetization = np.nan_to_num(magnetization, nan=0.0) / natoms
        if spec.magnetic_state == "non-magnetic":
            magnetization[:] = 0.0

    features = np.column_stack(
        [
            np.full(valid.sum(), volume_per_atom),
            energy_per_atom,
            temperature,
            pressure,
            force_rms,
            displacement_rms,
            magnetization,
        ]
    )
    finite_features = np.isfinite(features).all(axis=1)
    features = features[finite_features]

    target_temperature = infer_target_temperature(spec.path, temperature)
    summary = {
        "phase": spec.phase,
        "magnetic_state": spec.magnetic_state,
        "source": spec.path.stem,
        "path": str(spec.path),
        "natoms": natoms,
        "total_frames": len(valid),
        "valid_frames": int(features.shape[0]),
        "volume_per_atom_A3": volume_per_atom,
        "target_temperature_K": target_temperature,
        "mean_temperature_K": float(np.mean(features[:, 2])),
        "min_temperature_K": float(np.min(features[:, 2])),
        "max_temperature_K": float(np.max(features[:, 2])),
        "mean_pressure_GPa": float(np.mean(features[:, 3])),
        "min_pressure_GPa": float(np.min(features[:, 3])),
        "max_pressure_GPa": float(np.max(features[:, 3])),
        "mean_energy_eV_atom": float(np.mean(features[:, 1])),
        "min_energy_eV_atom": float(np.min(features[:, 1])),
        "max_energy_eV_atom": float(np.max(features[:, 1])),
        "mean_force_rms_eV_A": float(np.mean(features[:, 4])),
        "mean_displacement_rms_A": float(np.mean(features[:, 5])),
        "mean_abs_magnetization_Bohr_atom": float(np.mean(features[:, 6])),
    }
    return features, summary


def robust_pca(features):
    lower = np.nanpercentile(features, 1.0, axis=0)
    upper = np.nanpercentile(features, 99.0, axis=0)
    sparse_features = upper <= lower
    lower[sparse_features] = np.nanmin(features[:, sparse_features], axis=0)
    upper[sparse_features] = np.nanmax(features[:, sparse_features], axis=0)
    clipped = np.clip(features, lower, upper)
    center = np.median(clipped, axis=0)
    scale = np.nanpercentile(clipped, 75.0, axis=0) - np.nanpercentile(
        clipped, 25.0, axis=0
    )
    fallback_scale = np.nanmax(clipped, axis=0) - np.nanmin(clipped, axis=0)
    scale[scale == 0.0] = fallback_scale[scale == 0.0]
    scale[scale == 0.0] = 1.0
    standardized = (clipped - center) / scale
    covariance = np.cov(standardized, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    scores = standardized @ eigenvectors[:, :2]
    explained = eigenvalues[:2] / eigenvalues.sum()
    return scores, eigenvectors[:, :2], explained


def set_publication_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.titlesize": 18,
            "axes.labelsize": 17,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def phase_state_label(phase, state):
    return phase if state == "non-magnetic" else f"{phase} — {state}"


def add_panel_label(axis, label):
    axis.text(
        -0.10,
        1.01,
        label,
        transform=axis.transAxes,
        fontsize=19,
        fontweight="bold",
        va="bottom",
    )


def plot_dataset_space(features, labels, summaries, output_base, dpi):
    set_publication_style()
    scores, loadings, explained = robust_pca(features)
    phases = np.asarray([item[0] for item in labels])
    states = np.asarray([item[1] for item in labels])

    figure = plt.figure(figsize=(15.8, 10.2), constrained_layout=True)
    figure.get_layout_engine().set(rect=(0.0, 0.0, 1.0, 0.92))
    grid = figure.add_gridspec(2, 2, width_ratios=(1.18, 1.0))
    axis_pca = figure.add_subplot(grid[0, 0])
    axis_pt = figure.add_subplot(grid[0, 1])
    axis_ev = figure.add_subplot(grid[1, 0])
    axis_count = figure.add_subplot(grid[1, 1])

    plotting_order = [
        ("BCC", "non-magnetic"),
        ("FCC", "non-magnetic"),
        ("HCP", "non-magnetic"),
        ("BCC", "collinear DLM"),
        ("BCC", "noncollinear"),
    ]
    for phase, state in plotting_order:
        selector = (phases == phase) & (states == state)
        if not selector.any():
            continue
        magnetic = state != "non-magnetic"
        axis_pca.scatter(
            scores[selector, 0],
            scores[selector, 1],
            s=34 if magnetic else 4,
            marker=STATE_MARKERS[state],
            c=PHASE_COLORS[phase],
            alpha=0.90 if magnetic else 0.18,
            linewidths=0.45 if magnetic else 0,
            edgecolors="white" if magnetic else "none",
            rasterized=not magnetic,
            zorder=4 if magnetic else 2,
            label=phase_state_label(phase, state),
        )

    loading_scale = np.percentile(np.linalg.norm(scores, axis=1), 85) * 0.62
    loading_offsets = [(-7, 8), (7, 12), (7, -9), (7, -12), (7, -8), (-7, -9), (7, 8)]
    for index, label in enumerate(FEATURE_LABELS):
        if np.linalg.norm(loadings[index]) < 0.12:
            continue
        dx, dy = loadings[index] * loading_scale
        axis_pca.annotate(
            "",
            xy=(dx, dy),
            xytext=(0, 0),
            arrowprops={"arrowstyle": "-|>", "color": "#333333", "lw": 0.9},
            zorder=5,
        )
        axis_pca.annotate(
            label,
            xy=(dx, dy),
            xytext=loading_offsets[index],
            textcoords="offset points",
            fontsize=11,
            color="#222222",
            ha="left" if loading_offsets[index][0] > 0 else "right",
            va="center",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.72},
        )
    axis_pca.axhline(0, color="#B7B7B7", lw=0.6, zorder=0)
    axis_pca.axvline(0, color="#B7B7B7", lw=0.6, zorder=0)
    axis_pca.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    axis_pca.set_ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
    axis_pca.set_title("Thermodynamic–configurational dataset space", loc="left")
    legend_handles, legend_labels = axis_pca.get_legend_handles_labels()
    phase_legend = axis_pca.legend(
        legend_handles[:3],
        legend_labels[:3],
        loc="lower left",
        bbox_to_anchor=(0.015, 0.02),
        frameon=True,
        framealpha=0.94,
        borderpad=0.6,
        markerscale=1.5,
        fontsize=13,
    )
    axis_pca.add_artist(phase_legend)
    axis_pca.legend(
        legend_handles[3:],
        legend_labels[3:],
        loc="lower right",
        bbox_to_anchor=(0.985, 0.02),
        frameon=True,
        framealpha=0.94,
        borderpad=0.6,
        markerscale=1.5,
        fontsize=13,
    )
    add_panel_label(axis_pca, "a")

    for summary in summaries:
        phase = summary["phase"]
        state = summary["magnetic_state"]
        magnetic = state != "non-magnetic"
        axis_pt.scatter(
            summary["mean_pressure_GPa"],
            summary["mean_temperature_K"],
            s=18 + 1.8 * np.sqrt(summary["valid_frames"]),
            c=PHASE_COLORS[phase],
            marker=STATE_MARKERS[state],
            alpha=0.92 if magnetic else 0.68,
            linewidths=0.65 if magnetic else 0.25,
            edgecolors="black" if magnetic else "white",
            zorder=4 if magnetic else 2,
        )
    axis_pt.set_xlabel("Mean pressure (GPa)")
    axis_pt.set_ylabel("Mean ionic temperature (K)")
    axis_pt.set_title("Pressure–temperature coverage", loc="left")
    axis_pt.grid(color="#D9D9D9", lw=0.55, alpha=0.65)
    add_panel_label(axis_pt, "b")

    global_minimum = np.min(features[:, 1])
    for phase, state in plotting_order:
        selector = (phases == phase) & (states == state)
        if not selector.any():
            continue
        magnetic = state != "non-magnetic"
        axis_ev.scatter(
            features[selector, 0],
            features[selector, 1] - global_minimum,
            s=28 if magnetic else 3,
            c=PHASE_COLORS[phase],
            marker=STATE_MARKERS[state],
            alpha=0.82 if magnetic else 0.12,
            linewidths=0.35 if magnetic else 0,
            edgecolors="white" if magnetic else "none",
            rasterized=not magnetic,
            zorder=4 if magnetic else 2,
        )
    axis_ev.set_xlabel(r"Volume per atom ($\mathrm{\AA^3/atom}$)")
    axis_ev.set_ylabel(r"$E - E_{\min}$ (eV/atom)")
    axis_ev.set_title("Energy–volume coverage", loc="left")
    axis_ev.grid(color="#D9D9D9", lw=0.55, alpha=0.65)
    add_panel_label(axis_ev, "c")

    composition = []
    for phase, state in plotting_order:
        selected = [
            summary
            for summary in summaries
            if summary["phase"] == phase and summary["magnetic_state"] == state
        ]
        if not selected:
            continue
        composition.append(
            (
                phase_state_label(phase, state),
                phase,
                state,
                sum(item["valid_frames"] for item in selected),
                len(selected),
            )
        )
    composition.reverse()
    positions = np.arange(len(composition))
    counts = np.asarray([item[3] for item in composition])
    bars = axis_count.barh(
        positions,
        counts,
        color=[PHASE_COLORS[item[1]] for item in composition],
        edgecolor="white",
        linewidth=0.8,
        alpha=0.88,
    )
    axis_count.set_xscale("log")
    axis_count.set_yticks(positions, [item[0] for item in composition])
    axis_count.set_xlabel("Valid configurations (log scale)")
    axis_count.set_title("Dataset composition", loc="left")
    axis_count.grid(axis="x", color="#D9D9D9", lw=0.55, alpha=0.65)
    for bar, item in zip(bars, composition):
        axis_count.text(
            item[3] * 1.13,
            bar.get_y() + bar.get_height() / 2,
            f"{item[3]:,} frames · {item[4]} source{'s' if item[4] != 1 else ''}",
            va="center",
            fontsize=11,
        )
    axis_count.set_xlim(1, max(counts) * 7.0)
    add_panel_label(axis_count, "d")

    total_frames = int(sum(item["valid_frames"] for item in summaries))
    figure.suptitle(
        "Dataset Atlas",
        fontsize=27,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.947,
        f"{total_frames:,} valid configurations",
        ha="center",
        va="top",
        fontsize=18,
        color="#444444",
    )
    for suffix in ("png", "pdf", "svg"):
        figure.savefig(
            output_base.with_suffix(f".{suffix}"),
            dpi=dpi if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(figure)


def write_summary(path, summaries):
    if not summaries:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)


def main():
    args = parse_args()
    specs = discover_datasets(args.dataset_root)
    if not specs:
        raise SystemExit(f"No NPZ datasets found below {args.dataset_root}")

    feature_blocks = []
    labels = []
    summaries = []
    for spec in specs:
        features, summary = load_archive(spec)
        if features.size == 0:
            continue
        feature_blocks.append(features)
        labels.extend([(spec.phase, spec.magnetic_state)] * len(features))
        summaries.append(summary)

    combined_features = np.vstack(feature_blocks)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_base = args.output_dir / "iron_all_phases_dataset_space"
    plot_dataset_space(combined_features, labels, summaries, output_base, args.dpi)
    write_summary(
        args.output_dir / "iron_all_phases_dataset_summary.csv", summaries
    )

    if args.poster_dir is not None:
        args.poster_dir.mkdir(parents=True, exist_ok=True)
        for suffix in ("png", "pdf"):
            source = output_base.with_suffix(f".{suffix}")
            shutil.copy2(source, args.poster_dir / source.name)

    print(f"Loaded {len(summaries)} trajectories")
    print(f"Valid configurations: {len(combined_features)}")
    print(f"Figure: {output_base.with_suffix('.png')}")


if __name__ == "__main__":
    main()
