#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from plot_all_phase_dataset_space import (
    BOHR_TO_ANG,
    PHASE_COLORS,
    RY_BOHR_TO_EV_ANG,
    STATE_MARKERS,
    cell_from_archive,
    discover_datasets,
    infer_target_temperature,
    read_metadata,
    set_publication_style,
)


CATEGORY_ORDER = [
    ("BCC", "non-magnetic"),
    ("BCC", "collinear DLM"),
    ("BCC", "noncollinear"),
    ("FCC", "non-magnetic"),
    ("HCP", "non-magnetic"),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create EOS coverage and quality-control figures for IronCoreMD."
    )
    parser.add_argument("--dataset-root", type=Path, default=Path("dataset"))
    parser.add_argument("--output-dir", type=Path, default=Path("dataset/figures"))
    parser.add_argument("--poster-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def category_label(phase, state):
    if state == "non-magnetic":
        return f"{phase} nonmagnetic"
    if state == "collinear DLM":
        return "BCC collinear DLM"
    return "BCC noncollinear"


def add_panel_label(axis, label):
    axis.text(
        -0.11,
        1.01,
        label,
        transform=axis.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
    )


def robust_force_outliers(force_rms):
    finite = force_rms[np.isfinite(force_rms)]
    if finite.size < 8:
        return np.zeros(force_rms.shape, dtype=bool), np.nan
    first, third = np.percentile(finite, [25.0, 75.0])
    threshold = third + 3.0 * (third - first)
    return force_rms > threshold, threshold


def load_quality_archive(spec):
    with np.load(spec.path, allow_pickle=False) as data:
        metadata = read_metadata(data)
        natoms = int(metadata.get("natoms", len(data["symbols"])))
        temperature = np.asarray(data["temperature_K"], dtype=np.float64)
        pressure = np.asarray(data["pressure_GPa"], dtype=np.float64)
        forces = np.asarray(data["forces_ry_au"], dtype=np.float64)
        energy_key = "internal_energy_ry"
        energy = np.asarray(data[energy_key], dtype=np.float64)
        if not np.isfinite(energy).any():
            energy_key = "energy_ry"
            energy = np.asarray(data[energy_key], dtype=np.float64)

        total_frames = len(temperature)
        finite_energy = np.isfinite(energy)
        finite_temperature = np.isfinite(temperature)
        finite_pressure = np.isfinite(pressure)
        finite_forces = np.isfinite(forces).all(axis=(1, 2))
        finite_positions = np.ones(total_frames, dtype=bool)
        if "position_frame_valid" in data:
            finite_positions = np.asarray(data["position_frame_valid"], dtype=bool)
        archive_valid = np.ones(total_frames, dtype=bool)
        if "frame_valid" in data:
            archive_valid = np.asarray(data["frame_valid"], dtype=bool)

        valid = (
            archive_valid
            & finite_positions
            & finite_energy
            & finite_temperature
            & finite_pressure
            & finite_forces
        )
        cell_ang, _, _ = cell_from_archive(data, metadata)
        volume = abs(np.linalg.det(cell_ang)) / natoms
        force_rms = np.full(total_frames, np.nan)
        force_rms[finite_forces] = np.sqrt(
            np.mean((forces[finite_forces] * RY_BOHR_TO_EV_ANG) ** 2, axis=(1, 2))
        )
        outliers, outlier_threshold = robust_force_outliers(force_rms[valid])
        valid_indices = np.flatnonzero(valid)
        outlier_full = np.zeros(total_frames, dtype=bool)
        outlier_full[valid_indices] = outliers
        target_temperature = infer_target_temperature(spec.path, temperature)

    valid_temperature = temperature[valid]
    valid_pressure = pressure[valid]
    record = {
        "phase": spec.phase,
        "magnetic_state": spec.magnetic_state,
        "category": category_label(spec.phase, spec.magnetic_state),
        "source": spec.path.stem,
        "path": str(spec.path),
        "natoms": natoms,
        "total_frames": total_frames,
        "valid_frames": int(valid.sum()),
        "invalid_frames": int(total_frames - valid.sum()),
        "valid_fraction": float(valid.mean()) if total_frames else np.nan,
        "finite_energy_fraction": float(finite_energy.mean()),
        "finite_temperature_fraction": float(finite_temperature.mean()),
        "finite_pressure_fraction": float(finite_pressure.mean()),
        "finite_force_fraction": float(finite_forces.mean()),
        "finite_position_fraction": float(finite_positions.mean()),
        "force_outlier_frames": int(outlier_full.sum()),
        "force_outlier_fraction": float(outlier_full.sum() / max(valid.sum(), 1)),
        "force_outlier_threshold_eV_A": float(outlier_threshold),
        "median_force_rms_eV_A": float(np.nanmedian(force_rms[valid])),
        "target_temperature_K": float(target_temperature),
        "median_temperature_K": float(np.nanmedian(valid_temperature)),
        "median_temperature_deviation_K": float(
            np.nanmedian(valid_temperature - target_temperature)
        ),
        "pressure_min_GPa": float(np.nanmin(valid_pressure)),
        "pressure_median_GPa": float(np.nanmedian(valid_pressure)),
        "pressure_max_GPa": float(np.nanmax(valid_pressure)),
        "volume_per_atom_A3": float(volume),
    }
    frames = {
        "phase": spec.phase,
        "state": spec.magnetic_state,
        "volume": np.full(valid.sum(), volume),
        "temperature": valid_temperature,
        "pressure": valid_pressure,
        "force_rms": force_rms[valid],
        "temperature_deviation": valid_temperature - target_temperature,
        "force_outlier": outlier_full[valid],
    }
    return record, frames


def split_curve(x_values):
    if len(x_values) < 3:
        return [np.arange(len(x_values))]
    gaps = np.diff(x_values)
    positive = gaps[gaps > 1.0e-8]
    typical = np.median(positive) if positive.size else 0.0
    threshold = max(0.45, 4.0 * typical)
    return np.split(np.arange(len(x_values)), np.flatnonzero(gaps > threshold) + 1)


def plot_eos(frame_sets, output_base, dpi):
    set_publication_style()
    all_temperatures = np.concatenate([item["temperature"] for item in frame_sets])
    finite_temperatures = all_temperatures[np.isfinite(all_temperatures)]
    lower, upper = np.percentile(finite_temperatures, [1.0, 99.0])
    normalization = Normalize(lower, upper)
    colormap = plt.get_cmap("plasma")

    figure, axes = plt.subplots(1, 3, figsize=(16.2, 5.6), constrained_layout=True)
    figure.get_layout_engine().set(rect=(0.0, 0.0, 0.985, 0.88))
    rng = np.random.default_rng(20260721)
    state_legend = {}

    for panel_index, (axis, phase) in enumerate(zip(axes, ("BCC", "FCC", "HCP"))):
        selected_sets = [item for item in frame_sets if item["phase"] == phase]
        volumes = np.concatenate([item["volume"] for item in selected_sets])
        pressures = np.concatenate([item["pressure"] for item in selected_sets])
        temperatures = np.concatenate([item["temperature"] for item in selected_sets])
        sample_size = min(6000, len(volumes))
        sample = rng.choice(len(volumes), sample_size, replace=False)
        axis.scatter(
            volumes[sample],
            pressures[sample],
            c=temperatures[sample],
            cmap=colormap,
            norm=normalization,
            s=5,
            alpha=0.10,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

        unique_volumes = np.unique(np.round(volumes, 6))
        medians = []
        lower_band = []
        upper_band = []
        for volume in unique_volumes:
            selector = np.isclose(volumes, volume, atol=1.0e-5)
            medians.append(np.nanmedian(pressures[selector]))
            lower_band.append(np.nanpercentile(pressures[selector], 16.0))
            upper_band.append(np.nanpercentile(pressures[selector], 84.0))
        medians = np.asarray(medians)
        lower_band = np.asarray(lower_band)
        upper_band = np.asarray(upper_band)
        for segment in split_curve(unique_volumes):
            axis.fill_between(
                unique_volumes[segment],
                lower_band[segment],
                upper_band[segment],
                color=PHASE_COLORS[phase],
                alpha=0.18,
                linewidth=0,
                zorder=2,
            )
            axis.plot(
                unique_volumes[segment],
                medians[segment],
                color=PHASE_COLORS[phase],
                lw=2.0,
                zorder=3,
            )

        for item in selected_sets:
            state = item["state"]
            marker = STATE_MARKERS[state]
            point = axis.scatter(
                np.nanmedian(item["volume"]),
                np.nanmedian(item["pressure"]),
                c=[np.nanmedian(item["temperature"])],
                cmap=colormap,
                norm=normalization,
                marker=marker,
                s=72 if state != "non-magnetic" else 26,
                edgecolors="black" if state != "non-magnetic" else "white",
                linewidths=0.8 if state != "non-magnetic" else 0.35,
                zorder=5,
            )
            state_legend[state] = point

        axis.set_title(f"{phase} iron", loc="left", fontweight="bold")
        axis.set_xlabel(r"Atomic volume ($\mathrm{\AA^3/atom}$)")
        if panel_index == 0:
            axis.set_ylabel("Pressure (GPa)")
        axis.grid(color="#D8D8D8", lw=0.55, alpha=0.7)
        add_panel_label(axis, chr(ord("a") + panel_index))

    colorbar = figure.colorbar(
        plt.cm.ScalarMappable(norm=normalization, cmap=colormap),
        ax=axes,
        location="right",
        shrink=0.89,
        pad=0.012,
    )
    colorbar.set_label("Instantaneous ionic temperature (K)")
    handles = [
        Line2D([0], [0], color=PHASE_COLORS["BCC"], lw=2.0, label="Median pressure"),
        Line2D(
            [0],
            [0],
            color=PHASE_COLORS["BCC"],
            lw=7,
            alpha=0.18,
            label="16th–84th percentile",
        ),
    ]
    for state in ("non-magnetic", "collinear DLM", "noncollinear"):
        if state not in state_legend:
            continue
        handles.append(
            Line2D(
                [0],
                [0],
                marker=STATE_MARKERS[state],
                color="none",
                markerfacecolor="#777777",
                markeredgecolor="black" if state != "non-magnetic" else "white",
                markersize=7 if state != "non-magnetic" else 6,
                label=state,
            )
        )
    axes[0].legend(handles=handles, loc="best", frameon=True, framealpha=0.94)
    figure.suptitle(
        "Equation-of-state coverage of the IronCoreMD dataset",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.925,
        "Frame-level pressure–volume sampling; curves summarize all temperatures at each sampled volume",
        ha="center",
        fontsize=10.5,
        color="#444444",
    )
    save_figure(figure, output_base, dpi)


def grouped_records(records):
    return [
        [
            record
            for record in records
            if (record["phase"], record["magnetic_state"]) == category
        ]
        for category in CATEGORY_ORDER
        if any(
            (record["phase"], record["magnetic_state"]) == category
            for record in records
        )
    ]


def plot_quality(records, frame_sets, output_base, dpi):
    set_publication_style()
    groups = grouped_records(records)
    labels = [group[0]["category"] for group in groups]
    colors = [PHASE_COLORS[group[0]["phase"]] for group in groups]
    positions = np.arange(len(groups))

    figure, axes = plt.subplots(2, 3, figsize=(16.2, 10.1), constrained_layout=True)
    figure.get_layout_engine().set(rect=(0.0, 0.0, 1.0, 0.91))
    axis_valid, axis_force, axis_temperature, axis_pressure, axis_length, axis_matrix = axes.flat

    valid_counts = np.asarray([sum(item["valid_frames"] for item in group) for group in groups])
    invalid_counts = np.asarray([sum(item["invalid_frames"] for item in group) for group in groups])
    axis_valid.bar(positions, valid_counts, color=colors, edgecolor="white", label="Valid")
    axis_valid.bar(
        positions,
        invalid_counts,
        bottom=valid_counts,
        color="#B8B8B8",
        edgecolor="white",
        label="Invalid",
    )
    for position, valid_count, invalid_count in zip(positions, valid_counts, invalid_counts):
        fraction = valid_count / max(valid_count + invalid_count, 1)
        axis_valid.text(
            position,
            valid_count + invalid_count,
            f"{fraction * 100:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    axis_valid.set_yscale("log")
    axis_valid.set_ylabel("Frames (log scale)")
    axis_valid.set_title("Frame validity", loc="left")
    axis_valid.legend(frameon=False, ncol=2)
    add_panel_label(axis_valid, "a")

    force_values = []
    outlier_rates = []
    temperature_deviations = []
    for group in groups:
        keys = {(item["phase"], item["magnetic_state"]) for item in group}
        matching_frames = [
            item for item in frame_sets if (item["phase"], item["state"]) in keys
        ]
        force = np.concatenate([item["force_rms"] for item in matching_frames])
        deviation = np.concatenate(
            [item["temperature_deviation"] for item in matching_frames]
        )
        force_values.append(force)
        temperature_deviations.append(deviation)
        outlier_rates.append(
            sum(item["force_outlier_frames"] for item in group)
            / max(sum(item["valid_frames"] for item in group), 1)
        )

    box = axis_force.boxplot(
        force_values,
        positions=positions,
        widths=0.62,
        patch_artist=True,
        showfliers=False,
        whis=(1, 99),
        medianprops={"color": "black", "lw": 1.2},
    )
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.68)
    for position, rate in zip(positions, outlier_rates):
        axis_force.text(
            position,
            axis_force.get_ylim()[1],
            f"{rate * 100:.2f}%",
            ha="center",
            va="top",
            fontsize=8,
        )
    axis_force.set_yscale("log")
    axis_force.set_ylabel(r"Frame force RMS ($\mathrm{eV/\AA}$)")
    axis_force.set_title("Force distributions and robust outlier rates", loc="left")
    add_panel_label(axis_force, "b")

    temperature_box = axis_temperature.boxplot(
        temperature_deviations,
        positions=positions,
        widths=0.62,
        patch_artist=True,
        showfliers=False,
        whis=(5, 95),
        medianprops={"color": "black", "lw": 1.2},
    )
    for patch, color in zip(temperature_box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.68)
    axis_temperature.axhline(0, color="#444444", lw=0.8, ls="--")
    axis_temperature.set_ylabel(r"$T-T_{\mathrm{target}}$ (K)")
    axis_temperature.set_title("Temperature control", loc="left")
    add_panel_label(axis_temperature, "c")

    for position, (group, color) in enumerate(zip(groups, colors)):
        for record in group:
            axis_pressure.plot(
                [record["pressure_min_GPa"], record["pressure_max_GPa"]],
                [position, position],
                color=color,
                alpha=0.20,
                lw=1.0,
            )
            axis_pressure.scatter(
                record["pressure_median_GPa"],
                position,
                color=color,
                s=10,
                alpha=0.45,
                linewidths=0,
            )
        all_min = np.percentile([item["pressure_min_GPa"] for item in group], 5)
        all_max = np.percentile([item["pressure_max_GPa"] for item in group], 95)
        all_median = np.median([item["pressure_median_GPa"] for item in group])
        axis_pressure.plot([all_min, all_max], [position, position], color=color, lw=4)
        axis_pressure.scatter(
            all_median,
            position,
            color="white",
            edgecolor=color,
            s=42,
            linewidth=1.5,
            zorder=4,
        )
    axis_pressure.set_yticks(positions, labels)
    axis_pressure.set_xlabel("Pressure (GPa)")
    axis_pressure.set_title("Trajectory pressure ranges", loc="left")
    axis_pressure.grid(axis="x", color="#D8D8D8", lw=0.55, alpha=0.7)
    add_panel_label(axis_pressure, "d")

    lengths = [[item["total_frames"] for item in group] for group in groups]
    length_box = axis_length.boxplot(
        lengths,
        positions=positions,
        widths=0.62,
        patch_artist=True,
        showfliers=True,
        medianprops={"color": "black", "lw": 1.2},
        flierprops={"marker": ".", "markersize": 3, "alpha": 0.5},
    )
    for patch, color in zip(length_box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.68)
    for position, group in enumerate(groups):
        axis_length.text(
            position,
            max(lengths[position]),
            f"n={len(group)}",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    axis_length.set_ylabel("Frames per trajectory")
    axis_length.set_title("Trajectory-length distribution", loc="left")
    add_panel_label(axis_length, "e")

    completeness_fields = [
        ("valid_fraction", "Fully valid"),
        ("finite_energy_fraction", "Energy"),
        ("finite_temperature_fraction", "Temperature"),
        ("finite_pressure_fraction", "Pressure"),
        ("finite_force_fraction", "Forces"),
        ("finite_position_fraction", "Positions"),
    ]
    matrix = np.asarray(
        [
            [np.mean([item[field] for item in group]) for field, _ in completeness_fields]
            for group in groups
        ]
    )
    image = axis_matrix.imshow(matrix * 100.0, cmap="YlGnBu", vmin=90.0, vmax=100.0)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column] * 100.0
            axis_matrix.text(
                column,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if value > 97.0 else "black",
            )
    axis_matrix.set_xticks(
        np.arange(len(completeness_fields)),
        [item[1] for item in completeness_fields],
        rotation=35,
        ha="right",
    )
    axis_matrix.set_yticks(np.arange(len(labels)), labels)
    axis_matrix.set_title("Mean field completeness (%)", loc="left")
    figure.colorbar(image, ax=axis_matrix, shrink=0.74, label="Completeness (%)")
    add_panel_label(axis_matrix, "f")

    for axis in (axis_valid, axis_force, axis_temperature, axis_length):
        axis.set_xticks(positions, labels, rotation=24, ha="right")
        axis.grid(axis="y", color="#D8D8D8", lw=0.55, alpha=0.7)

    total_frames = sum(record["total_frames"] for record in records)
    valid_frames = sum(record["valid_frames"] for record in records)
    figure.suptitle(
        "IronCoreMD dataset quality and reproducibility dashboard",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.945,
        f"{len(records)} trajectories · {valid_frames:,}/{total_frames:,} fully valid frames · force outliers use Q3 + 3×IQR per trajectory",
        ha="center",
        fontsize=10.5,
        color="#444444",
    )
    save_figure(figure, output_base, dpi)


def save_figure(figure, output_base, dpi):
    for suffix in ("png", "pdf", "svg"):
        figure.savefig(
            output_base.with_suffix(f".{suffix}"),
            dpi=dpi if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(figure)


def write_quality_csv(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    args = parse_args()
    specs = discover_datasets(args.dataset_root)
    if not specs:
        raise SystemExit(f"No NPZ datasets found below {args.dataset_root}")

    records = []
    frame_sets = []
    for spec in specs:
        record, frames = load_quality_archive(spec)
        if record["valid_frames"] == 0:
            continue
        records.append(record)
        frame_sets.append(frames)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    eos_base = args.output_dir / "iron_all_phases_eos_dataset_map"
    quality_base = args.output_dir / "iron_all_phases_dataset_quality_dashboard"
    plot_eos(frame_sets, eos_base, args.dpi)
    plot_quality(records, frame_sets, quality_base, args.dpi)
    quality_csv = args.output_dir / "iron_all_phases_dataset_quality.csv"
    write_quality_csv(quality_csv, records)

    if args.poster_dir is not None:
        args.poster_dir.mkdir(parents=True, exist_ok=True)
        for base in (eos_base, quality_base):
            for suffix in ("png", "pdf"):
                source = base.with_suffix(f".{suffix}")
                shutil.copy2(source, args.poster_dir / source.name)

    print(f"Loaded {len(records)} trajectories")
    print(f"EOS map: {eos_base.with_suffix('.png')}")
    print(f"Quality dashboard: {quality_base.with_suffix('.png')}")
    print(f"Quality table: {quality_csv}")


if __name__ == "__main__":
    main()
