#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


BACKGROUND = np.array([7, 17, 30], dtype=np.uint8)
COLOR_STOPS = np.array(
    [
        [23, 55, 94],
        [8, 126, 139],
        [86, 207, 225],
        [255, 209, 102],
        [240, 93, 94],
    ],
    dtype=np.float64,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Animate a large BCC Fe supercell with collective wave motion."
    )
    parser.add_argument("--cells", type=int, default=18)
    parser.add_argument("--lattice-parameter", type=float, default=2.55)
    parser.add_argument("--frames", type=int, default=60)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--output-dir", type=Path, default=Path("dataset/figures"))
    parser.add_argument("--poster-dir", type=Path, default=None)
    return parser.parse_args()


def bcc_supercell(cells, lattice_parameter):
    grid = np.indices((cells, cells, cells), dtype=np.float64).reshape(3, -1).T
    positions = np.vstack((grid, grid + 0.5)) * lattice_parameter
    box_length = cells * lattice_parameter
    return positions - 0.5 * box_length, box_length


def rotation_matrix():
    azimuth = np.deg2rad(43.0)
    elevation = np.deg2rad(24.0)
    rotation_z = np.array(
        [
            [np.cos(azimuth), -np.sin(azimuth), 0.0],
            [np.sin(azimuth), np.cos(azimuth), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    rotation_x = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(elevation), -np.sin(elevation)],
            [0.0, np.sin(elevation), np.cos(elevation)],
        ]
    )
    return rotation_z.T @ rotation_x.T


def wave_displacement(reference, box_length, phase):
    normalized = reference / box_length
    wave_one = np.sin(2.0 * np.pi * (1.10 * normalized[:, 0] + 0.48 * normalized[:, 2]) - phase)
    wave_two = np.sin(2.0 * np.pi * (0.42 * normalized[:, 1] - 0.85 * normalized[:, 2]) + 1.35 * phase)
    radius = np.sqrt(normalized[:, 0] ** 2 + normalized[:, 1] ** 2)
    radial_wave = np.sin(2.0 * np.pi * 2.4 * radius - 1.7 * phase)
    envelope = np.exp(-1.8 * radius**2)

    amplitude = 0.23
    displacement = np.column_stack(
        (
            amplitude * (0.75 * wave_one + 0.25 * envelope * radial_wave),
            amplitude * (0.55 * wave_two + 0.25 * envelope * radial_wave),
            amplitude * (0.50 * wave_one - 0.35 * wave_two + 0.45 * envelope * radial_wave),
        )
    )
    return displacement, wave_one + 0.65 * wave_two + envelope * radial_wave


def interpolate_colors(values):
    normalized = np.clip((values + 2.0) / 4.0, 0.0, 1.0)
    scaled = normalized * (len(COLOR_STOPS) - 1)
    lower = np.floor(scaled).astype(int)
    upper = np.minimum(lower + 1, len(COLOR_STOPS) - 1)
    fraction = (scaled - lower)[:, None]
    colors = COLOR_STOPS[lower] * (1.0 - fraction) + COLOR_STOPS[upper] * fraction
    return np.asarray(colors, dtype=np.uint8)


def project(positions, rotation, width, height, scale):
    rotated = positions @ rotation
    depth = rotated[:, 2]
    perspective = 1.0 + 0.0040 * (depth - depth.min())
    x = width * 0.5 + rotated[:, 0] * perspective * scale
    y = height * 0.5 - rotated[:, 1] * perspective * scale
    return x, y, depth


def draw_frame(reference, box_length, phase, rotation, width, height, scale):
    displacement, wave_value = wave_displacement(reference, box_length, phase)
    positions = reference + displacement
    x, y, depth = project(positions, rotation, width, height, scale)
    order = np.argsort(depth)
    x = x[order]
    y = y[order]
    depth = depth[order]
    wave_value = wave_value[order]
    normalized_depth = (depth - depth.min()) / (np.ptp(depth) + 1.0e-12)
    colors = interpolate_colors(wave_value + 0.55 * (normalized_depth - 0.5))
    radii = 1.15 + 1.55 * normalized_depth**1.6

    image = Image.new("RGB", (width, height), tuple(BACKGROUND.tolist()))
    draw = ImageDraw.Draw(image, "RGB")
    for atom_x, atom_y, radius, color in zip(x, y, radii, colors):
        if atom_x < -4 or atom_x > width + 4 or atom_y < -4 or atom_y > height + 4:
            continue
        fill = tuple(int(value) for value in color)
        draw.ellipse(
            (atom_x - radius, atom_y - radius, atom_x + radius, atom_y + radius),
            fill=fill,
        )
    return image


def adaptive_palette(frames):
    sample = frames[0].quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    palette = sample.getpalette()
    output = []
    for frame in frames:
        quantized = frame.quantize(palette=sample, dither=Image.Dither.NONE)
        quantized.putpalette(palette)
        output.append(quantized)
    return output


def make_animation(cells, lattice_parameter, frame_count, fps, width, height, output):
    reference, box_length = bcc_supercell(cells, lattice_parameter)
    rotation = rotation_matrix()
    corners = np.array(
        [
            [x, y, z]
            for x in (-0.5 * box_length, 0.5 * box_length)
            for y in (-0.5 * box_length, 0.5 * box_length)
            for z in (-0.5 * box_length, 0.5 * box_length)
        ]
    )
    rotated_corners = corners @ rotation
    scale = 0.82 * min(
        width / np.ptp(rotated_corners[:, 0]),
        height / np.ptp(rotated_corners[:, 1]),
    )

    frames = []
    for frame_index in range(frame_count):
        phase = 2.0 * np.pi * frame_index / frame_count
        frames.append(
            draw_frame(reference, box_length, phase, rotation, width, height, scale)
        )
        if (frame_index + 1) % 10 == 0 or frame_index + 1 == frame_count:
            print(f"Rendered {frame_index + 1}/{frame_count} frames")

    frames = adaptive_palette(frames)
    duration_ms = int(round(1000.0 / fps))
    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
        disposal=2,
    )
    return len(reference)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "iron_ml_large_scale_wave_dynamics.gif"
    number_atoms = make_animation(
        args.cells,
        args.lattice_parameter,
        args.frames,
        args.fps,
        args.width,
        args.height,
        output,
    )

    if args.poster_dir is not None:
        args.poster_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(output, args.poster_dir / output.name)

    print(f"Animated {number_atoms:,} BCC Fe atoms")
    print(f"GIF: {output}")


if __name__ == "__main__":
    main()
