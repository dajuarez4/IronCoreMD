"""Compute the hcp paramagnetic Birch-Murnaghan pressure-volume curve."""

from pathlib import Path


V0_BOHR3 = 70.1
B0_GPA = 217.4
B0_PRIME = 9.9
BOHR3_TO_ANGSTROM3 = 0.14818471147216278


def pressure_gpa(volume_bohr3: float) -> float:
    """Third-order Birch-Murnaghan pressure in GPa."""
    ratio = V0_BOHR3 / volume_bohr3
    return (
        1.5
        * B0_GPA
        * (ratio ** (7.0 / 3.0) - ratio ** (5.0 / 3.0))
        * (1.0 + 0.75 * (B0_PRIME - 4.0) * (ratio ** (2.0 / 3.0) - 1.0))
    )


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    volumes = [V0_BOHR3 * (0.70 + i * 0.005) for i in range(121)]
    pressures = [pressure_gpa(volume) for volume in volumes]

    csv_path = output_dir / "hcp_paramagnetic_p_vs_v.csv"
    rows = ["volume_bohr3,volume_angstrom3,pressure_GPa"]
    rows.extend(
        f"{v:.8f},{v * BOHR3_TO_ANGSTROM3:.8f},{p:.8f}"
        for v, p in zip(volumes, pressures)
    )
    csv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    # Make a portable SVG plot without third-party dependencies.
    width, height = 900, 600
    left, right, top, bottom = 100, 35, 55, 80
    xmin, xmax = min(volumes), max(volumes)
    ymin, ymax = min(pressures), max(pressures)
    xscale = lambda x: left + (x - xmin) / (xmax - xmin) * (width - left - right)
    yscale = lambda y: top + (ymax - y) / (ymax - ymin) * (height - top - bottom)
    points = " ".join(f"{xscale(v):.2f},{yscale(p):.2f}" for v, p in zip(volumes, pressures))
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="20">hcp paramagnetic: third-order Birch-Murnaghan EOS</text>',
    ]
    for i in range(7):
        x = xmin + i * (xmax - xmin) / 6
        sx = xscale(x)
        svg.extend([
            f'<line x1="{sx:.2f}" y1="{top}" x2="{sx:.2f}" y2="{height-bottom}" stroke="#ddd"/>',
            f'<text x="{sx:.2f}" y="{height-bottom+25}" text-anchor="middle" font-family="sans-serif" font-size="14">{x:.1f}</text>',
        ])
    for i in range(7):
        y = ymin + i * (ymax - ymin) / 6
        sy = yscale(y)
        svg.extend([
            f'<line x1="{left}" y1="{sy:.2f}" x2="{width-right}" y2="{sy:.2f}" stroke="#ddd"/>',
            f'<text x="{left-12}" y="{sy+5:.2f}" text-anchor="end" font-family="sans-serif" font-size="14">{y:.0f}</text>',
        ])
    svg.extend([
        f'<line x1="{left}" y1="{yscale(0):.2f}" x2="{width-right}" y2="{yscale(0):.2f}" stroke="black"/>',
        f'<line x1="{xscale(V0_BOHR3):.2f}" y1="{top}" x2="{xscale(V0_BOHR3):.2f}" y2="{height-bottom}" stroke="#d62728" stroke-dasharray="7 5"/>',
        f'<polyline points="{points}" fill="none" stroke="#1f77b4" stroke-width="3"/>',
        f'<text x="{width/2}" y="{height-18}" text-anchor="middle" font-family="sans-serif" font-size="17">Volume V (Bohr³)</text>',
        f'<text x="22" y="{height/2}" text-anchor="middle" font-family="sans-serif" font-size="17" transform="rotate(-90 22 {height/2})">Pressure P (GPa)</text>',
        f'<text x="{xscale(V0_BOHR3)+8:.2f}" y="{top+22}" font-family="sans-serif" font-size="14" fill="#d62728">V₀ = {V0_BOHR3:.1f} Bohr³</text>',
        '</svg>',
    ])
    (output_dir / "hcp_paramagnetic_p_vs_v.svg").write_text("\n".join(svg), encoding="utf-8")

    for i in range(13):
        fraction = 0.70 + 0.05 * i
        volume = fraction * V0_BOHR3
        volume_angstrom3 = volume * BOHR3_TO_ANGSTROM3
        print(
            f"{fraction:4.2f} V0  {volume:8.3f} Bohr^3  "
            f"{volume_angstrom3:8.3f} Angstrom^3  {pressure_gpa(volume):10.3f} GPa"
        )


if __name__ == "__main__":
    main()
