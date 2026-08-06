from pathlib import Path
import textwrap

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
ORANGE = "#F58025"
DARK = "#562B12"
SOFT = "#FFF7EF"
BLUE = "#004E7A"
GRAY = "#404040"
PALE_BLUE = "#EEF6FA"


def panel(fig, x, y, width, height, title):
    outer = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.004,rounding_size=0.008",
        transform=fig.transFigure, facecolor="white", edgecolor=ORANGE,
        linewidth=2.3, zorder=-10,
    )
    fig.patches.append(outer)
    header_height = 0.038
    fig.patches.append(FancyBboxPatch(
        (x, y + height - header_height), width, header_height,
        boxstyle="round,pad=0.004,rounding_size=0.007",
        transform=fig.transFigure, facecolor=DARK, edgecolor=DARK,
        linewidth=0, zorder=-5,
    ))
    fig.text(x + 0.010, y + height - header_height / 2, title,
             color="white", fontsize=22, weight="bold", va="center")
    return x + 0.012, y + 0.014, width - 0.024, height - header_height - 0.022


def wrapped(text, width):
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def bullet_block(fig, x, y, items, width=54, fontsize=17, color=GRAY, spacing=0.042):
    cursor = y
    for item in items:
        lines = textwrap.wrap(item, width=width, break_long_words=False)
        fig.text(x, cursor, "•", fontsize=fontsize + 2, color=ORANGE, weight="bold", va="top")
        fig.text(x + 0.012, cursor, "\n".join(lines), fontsize=fontsize, color=color,
                 va="top", linespacing=1.22)
        cursor -= spacing * max(1, len(lines))
    return cursor


def add_image(fig, bounds, path, pad=0.004):
    x, y, width, height = bounds
    ax = fig.add_axes([x + pad, y + pad, width - 2 * pad, height - 2 * pad], zorder=3)
    ax.imshow(Image.open(path))
    ax.axis("off")
    return ax


def metric(fig, x, y, width, value, label):
    fig.patches.append(FancyBboxPatch(
        (x, y), width, 0.064,
        boxstyle="round,pad=0.003,rounding_size=0.006",
        transform=fig.transFigure, facecolor=SOFT, edgecolor=ORANGE,
        linewidth=1.7, zorder=-5,
    ))
    fig.text(x + width / 2, y + 0.041, value, ha="center", va="center",
             color=DARK, fontsize=25, weight="bold")
    fig.text(x + width / 2, y + 0.015, label, ha="center", va="center",
             color=GRAY, fontsize=12.5)


def workflow(fig, bounds):
    x, y, width, height = bounds
    labels_top = [
        "BCC | FCC | HCP\nP–T–magnetic states",
        "Quantum ESPRESSO\nDFT + AIMD",
        "Validated NPZ\nenergies + forces",
        "TDEP\nphonons + free energy",
    ]
    labels_bottom = [
        "Gaussian-process\ninteratomic model",
        "Graph kernel\nKᵢⱼ = k(Gᵢ,Gⱼ)",
        "Atomic graphs\nnodes + edges + labels",
    ]
    box_width = width * 0.205
    box_height = height * 0.32
    gap = (width - 4 * box_width) / 3
    top_y = y + height * 0.56
    centers = []
    for index, label in enumerate(labels_top):
        left = x + index * (box_width + gap)
        fig.patches.append(FancyBboxPatch(
            (left, top_y), box_width, box_height,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            transform=fig.transFigure, facecolor=SOFT, edgecolor=ORANGE,
            linewidth=1.8, zorder=-5,
        ))
        fig.text(left + box_width / 2, top_y + box_height / 2, label,
                 ha="center", va="center", fontsize=15.5, color=DARK, weight="bold")
        centers.append((left + box_width / 2, top_y + box_height / 2))
        if index:
            previous = centers[index - 1]
            fig.patches.append(FancyArrowPatch(
                (previous[0] + box_width / 2, previous[1]),
                (centers[index][0] - box_width / 2, centers[index][1]),
                transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
                linewidth=2, color=BLUE,
            ))
    bottom_y = y + height * 0.05
    start_x = x + width * 0.095
    bottom_gap = width * 0.085
    bottom_centers = []
    for index, label in enumerate(labels_bottom):
        left = start_x + index * (box_width + bottom_gap)
        fig.patches.append(FancyBboxPatch(
            (left, bottom_y), box_width, box_height,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            transform=fig.transFigure, facecolor=PALE_BLUE, edgecolor=BLUE,
            linewidth=1.8, zorder=-5,
        ))
        fig.text(left + box_width / 2, bottom_y + box_height / 2, label,
                 ha="center", va="center", fontsize=15.5, color=BLUE, weight="bold")
        bottom_centers.append((left + box_width / 2, bottom_y + box_height / 2))
        if index:
            fig.patches.append(FancyArrowPatch(
                (bottom_centers[index][0] - box_width / 2, bottom_centers[index][1]),
                (bottom_centers[index - 1][0] + box_width / 2, bottom_centers[index - 1][1]),
                transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
                linewidth=2, color=BLUE,
            ))
    fig.patches.append(FancyArrowPatch(
        (centers[2][0], top_y), (bottom_centers[2][0], bottom_y + box_height),
        transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
        linewidth=2, color=BLUE,
    ))


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(48, 36), facecolor="#FCFCFC")

    fig.patches.append(Rectangle((0, 0.885), 1, 0.115, transform=fig.transFigure,
                                 facecolor=SOFT, edgecolor="none", zorder=-10))
    fig.patches.append(Rectangle((0, 0.885), 1, 0.008, transform=fig.transFigure,
                                 facecolor=ORANGE, edgecolor="none", zorder=-9))

    logo = FIGURES / "utep+loog.png"
    if logo.exists():
        add_image(fig, (0.018, 0.905, 0.115, 0.075), logo)
    else:
        fig.text(0.075, 0.944, "UTEP\nLOGO", ha="center", va="center",
                 fontsize=28, weight="bold", color=DARK,
                 bbox={"boxstyle": "round,pad=0.55", "facecolor": "white",
                       "edgecolor": ORANGE, "linewidth": 2.2})

    title = ("Earth's Core: Exploring Iron Allotropes under Extreme Conditions through\n"
             "First-Principles Simulations and Graph-Kernel-Based Machine Learning")
    fig.text(0.51, 0.956, title, ha="center", va="center", fontsize=37,
             color=DARK, weight="bold", linespacing=1.08)
    fig.text(0.51, 0.910,
             "Diego Juárez  |  Department of Physics, The University of Texas at El Paso  |  IronCoreMD",
             ha="center", va="center", fontsize=20, color=BLUE)
    fig.text(0.945, 0.944, "SC26\nRESEARCH\nPOSTERS", ha="center", va="center",
             fontsize=19, color=DARK, weight="bold",
             bbox={"boxstyle": "round,pad=0.55", "facecolor": "white",
                   "edgecolor": ORANGE, "linewidth": 2.2})

    left_x, center_x, right_x = 0.018, 0.342, 0.666
    column_width = 0.316

    panel(fig, left_x, 0.485, column_width, 0.385, "INTRODUCTION")
    panel(fig, left_x, 0.060, column_width, 0.410, "MOTIVATION")

    bx, by, bw, bh = panel(fig, center_x, 0.585, column_width, 0.285, "METHODOLOGY")
    workflow_labels = [
        "BCC | FCC | HCP state points",
        "Quantum ESPRESSO: DFT + AIMD",
        "Validated NPZ dataset",
        "Atomic graphs + graph kernel",
        "Gaussian-process model",
        "EOS | TDEP | phonons | free energies",
    ]
    box_height = bh * 0.115
    gap = bh * 0.032
    top = by + bh - box_height
    for index, label in enumerate(workflow_labels):
        box_y = top - index * (box_height + gap)
        color = BLUE if index == len(workflow_labels) - 1 else ORANGE
        fill = PALE_BLUE if index == len(workflow_labels) - 1 else SOFT
        fig.patches.append(FancyBboxPatch(
            (bx + bw * 0.08, box_y), bw * 0.84, box_height,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            transform=fig.transFigure, facecolor=fill, edgecolor=color,
            linewidth=1.8, zorder=-5,
        ))
        fig.text(bx + bw * 0.50, box_y + box_height / 2, label,
                 ha="center", va="center", fontsize=15.5,
                 color=DARK if index < len(workflow_labels) - 1 else BLUE,
                 weight="bold")
        if index:
            previous_y = top - (index - 1) * (box_height + gap)
            fig.patches.append(FancyArrowPatch(
                (bx + bw * 0.50, previous_y),
                (bx + bw * 0.50, box_y + box_height),
                transform=fig.transFigure, arrowstyle="-|>", mutation_scale=16,
                linewidth=1.8, color=BLUE,
            ))

    rx, ry, rw, rh = panel(fig, center_x, 0.060, column_width, 0.510, "RESULTS")
    add_image(fig, (rx + rw * 0.06, ry + rh * 0.64, rw * 0.88, rh * 0.33),
              FIGURES / "graph_kernel_parity.png")
    metric(fig, rx, ry + rh * 0.48, rw * 0.30, "7.26", "MAE (meV/atom)")
    metric(fig, rx + rw * 0.35, ry + rh * 0.48, rw * 0.30, "9.61", "RMSE (meV/atom)")
    metric(fig, rx + rw * 0.70, ry + rh * 0.48, rw * 0.30, "156", "evaluation structures")
    add_image(fig, (rx + rw * 0.02, ry + rh * 0.05, rw * 0.96, rh * 0.39),
              FIGURES / "bcc_magnetic_phonons.png")

    bx, by, bw, bh = panel(fig, right_x, 0.635, column_width, 0.235, "DISCUSSION")
    bullet_block(fig, bx, by + bh - 0.006, [
        "Graph kernels provide a flexible representation across structurally diverse iron phases.",
        "Energy accuracy is encouraging, but force, EOS, phonon, and free-energy validation remain essential.",
        "Magnetic preparation changes individual branches and requires systematic convergence studies.",
    ], width=72, fontsize=16.5, spacing=0.037)

    bx, by, bw, bh = panel(fig, right_x, 0.335, column_width, 0.285, "FUTURE WORK")
    bullet_block(fig, bx, by + bh - 0.006, [
        "Train force components with trajectory-grouped splits.",
        "Expand magnetic sampling across volume and temperature.",
        "Compare BCC, FCC, and HCP free energies at matched conditions.",
        "Benchmark thousand-atom ML simulations and release reproducibility artifacts.",
    ], width=72, fontsize=16.5, spacing=0.035)
    add_image(fig, (bx + bw * 0.12, by + 0.010, bw * 0.76, bh * 0.42),
              FIGURES / "iron_ml_large_scale_simulation_space.png")

    bx, by, bw, bh = panel(fig, right_x, 0.225, column_width, 0.105, "ACKNOWLEDGMENTS")
    fig.text(bx, by + bh - 0.006,
             "Computational resources, collaborators, funding sources, and award numbers will be added before final submission.",
             va="top", fontsize=16.5, color=GRAY, wrap=True)

    bx, by, bw, bh = panel(fig, right_x, 0.060, column_width, 0.150, "REFERENCES")
    references = [
        "[1] Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009).",
        "[2] Hellman et al., Phys. Rev. B 84, 180301 (2011).",
        "[3] Vishwanathan et al., JMLR 11, 1201–1242 (2010).",
        "[4] Rasmussen and Williams, Gaussian Processes for Machine Learning (2006).",
    ]
    fig.text(bx, by + bh - 0.006, "\n\n".join(references),
             va="top", fontsize=14.5, color=GRAY, linespacing=1.15)

    fig.savefig(ROOT / "poster_preview.pdf")
    fig.savefig(ROOT / "poster_preview.png", dpi=75)
    plt.close(fig)


if __name__ == "__main__":
    main()
