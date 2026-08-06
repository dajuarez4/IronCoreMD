from pathlib import Path
import re
import textwrap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path(__file__).resolve().parent
TITLE = ("Earth's Core: Exploring Iron Allotropes under Extreme Conditions through "
         "First-Principles Simulations and Graph-Kernel-Based Machine Learning")


def clean_latex(text):
    substitutions = {
        "--": "–",
        "\\&": "&",
        "\\textbar{}": "|",
        "\\,": " ",
    }
    for old, new in substitutions.items():
        text = text.replace(old, new)
    text = re.sub(r"\\textbf\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\emph\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\$([^$]*)\$", r"\1", text)
    text = text.replace("\\%", "%")
    return " ".join(text.split())


def extract_blocks(source):
    blocks = []
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", source, re.S)
    if abstract:
        blocks.append(("ABSTRACT", clean_latex(abstract.group(1))))
    body_start = source.find("\\section{Scientific Motivation}")
    body_end = source.find("\\end{document}")
    body = source[body_start:body_end]
    matches = list(re.finditer(r"\\section\*?\{([^}]*)\}", body))
    for index, match in enumerate(matches):
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(body)
        blocks.append((match.group(1).upper(), clean_latex(body[start:end])))
    return blocks


def wrapped_lines(text, width):
    paragraphs = text.split("\n")
    lines = []
    for paragraph in paragraphs:
        lines.extend(textwrap.wrap(paragraph, width=width, break_long_words=False))
        lines.append("")
    return lines


def main():
    source = (ROOT / "summary.tex").read_text()
    blocks = extract_blocks(source)
    with PdfPages(ROOT / "summary_preview.pdf") as pdf:
        fig = plt.figure(figsize=(8.5, 11), facecolor="white")
        fig.text(0.5, 0.965, "SC26 RESEARCH POSTERS — 800-WORD SUMMARY", ha="center",
                 va="top", color="#2CB1BC", fontsize=9, weight="bold")
        fig.text(0.5, 0.925, "\n".join(textwrap.wrap(TITLE, 70)), ha="center", va="top",
                 color="#102A43", fontsize=15, weight="bold", linespacing=1.12)
        fig.text(0.5, 0.855, "Diego Juarez  |  IronCoreMD Project", ha="center", va="top",
                 color="#243B53", fontsize=10)
        columns = [0.07, 0.525]
        column = 0
        y = 0.82
        line_height = 0.0108
        for heading, text in blocks:
            lines = wrapped_lines(text, 64)
            required = 0.025 + len(lines) * line_height
            if y - required < 0.055:
                column += 1
                if column == 2:
                    pdf.savefig(fig)
                    plt.close(fig)
                    fig = plt.figure(figsize=(8.5, 11), facecolor="white")
                    fig.text(0.5, 0.965, TITLE, ha="center", va="top",
                             color="#102A43", fontsize=9, weight="bold")
                    column = 0
                    y = 0.935
                else:
                    y = 0.82
            x = columns[column]
            fig.text(x, y, heading, ha="left", va="top", color="#102A43",
                     fontsize=9.4, weight="bold")
            y -= 0.018
            for line in lines:
                fig.text(x, y, line, ha="left", va="top", color="#111111", fontsize=8.3)
                y -= line_height
            y -= 0.006
        pdf.savefig(fig)
        plt.close(fig)


if __name__ == "__main__":
    main()
