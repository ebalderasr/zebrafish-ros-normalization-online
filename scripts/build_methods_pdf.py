#!/usr/bin/env python3
"""Typeset METHODS.md into docs/methods.pdf.

Requires pandoc and xelatex. Both ship with TeX Live::

    sudo apt install pandoc texlive-xetex texlive-fonts-recommended

Markdown written for GitHub uses HTML subscript tags, which LaTeX output drops.
They are rewritten to pandoc's own subscript syntax before the conversion.

Usage::

    python scripts/build_methods_pdf.py
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PREAMBLE = r"""
\usepackage{microtype}
\setlength{\emergencystretch}{2em}
\usepackage{titlesec}
\titleformat{\section}{\normalfont\Large\bfseries}{}{0pt}{}
\titleformat{\subsection}{\normalfont\large\bfseries}{}{0pt}{}
\titlespacing*{\section}{0pt}{2.2em}{0.8em}
\titlespacing*{\subsection}{0pt}{1.6em}{0.5em}
\usepackage{enumitem}
\setlist{itemsep=2pt, topsep=4pt}
"""


def to_pandoc_markdown(text: str) -> str:
    """Rewrite the HTML bits that only exist for the GitHub renderer."""
    text = re.sub(r"<sub>(.*?)</sub>", r"~\1~", text, flags=re.S)
    text = re.sub(r"<sup>(.*?)</sup>", r"^\1^", text, flags=re.S)
    return text


def build(source: Path, output: Path) -> int:
    for tool in ("pandoc", "xelatex"):
        if shutil.which(tool) is None:
            print(f"{tool} not found on PATH.", file=sys.stderr)
            return 1

    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        prepared = Path(tmp) / "methods.md"
        prepared.write_text(to_pandoc_markdown(source.read_text(encoding="utf-8")),
                            encoding="utf-8")
        preamble = Path(tmp) / "preamble.tex"
        preamble.write_text(PREAMBLE, encoding="utf-8")

        command = [
            "pandoc", str(prepared), "-o", str(output),
            "--pdf-engine=xelatex",
            "--from=markdown+tex_math_dollars+raw_tex",
            "--metadata", "title=Zebrafish DCF ROS Normalization Methods",
            "--metadata", "subtitle=Within-session normalization, outlier handling, and nested statistical analysis",
            "--metadata", "author=zebrafish-ros v1.0.0 · github.com/ebalderasr/zebrafish-ros-normalization-online",
            "--variable", "mainfont=Latin Modern Roman",
            "--variable", "sansfont=Latin Modern Sans",
            "--variable", "monofont=Latin Modern Mono",
            "--variable", "fontsize=11pt",
            "--variable", "geometry:a4paper,margin=2.6cm",
            "--variable", "linkcolor=black",
            "--variable", "urlcolor=black",
            "--variable", "colorlinks=true",
            "--include-in-header", str(preamble),
            "--toc", "--toc-depth=2",
        ]
        result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return result.returncode

    print(f"written: {output} ({output.stat().st_size // 1024} KB)")
    return 0


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, default=root / "METHODS.md")
    parser.add_argument("--output", type=Path, default=root / "docs" / "methods.pdf")
    args = parser.parse_args(argv)
    return build(args.source, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
