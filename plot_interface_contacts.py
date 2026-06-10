#!/usr/bin/env python3
"""
plot_interface_contacts.py — Heatmap visualization for interface_mdmat.py output.

Plots residue-residue contact frequency across one or more replicate NPZ files.
Residues are filtered to show only those that exceed a user-specified contact
frequency threshold in at least one replicate. Single-letter amino acid codes
are used on the axes.

Usage examples:
    # Single replicate
    python plot_interface_contacts.py \
        rep1.npz \
        --threshold 0.3 \
        --title "ELFN2–mGluR7 Interface (last 500 ns)" \
        -o interface_contacts.png

    # Three replicates — averaged
    python plot_interface_contacts.py \
        rep1.npz rep2.npz rep3.npz \
        --threshold 0.25 \
        --title "ELFN2–mGluR7 Interface (3 replicates)" \
        -o interface_contacts_avg.png

    # Show individual replicate panels side by side
    python plot_interface_contacts.py \
        rep1.npz rep2.npz rep3.npz \
        --threshold 0.25 \
        --panel \
        -o interface_contacts_panel.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize

matplotlib.rcParams["pdf.fonttype"] = 42   # editable text in PDFs
matplotlib.rcParams["ps.fonttype"] = 42


# ---------------------------------------------------------------------------
# Amino acid three-letter → one-letter lookup
# ---------------------------------------------------------------------------

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common non-standard / protonation variants
    "HIE": "H", "HID": "H", "HIP": "H",
    "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E",
    "LYN": "K",
    "MSE": "M",
}


def convert_label(label: str) -> str:
    """
    Convert a residue label like 'A:GLU42' or 'GLU42' to 'E42'.
    Preserves the chain/segid prefix as 'A:E42' if present.
    Falls back gracefully for unknown residue names.
    """
    prefix = ""
    core = label

    if ":" in label:
        prefix, core = label.split(":", 1)
        prefix = prefix + ":"

    # split residue name from residue number, e.g. 'GLU42' → ('GLU', '42')
    i = 0
    while i < len(core) and core[i].isalpha():
        i += 1
    resname = core[:i].upper()
    resnum = core[i:]

    one_letter = THREE_TO_ONE.get(resname, resname[0] if resname else "?")
    return f"{prefix}{one_letter}{resnum}"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot residue-residue interface contact frequency heatmaps "
            "from one or more interface_mdmat.py NPZ output files."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="NPZ",
        help="One or more .npz files produced by interface_mdmat.py.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.3,
        help=(
            "Minimum contact frequency (0.0–1.0) a residue pair must reach "
            "in at least one replicate to be shown. E.g. 0.3 = 30%%."
        ),
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Title to display on the figure. Omit to show no title.",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Output image path. Supported formats: png, pdf, svg, eps. "
            "If omitted, the figure is shown interactively."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolution in dots per inch for raster output (png).",
    )
    parser.add_argument(
        "--panel",
        action="store_true",
        help=(
            "When multiple files are provided, plot each replicate as a "
            "separate panel instead of averaging them together."
        ),
    )
    parser.add_argument(
        "--figwidth",
        type=float,
        default=None,
        help="Figure width in inches. Auto-sized if omitted.",
    )
    parser.add_argument(
        "--figheight",
        type=float,
        default=None,
        help="Figure height in inches. Auto-sized if omitted.",
    )
    parser.add_argument(
        "--fontsize",
        type=float,
        default=7.0,
        help="Tick label font size in points.",
    )
    parser.add_argument(
        "--no-chain-prefix",
        action="store_true",
        help="Strip chain/segid prefix from tick labels (e.g. 'A:E42' → 'E42').",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_npz(path: str) -> dict:
    """Load and validate an interface_mdmat NPZ file."""
    try:
        data = np.load(path, allow_pickle=True)
    except FileNotFoundError:
        raise SystemExit(f"ERROR: File not found: {path}")
    except Exception as exc:
        raise SystemExit(f"ERROR: Could not read {path}: {exc}")

    required = {"contact_frequency", "x_labels", "y_labels"}
    missing = required - set(data.files)
    if missing:
        raise SystemExit(
            f"ERROR: {path} is missing expected keys: {missing}\n"
            "Make sure this file was produced by interface_mdmat.py."
        )

    return {
        "contact_frequency": data["contact_frequency"],           # (n_x, n_y)
        "x_labels":          data["x_labels"].tolist(),
        "y_labels":          data["y_labels"].tolist(),
        "path":              path,
    }


def align_replicates(datasets: list[dict]) -> tuple[list[dict], list[str], list[str]]:
    """
    Align multiple replicates to a common set of residue labels.

    Replicates are expected to share the same residue order, but this function
    handles minor mismatches by intersecting labels and reindexing matrices.
    Returns reindexed datasets and the shared x/y label lists.
    """
    if len(datasets) == 1:
        d = datasets[0]
        return datasets, d["x_labels"], d["y_labels"]

    # Determine shared labels (preserving order from first replicate)
    ref_x = datasets[0]["x_labels"]
    ref_y = datasets[0]["y_labels"]

    for d in datasets[1:]:
        # Find labels present in both, keeping ref order
        set_x = set(d["x_labels"])
        set_y = set(d["y_labels"])
        ref_x = [l for l in ref_x if l in set_x]
        ref_y = [l for l in ref_y if l in set_y]

    if not ref_x or not ref_y:
        raise SystemExit(
            "ERROR: After aligning replicates, no shared residue labels remain. "
            "Check that all input files come from the same system."
        )

    dropped_x = len(datasets[0]["x_labels"]) - len(ref_x)
    dropped_y = len(datasets[0]["y_labels"]) - len(ref_y)
    if dropped_x or dropped_y:
        print(
            f"[WARN] Alignment dropped {dropped_x} x-residues and "
            f"{dropped_y} y-residues not shared across all replicates."
        )

    aligned = []
    for d in datasets:
        xi = [d["x_labels"].index(l) for l in ref_x]
        yi = [d["y_labels"].index(l) for l in ref_y]
        freq = d["contact_frequency"][np.ix_(xi, yi)]
        aligned.append({**d, "contact_frequency": freq, "x_labels": ref_x, "y_labels": ref_y})

    return aligned, ref_x, ref_y


def apply_threshold(
    matrices: list[np.ndarray],
    threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return boolean index arrays for x and y residues that exceed `threshold`
    in at least one replicate matrix.
    """
    stacked = np.stack(matrices, axis=0)          # (n_reps, n_x, n_y)
    max_over_reps = stacked.max(axis=0)           # (n_x, n_y)

    x_mask = max_over_reps.max(axis=1) >= threshold   # any y partner exceeds threshold
    y_mask = max_over_reps.max(axis=0) >= threshold

    return x_mask, y_mask


# ---------------------------------------------------------------------------
# Label formatting
# ---------------------------------------------------------------------------

def format_labels(labels: list[str], strip_prefix: bool) -> list[str]:
    converted = [convert_label(l) for l in labels]
    if strip_prefix:
        converted = [l.split(":", 1)[-1] if ":" in l else l for l in converted]
    return converted


# ---------------------------------------------------------------------------
# Figure sizing helpers
# ---------------------------------------------------------------------------

def auto_figsize(
    n_x: int,
    n_y: int,
    fontsize: float,
    n_panels: int = 1,
) -> tuple[float, float]:
    """
    Estimate a figure size that keeps tick labels readable.
    Each residue gets ~0.22 inches; panels are laid out in a row.
    """
    points_per_inch = 72.0
    char_width_in = fontsize / points_per_inch * 3.5   # rough char width
    cell_size = max(0.22, char_width_in)

    base_w = n_y * cell_size + 3.0   # +3 for colorbar, margins, y-axis labels
    base_h = n_x * cell_size + 2.5   # +2.5 for x-axis labels, title

    w = base_w * n_panels + (n_panels - 1) * 0.5
    return (w, base_h)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _draw_heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    x_tick_labels: list[str],
    y_tick_labels: list[str],
    fontsize: float,
    vmin: float = 0.0,
    vmax: float = 1.0,
    cmap: str = "Blues",
    norm: Normalize | None = None,
) -> matplotlib.image.AxesImage:
    im = ax.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        norm=norm if norm is not None else Normalize(vmin=vmin, vmax=vmax),
        interpolation="nearest",
        origin="upper",
    )

    ax.set_xticks(np.arange(len(y_tick_labels)))
    ax.set_xticklabels(y_tick_labels, rotation=90, fontsize=fontsize, fontfamily="monospace")

    ax.set_yticks(np.arange(len(x_tick_labels)))
    ax.set_yticklabels(x_tick_labels, fontsize=fontsize, fontfamily="monospace")

    ax.tick_params(axis="both", which="both", length=0)

    # Light grid lines between cells
    ax.set_xticks(np.arange(-0.5, len(y_tick_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(x_tick_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im


def plot_averaged(
    datasets: list[dict],
    x_mask: np.ndarray,
    y_mask: np.ndarray,
    x_labels_fmt: list[str],
    y_labels_fmt: list[str],
    args: argparse.Namespace,
) -> plt.Figure:
    matrices = [d["contact_frequency"] for d in datasets]

    if len(matrices) > 1:
        avg_matrix = np.mean(np.stack(matrices, axis=0), axis=0)
        subtitle = f"Mean of {len(matrices)} replicates"
    else:
        avg_matrix = matrices[0]
        subtitle = None

    sub = avg_matrix[np.ix_(x_mask, y_mask)]
    x_ticks = [x_labels_fmt[i] for i, v in enumerate(x_mask) if v]
    y_ticks = [y_labels_fmt[i] for i, v in enumerate(y_mask) if v]

    n_x, n_y = sub.shape

    if args.figwidth and args.figheight:
        figsize = (args.figwidth, args.figheight)
    else:
        figsize = auto_figsize(n_x, n_y, args.fontsize)

    fig, ax = plt.subplots(figsize=figsize)

    norm = Normalize(vmin=0.0, vmax=1.0)
    im = _draw_heatmap(ax, sub, x_ticks, y_ticks, args.fontsize, norm=norm)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, aspect=30)
    cbar.set_label("Contact frequency", fontsize=args.fontsize + 1)
    cbar.ax.tick_params(labelsize=args.fontsize)

    # Threshold line on colorbar
    cbar.ax.axhline(args.threshold, color="crimson", linewidth=1.2, linestyle="--")
    cbar.ax.text(
        1.6, args.threshold, f" {args.threshold:.0%}",
        va="center", ha="left",
        fontsize=args.fontsize - 0.5, color="crimson",
        transform=cbar.ax.transData,
    )

    ax.set_xlabel("Group B residues", fontsize=args.fontsize + 1, labelpad=4)
    ax.set_ylabel("Group A residues", fontsize=args.fontsize + 1, labelpad=4)

    title_parts = []
    if args.title:
        title_parts.append(args.title)
    if subtitle:
        title_parts.append(subtitle)
    if title_parts:
        fig.suptitle("\n".join(title_parts), fontsize=args.fontsize + 3, y=1.01)

    _add_stats_text(fig, sub, args.threshold, args.fontsize)

    fig.tight_layout()
    return fig


def plot_panel(
    datasets: list[dict],
    x_mask: np.ndarray,
    y_mask: np.ndarray,
    x_labels_fmt: list[str],
    y_labels_fmt: list[str],
    args: argparse.Namespace,
) -> plt.Figure:
    n_reps = len(datasets)

    x_ticks = [x_labels_fmt[i] for i, v in enumerate(x_mask) if v]
    y_ticks = [y_labels_fmt[i] for i, v in enumerate(y_mask) if v]
    n_x = sum(x_mask)
    n_y = sum(y_mask)

    if args.figwidth and args.figheight:
        figsize = (args.figwidth, args.figheight)
    else:
        figsize = auto_figsize(n_x, n_y, args.fontsize, n_panels=n_reps)

    fig, axes = plt.subplots(1, n_reps, figsize=figsize, sharey=True)
    if n_reps == 1:
        axes = [axes]

    norm = Normalize(vmin=0.0, vmax=1.0)

    ims = []
    for ax, d in zip(axes, datasets):
        sub = d["contact_frequency"][np.ix_(x_mask, y_mask)]
        rep_name = Path(d["path"]).stem
        im = _draw_heatmap(ax, sub, x_ticks, y_ticks, args.fontsize, norm=norm)
        ax.set_title(rep_name, fontsize=args.fontsize + 1, pad=4)
        ax.set_xlabel("Group B residues", fontsize=args.fontsize + 1, labelpad=4)
        ims.append(im)

    axes[0].set_ylabel("Group A residues", fontsize=args.fontsize + 1, labelpad=4)

    # Shared colorbar on the right
    cbar = fig.colorbar(ims[-1], ax=axes, fraction=0.02, pad=0.02, aspect=30)
    cbar.set_label("Contact frequency", fontsize=args.fontsize + 1)
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.ax.axhline(args.threshold, color="crimson", linewidth=1.2, linestyle="--")
    cbar.ax.text(
        1.6, args.threshold, f" {args.threshold:.0%}",
        va="center", ha="left",
        fontsize=args.fontsize - 0.5, color="crimson",
        transform=cbar.ax.transData,
    )

    if args.title:
        fig.suptitle(args.title, fontsize=args.fontsize + 3, y=1.01)

    fig.tight_layout()
    return fig


def _add_stats_text(fig: plt.Figure, matrix: np.ndarray, threshold: float, fontsize: float) -> None:
    """Print a small summary below the figure."""
    n_pairs = matrix.size
    n_above = int((matrix >= threshold).sum())
    pct = 100.0 * n_above / n_pairs if n_pairs else 0.0
    fig.text(
        0.5, -0.01,
        f"Showing {n_above} / {n_pairs} residue pairs ≥ {threshold:.0%} contact frequency",
        ha="center", va="top",
        fontsize=fontsize - 0.5,
        color="#555555",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if not 0.0 <= args.threshold <= 1.0:
        raise SystemExit("ERROR: --threshold must be between 0.0 and 1.0")

    # Load all input files
    datasets = [load_npz(p) for p in args.inputs]
    print(f"[INFO] Loaded {len(datasets)} file(s).")

    # Align residue labels across replicates
    datasets, shared_x_labels, shared_y_labels = align_replicates(datasets)

    n_x = len(shared_x_labels)
    n_y = len(shared_y_labels)
    print(f"[INFO] Shared residues: {n_x} (x-axis) × {n_y} (y-axis)")

    # Apply threshold filter
    matrices = [d["contact_frequency"] for d in datasets]
    x_mask, y_mask = apply_threshold(matrices, args.threshold)

    n_x_shown = int(x_mask.sum())
    n_y_shown = int(y_mask.sum())

    if n_x_shown == 0 or n_y_shown == 0:
        raise SystemExit(
            f"ERROR: No residues pass the threshold of {args.threshold:.0%}. "
            "Try lowering --threshold."
        )

    print(
        f"[INFO] Residues above threshold ({args.threshold:.0%}): "
        f"{n_x_shown} / {n_x} (x), {n_y_shown} / {n_y} (y)"
    )

    # Format labels
    x_labels_fmt = format_labels(shared_x_labels, args.no_chain_prefix)
    y_labels_fmt = format_labels(shared_y_labels, args.no_chain_prefix)

    # Draw figure
    if args.panel and len(datasets) > 1:
        fig = plot_panel(datasets, x_mask, y_mask, x_labels_fmt, y_labels_fmt, args)
    else:
        if args.panel and len(datasets) == 1:
            print("[INFO] Only one input file; --panel has no effect.")
        fig = plot_averaged(datasets, x_mask, y_mask, x_labels_fmt, y_labels_fmt, args)

    # Save or show
    if args.output:
        out_path = Path(args.output)
        fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")
        print(f"[DONE] Saved figure to {out_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
