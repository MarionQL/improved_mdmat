# plot_interface_contacts.py

A Python script for visualizing residue-residue interface contact frequency as a heatmap. It reads the `.npz` output files produced by `interface_mdmat.py` and generates a publication-ready figure showing which residue pairs spend the most time in contact at the interface.

It supports single runs and multi-replicate analysis — when given multiple files it can either average them into one heatmap or display them side by side as separate panels.

---

## Table of Contents

- [What the plot shows](#what-the-plot-shows)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [All Options](#all-options)
- [Working with Multiple Replicates](#working-with-multiple-replicates)
- [How the Threshold Works](#how-the-threshold-works)
- [Reading the Plot](#reading-the-plot)
- [Output Formats](#output-formats)
- [Common Recipes](#common-recipes)
- [Troubleshooting](#troubleshooting)

---

## What the plot shows

The heatmap is a grid where:

- **Rows** = residues from your first selection (Group A / x-axis group)
- **Columns** = residues from your second selection (Group B / y-axis group)
- **Color intensity** = how often that pair of residues was in contact, from 0% (white) to 100% (dark blue)

Only residues that exceed your chosen threshold in at least one replicate are shown — everything else is filtered out to keep the figure focused on the real interface.

Residue labels use **single-letter amino acid codes** (e.g. `E42` for GLU42, `K107` for LYS107) with an optional chain prefix (e.g. `A:E42`).

---

## Requirements

- Python 3.9 or newer
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)

Both are likely already installed if you followed the `interface_mdmat.py` setup. If not:

```bash
conda activate mdanalysis
conda install -c conda-forge matplotlib
# or
pip install matplotlib
```

---

## Installation

No installation needed — just place `plot_interface_contacts.py` in the same folder as your `.npz` files (or anywhere on your path) and run it with Python.

To confirm your environment is ready:

```bash
python -c "import numpy; import matplotlib; print('All good!')"
```

---

## Quick Start

**Single replicate, save to PNG:**

```bash
python plot_interface_contacts.py \
    rep1.npz \
    --threshold 0.3 \
    -o contact_map.png
```

**Single replicate with a title:**

```bash
python plot_interface_contacts.py \
    rep1.npz \
    --threshold 0.3 \
    --title "ELFN2–mGluR7 Interface (last 500 ns)" \
    -o contact_map.png
```

**Three replicates averaged together:**

```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.3 \
    --title "ELFN2–mGluR7 Interface" \
    -o contact_map_avg.png
```

**Preview without saving (interactive window):**

```bash
python plot_interface_contacts.py rep1.npz --threshold 0.3
```

---

## All Options

| Flag | Short | Required | Default | Description |
|---|---|---|---|---|
| `NPZ [NPZ ...]` | | ✅ | — | One or more `.npz` files from `interface_mdmat.py` |
| `--threshold` | | | `0.3` | Minimum contact frequency (0.0–1.0) to include a residue |
| `--title` | | | none | Title text shown at the top of the figure |
| `--output` | `-o` | | interactive | Output file path (see [Output Formats](#output-formats)) |
| `--dpi` | | | `150` | Resolution for PNG output. Use `300` for publication. |
| `--panel` | | | off | Show each replicate as its own panel instead of averaging |
| `--figwidth` | | | auto | Figure width in inches |
| `--figheight` | | | auto | Figure height in inches |
| `--fontsize` | | | `7.0` | Tick label size in points |
| `--no-chain-prefix` | | | off | Strip chain/segid prefix from labels (`A:E42` → `E42`) |

> **Note on `--threshold`:** The value is a fraction, not a percentage. Use `0.3` for 30%, `0.5` for 50%, etc.

---

## Working with Multiple Replicates

When you pass more than one `.npz` file, the script aligns the residue labels across all files before plotting. As long as your replicates come from the same system with the same selections, they will align automatically.

### Averaged heatmap (default)

Contact frequencies are averaged element-wise across replicates. The subtitle "Mean of N replicates" is added automatically.

```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.3 \
    -o averaged.png
```

### Side-by-side panel view (`--panel`)

Each replicate is shown as its own panel with a shared colorbar and shared y-axis scale, making it easy to compare consistency between runs. Each panel is labeled with the filename stem.

```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.3 \
    --panel \
    -o panel_view.png
```

### What happens if replicates don't match exactly?

The script takes the **intersection** of residue labels across all files, dropping any residues that don't appear in every replicate. A warning is printed if any are dropped:

```
[WARN] Alignment dropped 2 x-residues and 0 y-residues not shared across all replicates.
```

This is expected if your replicates used slightly different selections. If many residues are dropped, double-check that all files came from the same system and used the same `--x-select` / `--y-select` arguments in `interface_mdmat.py`.

---

## How the Threshold Works

The threshold filters out residue pairs that never form a meaningful contact. A residue is **kept** if it has at least one partner on the other side of the interface with a contact frequency at or above the threshold in **at least one replicate**.

For example, with `--threshold 0.3`:

- A residue pair that is in contact 35% of the time in rep1 but 10% in rep2 and rep3 → **shown** (exceeds threshold in at least one replicate)
- A residue pair that never exceeds 20% in any replicate → **hidden**

This approach keeps interface-relevant residues visible even when they're not consistently present in all replicates.

The threshold is also drawn as a **dashed red line** on the colorbar so you can always see where the cutoff falls on the color scale.

> **Tip:** Start with `--threshold 0.3` (30%) and adjust up or down depending on how crowded or sparse the result looks. If you see hundreds of residues, try `0.5`. If you see very few, try `0.1`.

---

## Reading the Plot

- **Dark blue cells** = residue pairs in contact most of the time (high frequency)
- **Light blue / white cells** = residues that were rarely or never in contact
- **Dashed red line on colorbar** = your chosen threshold
- **Stats line below the plot** = how many residue pairs are shown vs. the total number screened
- **Axis labels** use one-letter amino acid codes with residue number (e.g. `E42` = Glu42)
- **Chain prefix** (e.g. `A:E42`) is shown by default; use `--no-chain-prefix` to remove it

---

## Output Formats

The `-o` / `--output` flag accepts any file extension that Matplotlib supports:

| Extension | Format | Best for |
|---|---|---|
| `.png` | Raster image | Presentations, quick sharing |
| `.pdf` | Vector PDF | Publications, editable in Illustrator |
| `.svg` | Scalable vector | Web, editable in Inkscape |
| `.eps` | Encapsulated PostScript | Some journal submission systems |

For publication figures, use `.pdf` or `.svg` (vector), or `.png` with `--dpi 300`.

```bash
# High-resolution PNG
python plot_interface_contacts.py rep1.npz --threshold 0.3 --dpi 300 -o figure.png

# Editable PDF (text remains editable in Illustrator)
python plot_interface_contacts.py rep1.npz --threshold 0.3 -o figure.pdf
```

If `-o` is omitted entirely, an interactive Matplotlib window opens instead.

---

## Common Recipes

**Publication-quality averaged figure with no chain prefix:**
```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.3 \
    --title "ELFN2–mGluR7 Interface (last 500 ns)" \
    --no-chain-prefix \
    --dpi 300 \
    -o figure_3A.pdf
```

**Compare all three replicates side by side:**
```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.25 \
    --panel \
    --title "Replicate comparison" \
    -o panel_comparison.png
```

**Stricter threshold to highlight only the most persistent contacts:**
```bash
python plot_interface_contacts.py \
    rep1.npz rep2.npz rep3.npz \
    --threshold 0.5 \
    --title "Persistent contacts (≥50%)" \
    -o persistent_contacts.png
```

**Fix figure size if auto-sizing isn't right for your number of residues:**
```bash
python plot_interface_contacts.py \
    rep1.npz \
    --threshold 0.3 \
    --figwidth 14 \
    --figheight 10 \
    -o contact_map.png
```

**Increase font size if labels are hard to read:**
```bash
python plot_interface_contacts.py \
    rep1.npz \
    --threshold 0.3 \
    --fontsize 9 \
    -o contact_map.png
```

---

## Troubleshooting

**`ERROR: File not found: rep1.npz`**
Check that the path to your `.npz` file is correct. If the file is in a subdirectory, include the path: `results/rep1.npz`.

**`ERROR: No residues pass the threshold of 30%`**
Every residue pair had a contact frequency below your threshold. Either your threshold is too high, or your data genuinely has no persistent contacts. Try lowering the value:
```bash
--threshold 0.1
```
You can also inspect your raw data first:
```python
import numpy as np
data = np.load("rep1.npz", allow_pickle=True)
freq = data["contact_frequency"]
print(f"Max contact frequency in file: {freq.max():.1%}")
print(f"Mean contact frequency: {freq.mean():.1%}")
```

**`ERROR: After aligning replicates, no shared residue labels remain`**
Your `.npz` files don't share any residue labels in common. This usually means they came from different systems or used different `--x-select` / `--y-select` arguments. Check that all files were produced with the same selections.

**The figure is very wide / very tall**
Use `--figwidth` and `--figheight` to override the auto-sizing. If you have many residues passing the threshold, also consider raising `--threshold` to reduce clutter, or increasing `--fontsize` slightly and providing explicit dimensions.

**Labels are overlapping**
Try reducing `--fontsize` (e.g. `--fontsize 5`) or increasing the figure size with `--figwidth` and `--figheight`.

**`ImportError: No module named 'matplotlib'`**
Matplotlib isn't installed in the active environment. Run:
```bash
conda activate mdanalysis
conda install -c conda-forge matplotlib
```

**Panel figure is too narrow**
With many replicates and residues the auto-size can be tight. Override with:
```bash
--figwidth 24 --figheight 8
```
