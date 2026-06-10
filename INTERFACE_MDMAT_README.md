# interface_mdmat.py

A Python script for computing **residue-residue interface contact matrices** between two groups of atoms across a molecular dynamics (MD) trajectory. Think of it as a between-group version of GROMACS's `gmx mdmat` tool.

Given two selections (e.g. chain A and chain B of a protein complex), it tells you:

- How close each residue in group A gets to each residue in group B on average
- The closest any pair of residues ever got across the whole trajectory
- How often each residue pair was in contact
- A snapshot of which residues were in contact in the first analyzed frame

---

## Table of Contents

- [What does "contact" mean?](#what-does-contact-mean)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [All Options](#all-options)
- [Selection Syntax](#selection-syntax)
- [Output Files](#output-files)
- [Loading Results in Python](#loading-results-in-python)
- [Common Recipes](#common-recipes)
- [Troubleshooting](#troubleshooting)

---

## What does "contact" mean?

Two residues are considered **in contact** when the closest distance between *any atom* in residue A and *any atom* in residue B is less than or equal to the cutoff distance (default: 4 Å). This is the standard heavy-atom contact definition used in structural biology.

---

## Requirements

- Python 3.9 or newer
- [MDAnalysis](https://www.mdanalysis.org/) — for reading trajectories
- [NumPy](https://numpy.org/) — for numerical arrays
- [pandas](https://pandas.pydata.org/) — for CSV output

---

## Installation

**Option 1 — conda (recommended)**

```bash
conda create -n mdanalysis python=3.11
conda activate mdanalysis
conda install -c conda-forge mdanalysis pandas
```

**Option 2 — pip**

```bash
pip install MDAnalysis pandas numpy
```

To verify everything is installed correctly:

```bash
python -c "import MDAnalysis; import pandas; import numpy; print('All good!')"
```

---

## Quick Start

The minimum you need to provide is a topology file, a trajectory file, the two atom selections you want to compare, and an output name:

```bash
python interface_mdmat.py \
    -s my_system.tpr \
    -f my_trajectory.xtc \
    --x-select "segid A and protein" \
    --y-select "segid B and protein" \
    -o results/my_analysis
```

This will create `results/my_analysis.npz` (the main data file), plus a metadata JSON and a labels CSV alongside it.

> **Tip:** The `-o` argument is a *prefix*, not a filename. The script appends suffixes automatically. If the `results/` folder doesn't exist, create it first with `mkdir results`.

---

## All Options

| Flag | Short | Required | Default | Description |
|---|---|---|---|---|
| `--topology` | `-s` | ✅ | — | Topology file (`.tpr`, `.pdb`, `.gro`, `.psf`) |
| `--trajectory` | `-f` | ✅ | — | Trajectory file (`.xtc`, `.trr`, `.dcd`) |
| `--output-prefix` | `-o` | ✅ | — | Prefix for all output files |
| `--x-select` | | ✅ | — | MDAnalysis selection for the x-axis group |
| `--y-select` | | ✅ | — | MDAnalysis selection for the y-axis group |
| `--cutoff` | | | `4.0` | Contact distance cutoff |
| `--cutoff-unit` | | | `angstrom` | Unit for the cutoff (`angstrom` or `nm`) |
| `--units` | | | `angstrom` | Unit for distances in output matrices |
| `--begin` | | | first frame | Start time in picoseconds |
| `--end` | | | last frame | End time in picoseconds |
| `--stride` | | | `1` | Analyze every Nth frame (e.g. `10` = every 10th) |
| `--backend` | | | `OpenMP` | Distance calculation backend (`OpenMP` or `serial`) |
| `--csv` | | | off | Also write CSV versions of all matrices |

### Time and frame filtering

`--begin` and `--end` use **picoseconds** (ps). Common conversions:

| You want | Use |
|---|---|
| Start at 500 ns | `--begin 500000` |
| Analyze only the last 100 ns of a 1 µs run | `--begin 900000` |
| Stop at 200 ns | `--end 200000` |

`--stride` skips frames relative to your time window, not the whole trajectory. `--stride 10` with `--begin 500000` will give you every 10th frame *starting from* 500 ns.

---

## Selection Syntax

Selections use [MDAnalysis selection language](https://docs.mdanalysis.org/stable/documentation_pages/selections.html), which is similar to VMD. Always wrap selections in quotes on the command line.

**By segment ID (common in GROMACS/CHARMM setups):**
```bash
--x-select "segid A and protein"
--y-select "segid B and protein"
```

**By chain ID (common in PDB files):**
```bash
--x-select "chainID A and protein"
--y-select "chainID B and protein"
```

**By residue number range:**
```bash
--x-select "resid 1:150"
--y-select "resid 151:300"
```

**Excluding water and ions (useful for getting just protein):**
```bash
--x-select "segid A and not (resname WAT SOL HOH NA CL)"
```

**Combining conditions:**
```bash
--x-select "segid A and protein and resid 50:100"
```

> **Not sure what segment IDs your system uses?** Load your topology in Python and check:
> ```python
> import MDAnalysis as mda
> u = mda.Universe("your_system.tpr")
> print(set(u.atoms.segids))     # all segment IDs
> print(set(u.atoms.chainIDs))   # all chain IDs
> ```

---

## Output Files

Every run always produces these three files:

### `<prefix>.npz` — Main data file

A compressed NumPy archive. Load it with `np.load(...)`. Contains:

| Key | Shape | Description |
|---|---|---|
| `mean_distance` | `(n_x_res, n_y_res)` | Average minimum distance per residue pair, over all frames |
| `min_distance` | `(n_x_res, n_y_res)` | Closest observed distance per residue pair across all frames. `NaN` if a pair was never within any calculable distance. |
| `contact_frequency` | `(n_x_res, n_y_res)` | Fraction of frames (0.0–1.0) each pair was in contact |
| `first_frame_contacts` | `(n_x_res, n_y_res)` | Binary map (0/1) of contacts in the first analyzed frame |
| `x_labels` | `(n_x_res,)` | Residue labels for rows (e.g. `"A:GLU42"`) |
| `y_labels` | `(n_y_res,)` | Residue labels for columns |
| `units` | scalar | Distance unit used (`"angstrom"` or `"nm"`) |
| `cutoff` | scalar | Cutoff value as provided |
| `n_frames` | scalar | Number of frames analyzed |

### `<prefix>_metadata.json` — Run parameters

A plain-text JSON file recording every parameter used in the run. Useful for reproducing results or keeping a record of your analysis.

### `<prefix>_labels.csv` — Residue labels

A two-column CSV listing all residue labels for both groups. Handy for cross-referencing with other tools.

### Optional CSV matrices (with `--csv`)

When `--csv` is passed, four additional human-readable files are written — one for each matrix above. These can be opened directly in Excel or loaded into R. They are larger than the `.npz` but easier to inspect manually.

---

## Loading Results in Python

```python
import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.load("results/my_analysis.npz", allow_pickle=True)

# Pull out the arrays
contact_freq = data["contact_frequency"]
x_labels     = data["x_labels"]
y_labels     = data["y_labels"]

print(f"Matrix shape: {contact_freq.shape}")  # (n_x_res, n_y_res)
print(f"Units: {data['units']}")
print(f"Frames analyzed: {data['n_frames']}")

# Find residue pairs in contact > 50% of the time
rows, cols = np.where(contact_freq > 0.5)
for r, c in zip(rows, cols):
    print(f"{x_labels[r]} — {y_labels[c]}: {contact_freq[r, c]:.1%}")

# Plot a heatmap of contact frequency
plt.figure(figsize=(12, 10))
plt.imshow(contact_freq, aspect="auto", cmap="hot_r", vmin=0, vmax=1)
plt.colorbar(label="Contact frequency")
plt.xlabel("Group B residues")
plt.ylabel("Group A residues")
plt.title("Interface contact frequency")
plt.tight_layout()
plt.savefig("contact_map.png", dpi=150)
plt.show()
```

---

## Common Recipes

**Analyze only the last 500 ns of a 1 µs trajectory:**
```bash
python interface_mdmat.py \
    -s system.tpr \
    -f traj_1us.xtc \
    --x-select "segid A and protein" \
    --y-select "segid B and protein" \
    --begin 500000 \
    -o output/last500ns
```

**Speed things up by analyzing every 10th frame:**
```bash
python interface_mdmat.py \
    -s system.tpr \
    -f traj.xtc \
    --x-select "segid A and protein" \
    --y-select "segid B and protein" \
    --stride 10 \
    -o output/strided
```

**Use a tighter 3.5 Å cutoff and output distances in nm:**
```bash
python interface_mdmat.py \
    -s system.tpr \
    -f traj.xtc \
    --x-select "segid A and protein" \
    --y-select "segid B and protein" \
    --cutoff 3.5 \
    --cutoff-unit angstrom \
    --units nm \
    -o output/tight_cutoff
```

**Write CSV files for inspection in Excel:**
```bash
python interface_mdmat.py \
    -s system.tpr \
    -f traj.xtc \
    --x-select "segid A and protein" \
    --y-select "segid B and protein" \
    --csv \
    -o output/my_run
```

**Use a PDB + DCD trajectory (NAMD/CHARMM users):**
```bash
python interface_mdmat.py \
    -s system.psf \
    -f traj.dcd \
    --x-select "chainID A" \
    --y-select "chainID B" \
    -o output/namd_run
```

---

## Troubleshooting

**`ERROR: x selection returned 0 atoms`**
Your selection string didn't match any atoms. Check the selection syntax section above. Use the Python snippet there to print the segment IDs and chain IDs in your system to make sure you're using the right identifiers.

**`ERROR: No frames selected. Check --begin, --end, and --stride.`**
The time window you specified (via `--begin` / `--end`) doesn't overlap with your trajectory. Times must be in **picoseconds**. Check the length of your trajectory with:
```python
import MDAnalysis as mda
u = mda.Universe("system.tpr", "traj.xtc")
print(f"Trajectory runs from {u.trajectory[0].time:.0f} to {u.trajectory[-1].time:.0f} ps")
```

**Script is very slow**
- Use `--stride` to skip frames (e.g. `--stride 10`)
- Make sure `--backend OpenMP` is active (it is by default)
- Narrow the time window with `--begin` / `--end`
- Check that your selections aren't accidentally including solvent atoms

**`ImportError` or `ModuleNotFoundError`**
The required packages aren't installed in the Python environment you're running. Make sure you've activated the right conda environment (`conda activate mdanalysis`) before running the script.

**Output prefix contains a dot (e.g. `run_v2.0`)**
This is handled correctly — the script uses string concatenation for output paths, so `run_v2.0` will produce `run_v2.0.npz`, not `run_v2.npz`.
