#!/usr/bin/env python3
"""
Residue-residue interface distance/contact matrix generator for MD trajectories.

Computes matrices between two MDAnalysis selections:
  - mean minimum residue-residue distance
  - minimum observed residue-residue distance
  - contact frequency over selected frames
  - binary contact map from the first analyzed frame

Designed as a between-group version of `gmx mdmat`.

Example:
    python interface_mdmat.py \
        -s cont2_1us_rep1.tpr \
        -f last500ns.xtc \
        --x-select "segid A and protein" \
        --y-select "segid B and protein" \
        --cutoff 4 \
        --cutoff-unit angstrom \
        --units angstrom \
        -o elfn2_mglur_last500
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array
except ImportError as exc:
    raise SystemExit(
        "ERROR: This script requires MDAnalysis.\n"
        "Install with: conda install -c conda-forge mdanalysis\n"
        "or: pip install MDAnalysis"
    ) from exc


ANGSTROM_TO_NM = 0.1
NM_TO_ANGSTROM = 10.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate residue-residue interface distance/contact matrices "
            "between two MDAnalysis atom selections."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-s", "--topology", required=True, help="Topology file, e.g. .tpr, .pdb, .gro, .psf")
    parser.add_argument("-f", "--trajectory", required=True, help="Trajectory file, e.g. .xtc, .trr, .dcd")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix")

    parser.add_argument(
        "--x-select",
        required=True,
        help='MDAnalysis selection for x-axis group, e.g. "segid A and protein"',
    )
    parser.add_argument(
        "--y-select",
        required=True,
        help='MDAnalysis selection for y-axis group, e.g. "segid B and protein"',
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=4.0,
        help="Contact cutoff distance.",
    )
    parser.add_argument(
        "--cutoff-unit",
        choices=["angstrom", "nm"],
        default="angstrom",
        help="Unit for --cutoff.",
    )
    parser.add_argument(
        "--units",
        choices=["angstrom", "nm"],
        default="angstrom",
        help="Distance unit for output matrices.",
    )

    parser.add_argument(
        "--begin",
        type=float,
        default=None,
        help="Start time in ps. If omitted, starts from first frame.",
    )
    parser.add_argument(
        "--end",
        type=float,
        default=None,
        help="End time in ps. If omitted, runs through final frame.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Analyze every Nth frame.",
    )

    parser.add_argument(
        "--backend",
        choices=["serial", "OpenMP"],
        default="OpenMP",
        help="MDAnalysis distance calculation backend.",
    )

    parser.add_argument(
        "--csv",
        action="store_true",
        help="Also write human-readable CSV matrices. NPZ output is always written.",
    )

    return parser.parse_args()


def cutoff_to_angstrom(cutoff: float, unit: str) -> float:
    if unit == "angstrom":
        return cutoff
    if unit == "nm":
        return cutoff * NM_TO_ANGSTROM
    raise ValueError(f"Unsupported cutoff unit: {unit}")


def scale_from_angstrom(unit: str) -> float:
    if unit == "angstrom":
        return 1.0
    if unit == "nm":
        return ANGSTROM_TO_NM
    raise ValueError(f"Unsupported output unit: {unit}")


def residue_labels(residues) -> list[str]:
    labels = []
    for res in residues:
        segid = getattr(res, "segid", "").strip()
        chain = getattr(res, "chainID", "").strip()
        resid = getattr(res, "resid", "")
        resname = getattr(res, "resname", "").strip()

        prefix = segid or chain
        if prefix:
            labels.append(f"{prefix}:{resname}{resid}")
        else:
            labels.append(f"{resname}{resid}")
    return labels


def selected_frame_indices(universe: mda.Universe, begin: float | None, end: float | None, stride: int) -> list[int]:
    indices = []

    for ts in universe.trajectory:
        time_ps = float(ts.time)

        if begin is not None and time_ps < begin:
            continue
        if end is not None and time_ps > end:
            continue
        if ts.frame % stride != 0:
            continue

        indices.append(ts.frame)

    return indices


def build_atom_residue_index(atomgroup, residues) -> np.ndarray:
    residue_to_idx = {res.ix: i for i, res in enumerate(residues)}
    return np.asarray([residue_to_idx[atom.residue.ix] for atom in atomgroup], dtype=np.int32)


def reduce_atom_distances_to_residue_min(
    atom_distances: np.ndarray,
    x_atom_residx: np.ndarray,
    y_atom_residx: np.ndarray,
    n_x_res: int,
    n_y_res: int,
) -> np.ndarray:
    """
    Convert atom-atom distances into residue-residue minimum distances.

    Uses np.minimum.at to avoid Python-level residue-pair loops.
    """
    flat_size = n_x_res * n_y_res
    residue_min_flat = np.full(flat_size, np.inf, dtype=np.float64)

    x_idx = np.repeat(x_atom_residx, len(y_atom_residx))
    y_idx = np.tile(y_atom_residx, len(x_atom_residx))
    flat_idx = x_idx * n_y_res + y_idx

    np.minimum.at(residue_min_flat, flat_idx, atom_distances.ravel())

    return residue_min_flat.reshape((n_x_res, n_y_res))


def write_csv_matrix(path: Path, matrix: np.ndarray, row_labels: list[str], col_labels: list[str]) -> None:
    df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
    df.index.name = "x_residue"
    df.to_csv(path)


def main() -> None:
    args = parse_args()

    if args.stride < 1:
        raise SystemExit("ERROR: --stride must be >= 1")

    out_prefix = Path(args.output_prefix)
    cutoff_angstrom = cutoff_to_angstrom(args.cutoff, args.cutoff_unit)
    output_scale = scale_from_angstrom(args.units)

    print("[INFO] Loading trajectory...")
    u = mda.Universe(args.topology, args.trajectory)

    x_group = u.select_atoms(args.x_select)
    y_group = u.select_atoms(args.y_select)

    if len(x_group) == 0:
        raise SystemExit(f"ERROR: x selection returned 0 atoms: {args.x_select}")
    if len(y_group) == 0:
        raise SystemExit(f"ERROR: y selection returned 0 atoms: {args.y_select}")

    x_residues = x_group.residues
    y_residues = y_group.residues

    x_labels = residue_labels(x_residues)
    y_labels = residue_labels(y_residues)

    x_atom_residx = build_atom_residue_index(x_group, x_residues)
    y_atom_residx = build_atom_residue_index(y_group, y_residues)

    n_x_res = len(x_residues)
    n_y_res = len(y_residues)

    print(f"[INFO] X group: {len(x_group)} atoms, {n_x_res} residues")
    print(f"[INFO] Y group: {len(y_group)} atoms, {n_y_res} residues")
    print(f"[INFO] Contact cutoff: {cutoff_angstrom:.3f} Å")

    frame_indices = selected_frame_indices(u, args.begin, args.end, args.stride)

    if not frame_indices:
        raise SystemExit("ERROR: No frames selected. Check --begin, --end, and --stride.")

    print(f"[INFO] Frames selected: {len(frame_indices)}")

    sum_distance = np.zeros((n_x_res, n_y_res), dtype=np.float64)
    min_distance = np.full((n_x_res, n_y_res), np.inf, dtype=np.float64)
    contact_counts = np.zeros((n_x_res, n_y_res), dtype=np.uint32)
    first_frame_contacts = None

    for i, frame_idx in enumerate(frame_indices, start=1):
        u.trajectory[frame_idx]

        atom_dist = distance_array(
            x_group.positions,
            y_group.positions,
            backend=args.backend,
        )

        residue_min = reduce_atom_distances_to_residue_min(
            atom_dist,
            x_atom_residx,
            y_atom_residx,
            n_x_res,
            n_y_res,
        )

        contacts = residue_min <= cutoff_angstrom

        sum_distance += residue_min
        min_distance = np.minimum(min_distance, residue_min)
        contact_counts += contacts.astype(np.uint32)

        if first_frame_contacts is None:
            first_frame_contacts = contacts.astype(np.uint8)

        if i % 25 == 0 or i == len(frame_indices):
            print(f"[INFO] Processed {i}/{len(frame_indices)} frames", flush=True)

    mean_distance = sum_distance / len(frame_indices)
    contact_frequency = contact_counts.astype(np.float64) / len(frame_indices)

    mean_distance_out = mean_distance * output_scale
    min_distance_out = min_distance * output_scale

    npz_path = out_prefix.with_suffix(".npz")

    np.savez_compressed(
        npz_path,
        mean_distance=mean_distance_out,
        min_distance=min_distance_out,
        contact_frequency=contact_frequency,
        first_frame_contacts=first_frame_contacts,
        x_labels=np.asarray(x_labels, dtype=object),
        y_labels=np.asarray(y_labels, dtype=object),
        units=args.units,
        cutoff=args.cutoff,
        cutoff_unit=args.cutoff_unit,
        cutoff_angstrom=cutoff_angstrom,
        n_frames=len(frame_indices),
    )

    metadata = {
        "topology": args.topology,
        "trajectory": args.trajectory,
        "x_selection": args.x_select,
        "y_selection": args.y_select,
        "x_n_atoms": int(len(x_group)),
        "y_n_atoms": int(len(y_group)),
        "x_n_residues": int(n_x_res),
        "y_n_residues": int(n_y_res),
        "cutoff": args.cutoff,
        "cutoff_unit": args.cutoff_unit,
        "cutoff_angstrom": cutoff_angstrom,
        "output_units": args.units,
        "begin_ps": args.begin,
        "end_ps": args.end,
        "stride": args.stride,
        "n_frames_analyzed": len(frame_indices),
        "distance_backend": args.backend,
        "outputs": {
            "npz": str(npz_path),
            "mean_distance": "mean_distance",
            "min_distance": "min_distance",
            "contact_frequency": "contact_frequency",
            "first_frame_contacts": "first_frame_contacts",
        },
    }

    json_path = out_prefix.with_name(out_prefix.name + "_metadata.json")
    json_path.write_text(json.dumps(metadata, indent=2))

    labels_path = out_prefix.with_name(out_prefix.name + "_labels.csv")
    labels_df = pd.concat(
        [
            pd.DataFrame({"axis": "x", "index": np.arange(n_x_res), "label": x_labels}),
            pd.DataFrame({"axis": "y", "index": np.arange(n_y_res), "label": y_labels}),
        ],
        ignore_index=True,
    )
    labels_df.to_csv(labels_path, index=False)

    if args.csv:
        write_csv_matrix(
            out_prefix.with_name(out_prefix.name + "_mean_distance.csv"),
            mean_distance_out,
            x_labels,
            y_labels,
        )
        write_csv_matrix(
            out_prefix.with_name(out_prefix.name + "_min_distance.csv"),
            min_distance_out,
            x_labels,
            y_labels,
        )
        write_csv_matrix(
            out_prefix.with_name(out_prefix.name + "_contact_frequency.csv"),
            contact_frequency,
            x_labels,
            y_labels,
        )
        write_csv_matrix(
            out_prefix.with_name(out_prefix.name + "_first_frame_contacts.csv"),
            first_frame_contacts,
            x_labels,
            y_labels,
        )

    print("[DONE] Wrote:")
    print(f"  {npz_path}")
    print(f"  {json_path}")
    print(f"  {labels_path}")

    if args.csv:
        print("  CSV matrices")


if __name__ == "__main__":
    main()
