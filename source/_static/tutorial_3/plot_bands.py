#!/usr/bin/env python3
"""plot_bands.py: Plot ONETEP band structure from an ONETEP input file and a
CASTEP-format eigenvalue file.

The band structure is plotted from the eigenvalues written by the ONETEP
properties module (set with :code:`do_properties : T`), which are stored in a
file with the same format as the :code:`.bands` file of a standard CASTEP
calculation.

Parsing is done with ASE:

  * input file        -> ase.io.onetep.read_onetep_in (lattice + kpoint list)
  * eigenvalue file   -> ase.io.castep.read_bands     (kpts, eigenvalues, E_F)

The inline high-symmetry labels (with `!`) are recovered from the input file
with a regex, because ASE strips them on read.

Usage
-----
    python plot_bands.py INPUT.dat VAL_BANDS [-n N] [-o OUT.PNG]
                         [--emin E] [--emax E]

Options
-------
input                  ONETEP input file (.dat).
bands                  CASTEP-format eigenvalue file (.val_bands).
-n, --skip <n>         Ignore the first N k-points; by default the number of
                       k-points with non-zero weight is auto-detected.
-o, --output <file>    Output image file. Default: bands.png.
--emin <value>         Minimum of the energy axis (eV).
--emax <value>         Maximum of the energy axis (eV).

Output
------
A PNG image of the band structure, with the energy relative to the Fermi
level (E - E_F) on the vertical axis and the high-symmetry points labelled
along the horizontal axis.
"""

import argparse
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from ase.io.onetep import read_onetep_in
from ase.io.castep import read_bands

matplotlib.use("Agg")

LABEL_MAP = {"GAMMA": "Γ"}


def read_input(path):
    """Parse an ONETEP input file.

    Returns
    -------
    cell : (3, 3) ndarray
        Lattice vectors in Angstrom (rows).
    kpoints : (N, 4) ndarray
        Columns [kx, ky, kz, weight] in fractional coordinates.
    labels : list of (int, str)
        Kpoint index (0-based) and high-symmetry label for each inline
        '!LABEL' comment in the kpoints_list block, in file order.
    """
    with open(path) as fd:
        parsed = read_onetep_in(fd)
    cell = np.asarray(parsed["atoms"].cell[:], dtype=float)
    klines = parsed["keywords"]["kpoints_list"]
    kpoints = np.loadtxt(klines)

    labels = []
    idx = 0
    in_block = False
    with open(path) as fd:
        for line in fd:
            low = line.lower()
            if not in_block and "%block" in low and "kpoints_list" in low:
                in_block = True
                continue
            if in_block and "%endblock" in low:
                break
            if in_block:
                m = re.search(r"!\s*(\S+)", line)
                if m:
                    labels.append((idx, LABEL_MAP.get(m.group(1).upper(), m.group(1).upper())))
                idx += 1
    return cell, kpoints, labels


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="ONETEP input file (.dat)")
    ap.add_argument("bands", help="CASTEP-format eigenvalue file (.val_bands)")
    ap.add_argument("-n", "--skip", type=int, default=None,
                    help="ignore the first N k-points; default: auto-detect "
                         "(number of k-points with nonzero weight)")
    ap.add_argument("-o", "--output", default="bands.png", help="output PNG")
    ap.add_argument("--emin", type=float, default=None, help="y-axis minimum (eV)")
    ap.add_argument("--emax", type=float, default=None, help="y-axis maximum (eV)")
    args = ap.parse_args()

    try:
        cell, kpoints, labels = read_input(args.input)
    except Exception as exc:
        sys.exit(f"Failed to read input file: {exc}")
    try:
        kpts, weights, eigenvalues, efermi = read_bands(args.bands)
    except Exception as exc:
        sys.exit(f"Failed to read eigenvalue file: {exc}")

    n_kpts = kpoints.shape[0]
    if kpts.shape[0] != n_kpts:
        sys.exit(f"k-point count mismatch: input has {n_kpts}, bands file has {kpts.shape[0]}")
    if not np.allclose(kpoints[:, :3], kpts, atol=1e-4):
        sys.exit("k-point coordinates in input file and bands file do not match")

    if args.skip is None:
        skip = int((kpoints[:, 3] > 0).sum())
        print(f"Auto-detected skip = {skip} (k-points with nonzero weight)")
    else:
        skip = args.skip
        if skip < 0:
            sys.exit("--skip must be >= 0")
    if skip >= n_kpts:
        sys.exit(f"--skip ({skip}) leaves no k-points (total {n_kpts})")

    eig = eigenvalues[:, skip:, :]
    frac = kpoints[skip:, :3]
    n_spin, n_band = eig.shape[0], eig.shape[2]

    recip = 2.0 * np.pi * np.linalg.inv(cell).T
    k_cart = frac @ recip
    dist = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(k_cart, axis=0), axis=1))))

    tick_pos, tick_lab = [], []
    for idx, label in labels:
        if not (skip <= idx < n_kpts):
            continue
        x = dist[idx - skip]
        if not tick_pos or abs(x - tick_pos[-1]) > 1e-8:
            tick_pos.append(x)
            tick_lab.append(label)

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    for spin in range(n_spin):
        for b in range(n_band):
            ax.plot(dist, eig[spin, :, b] - efermi, lw=1.0, color=f"C{spin}")
    ax.axhline(0.0, color="k", ls="--", lw=0.8)
    for x in tick_pos:
        ax.axvline(x, color="k", ls=":", lw=0.6)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab)
    ax.set_xlim(dist[0], dist[-1])
    ax.set_ylabel("$E - E_F$ (eV)")
    if args.emin is not None or args.emax is not None:
        ax.set_ylim(args.emin, args.emax)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
