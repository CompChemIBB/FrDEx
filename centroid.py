#!/usr/bin/env python3
"""
centroid.py

Utilities for reading fragment XYZ files and computing fragment centroids.
"""

from pathlib import Path
import numpy as np


def read_xyz(filename):
    """
    Read an XYZ file.

    Parameters
    ----------
    filename : str or Path

    Returns
    -------
    atoms : list[str]
    coords : ndarray (N,3)
    """

    with open(filename) as f:
        lines = f.readlines()

    natoms = int(lines[0])

    atoms = []
    coords = []

    for line in lines[2:2 + natoms]:
        toks = line.split()

        atoms.append(toks[0])
        coords.append([float(x) for x in toks[1:4]])

    return atoms, np.array(coords)


def compute_centroid(coords):
    """
    Compute the geometric centroid.
    """
    return coords.mean(axis=0)


def compute_fragment_centroids(
    directory=".",
    xyz_pattern="frag*.xyz",
    outfile="fragment_centroids.dat",
):
    """
    Compute centroids for every fragment XYZ file.

    Parameters
    ----------
    directory : str or Path
        Directory containing frag*.xyz files.

    xyz_pattern : str
        Pattern used to locate fragment xyz files.

    outfile : str
        Output filename written inside `directory`.

    Returns
    -------
    dict

        {
            1 : np.array([x,y,z]),
            2 : np.array([x,y,z]),
            ...
        }
    """

    directory = Path(directory)

    xyz_files = sorted(directory.glob(xyz_pattern))

    if not xyz_files:
        raise RuntimeError(
            f"No fragment XYZ files matching '{xyz_pattern}' "
            f"found in '{directory}'."
        )

    outfile = directory / outfile

    centroids = {}

    with open(outfile, "w") as out:

        out.write(
            "#Fragment        X              Y              Z\n"
        )

        for xyz in xyz_files:

            atoms, coords = read_xyz(xyz)

            c = compute_centroid(coords)

            frag = int(
                xyz.stem.replace("frag", "")
            )

            centroids[frag] = c

            out.write(
                f"{frag:4d}"
                f"{c[0]:15.6f}"
                f"{c[1]:15.6f}"
                f"{c[2]:15.6f}\n"
            )

    return centroids


if __name__ == "__main__":

    centroids = compute_fragment_centroids()

    print()

    for frag in sorted(centroids):
        print(frag, centroids[frag])