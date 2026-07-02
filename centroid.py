#!/usr/bin/env python3
"""
geometry.py

Utilities for reading fragment XYZ files and computing fragment centroids.
"""

from pathlib import Path
import numpy as np


def read_xyz(filename):
    """
    Read an XYZ file.

    Parameters
    ----------
    filename : str

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

        coords.append(
            [float(x) for x in toks[1:4]]
        )

    return atoms, np.array(coords)


def compute_centroid(coords):
    """
    Geometric centroid.
    """

    return coords.mean(axis=0)

def compute_fragment_centroids(
        xyz_pattern="frag*.xyz",
        outfile="fragment_centroids.dat"
):
    """
    Compute centroids for every fragment xyz.

    Returns
    -------
    dict

        {
            1 : np.array([x,y,z]),
            2 : np.array([x,y,z]),
            ...
        }
    """

    xyz_files = sorted(
        Path(".").glob(xyz_pattern)
    )

    if len(xyz_files) == 0:

        raise RuntimeError(
            "No fragment xyz files found."
        )

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

        print(
            frag,
            centroids[frag]
        )