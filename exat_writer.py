#!/usr/bin/env python3

import numpy as np


def write_site_file(system, filename="site.in"):

    with open(filename, "w") as f:

        for i, energy in enumerate(system.site_energies, start=1):

            f.write(f"{i:4d} {energy:14.8f}\n")


def write_dipo_file(system, filename="dipo.in"):

    np.savetxt(
        filename,
        system.mu,
        fmt="%16.7f"
    )


def write_dipovel_file(system, filename="dipovel.in"):

    np.savetxt(
        filename,
        system.vel,
        fmt="%16.7f"
    )


def write_dipomag_file(system, filename="dipomag.in"):

    np.savetxt(
        filename,
        system.mag,
        fmt="%16.7f"
    )


def write_couplings(system, filename="coup.in"):
    """
    Write lower triangular Hamiltonian
    (without diagonal), exactly as EXAT expects.
    """

    H = system.H

    with open(filename, "w") as f:

        n = H.shape[0]

        for i in range(1, n):

            for j in range(i):

                #
                # Convert eV -> cm^-1
                #

                coupling = H[i, j] * 8065.54429

                f.write(f"{coupling:16.8f}\n")


def write_cent_file(system, filename="cent.in"):

    with open(filename, "w") as f:

        for center in system.state_centers:

            f.write(
                f"{center[0]:16.10f} "
                f"{center[1]:16.10f} "
                f"{center[2]:16.10f}\n"
            )


def write_chromlist(system, filename="chromlist.in"):

    with open(filename, "w") as f:

        for i in range(len(system.registry)):

            f.write(f"{i+1:d} 1\n")


def write_all(system, output_dir):

    write_site_file(system, filename=f"{output_dir}/site.in")

    write_couplings(system, filename=f"{output_dir}/coup.in")

    write_dipo_file(system, filename=f"{output_dir}/dipo.in")

    write_dipovel_file(system, filename=f"{output_dir}/dipovel.in")

    write_dipomag_file(system, filename=f"{output_dir}/dipomag.in")

    write_cent_file(system, filename=f"{output_dir}/cent.in")

    write_chromlist(system, filename=f"{output_dir}/chromlist.in")

    print()

    print("EXAT input files written successfully.")