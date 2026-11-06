#!/usr/bin/env python3

from pathlib import Path

from parse_input import parse_frdex_input
from parse_output import parse_frdex_output
from models import FrDExPairResult


def parse_full_system(inpfile, outfile):
    """
    Parse a full OVERDIA calculation (all fragments together).

    Unlike pair calculations there is no fragment pair encoded
    in the filename.
    """

    inpfile = Path(inpfile)
    outfile = Path(outfile)

    #
    # Parse diabatic-state definitions
    #
    states = parse_frdex_input(inpfile)

    #
    # Parse OVERDIA output
    #
    H, mu, vel, mag = parse_frdex_output(outfile)

    #
    # Number of fragments
    #
    nfrag = None

    with open(inpfile) as f:

        for line in f:

            toks = line.split()

            if len(toks) >= 2:

                try:
                    value = int(toks[0])
                    print(value, toks[1])

                    if "number of fragments" in line.lower():

                        nfrag = value
                        break

                except ValueError:
                    pass

    if nfrag is None:
        raise RuntimeError(
            "Cannot determine number of fragments."
        )

    #
    # Return the same object as the pair parser
    #
    return FrDExPairResult(

        pair=None,

        states=states,

        hamiltonian=H,

        dipoles=mu,

        velocity_dipoles=vel,

        magnetic_dipoles=mag,

    )


if __name__ == "__main__":

    result = parse_full_system(

        "diab_fratd_4ct_meGCGC.inp",

        "diab_fratd_4ct_meGCGC.out"

    )

    print(result.hamiltonian.shape)
    print(result.hamiltonian)