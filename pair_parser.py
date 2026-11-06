#!/usr/bin/env python3

import argparse

from parse_input import parse_frdex_input
from parse_output import parse_frdex_output
from models import FrDExPairResult
from pathlib import Path

import re


def parse_pair(inpfile, outfile):
    
    inpfile = Path(inpfile)
    outfile = Path(outfile)

    states = parse_frdex_input(inpfile)

    H, mu, vel, mag = parse_frdex_output(outfile)

    m = re.search(
    r"frag(\d+)_frag(\d+)",
    inpfile.name
    )

    print(f"group 1: {m.group(1)}")
    print(f"group 2: {m.group(2)}")

    pair = (
        int(m.group(1)),
        int(m.group(2))
    )

    return FrDExPairResult(

        pair=pair,

        states=states,

        hamiltonian=H,

        dipoles=mu,

        velocity_dipoles=vel,

        magnetic_dipoles=mag

    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Parser Test input/output files from FrDEx pair calculations."
    )

    parser.add_argument(
        "-i",
        "--input",
        default=".",
        help="FrDEx *.inp files (default: current directory)"
    )

    parser.add_argument(
        "-o",
        "--output",
        default=".",
        help="FrDEx *.out files (default: current directory)"
    )

    args = parser.parse_args()


    result = parse_pair(args.input, args.output)

    print(result.pair)
    
    print(result.states)
    
    print(f' Hamiltonian shape: {result.hamiltonian.shape}')

    print(result.hamiltonian)