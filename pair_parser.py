#!/usr/bin/env python3

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

    result = parse_pair(

    "diab_fratd_resp_2ct_imot_frag1_frag2.inp",

    "diab_fratd_resp_2ct_imot_frag1_frag2.out"

    )
    
    print(result.pair)
    
    print(result.states)
    
    print(f' Hamiltonian shape: {result.hamiltonian.shape}')

    print(result.hamiltonian)