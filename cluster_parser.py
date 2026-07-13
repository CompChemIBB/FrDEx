#!/usr/bin/env python3

from pathlib import Path
from parse_input import parse_frdex_input
from parse_output import parse_frdex_output
from models import FrDExResult


def parse_cluster(inpfile, outfile):

    inpfile = Path(inpfile)
    outfile = Path(outfile)

    states = parse_frdex_input(inpfile)

    H, mu, vel, mag = parse_frdex_output(outfile)

    #
    # Recover fragments from the LE states
    #

    fragments = sorted({

        s.fragment

        for s in states

        if s.state_type == "LE"

    })

    return FrDExResult(

        calculation_type="cluster",

        fragments=fragments,

        states=states,

        hamiltonian=H,

        dipoles=mu,

        velocity_dipoles=vel,

        magnetic_dipoles=mag

    )