#!/usr/bin/env python3

"""
Generate OVERDIA input files for all unique fragment pairs
and create a run_diab.sh script to execute them.

Equivalent to the original Bash script.
"""

from pathlib import Path
import argparse

# ------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Generate OVERDIA input files for all unique fragment pairs."
)
parser.add_argument("NFRAGS", type=int, help="number of fragments")
parser.add_argument('-s', '--system-name', default='G4-1qdf', help='System name for output files (default: G4-1qdf)')
args = parser.parse_args()
NFRAGS = args.NFRAGS

OVERDIA_EXE = (
    "/home/roberto/JAMES/Programs/OVERDIA-REF/FRAG-ORB/overdia-par.e"
)

SYSTEM_NAME = args.system_name

# ------------------------------------------------------------------
# Create run script
# ------------------------------------------------------------------

run_script = Path("run_diab.sh")

with run_script.open("w") as runfile:

    runfile.write("#!/bin/bash\n\n")

    for frag1 in range(1, NFRAGS + 1):

        for frag2 in range(frag1 + 1, NFRAGS + 1):

            inp_name = (
                f"diab_{SYSTEM_NAME}_frag{frag1}_frag{frag2}.inp"
            )

            out_name = (
                f"diab_{SYSTEM_NAME}_frag{frag1}_frag{frag2}.out"
            )

            pair_log = (
                f"{SYSTEM_NAME}_frag{frag1}_frag{frag2}-mon_sub-diab.log"
            )

            frag1_fchk = (
                f"{SYSTEM_NAME}_frag{frag1}-mon_sub-diab.fchk"
            )

            frag2_fchk = (
                f"{SYSTEM_NAME}_frag{frag2}-mon_sub-diab.fchk"
            )

            lines = [
                "'CONCISE'                        printing option",
                f"'{pair_log}'  log file of pair",
                "40                                 number of TD states",
                "'FRAMO'     possible ops KSORB/MXORB/FRAMO/FRATD/",
                "8                 number of diabatic states from input",
                f"'{frag1_fchk}' '{frag2_fchk}'     fchk files",
                "2   4             number of occ and virtual orb from input",
                "0 'FR1' 0 'FR2'             occ orbitals expressed on these KS occ orbitals 0 is HOMO",
                "1 'FR1' 2 'FR1' 1 'FR2' 2 'FR2'      virt orbitals expressed on these KS virtual orb   1 is LUMO",
                "1.0 0.0     coef  Occ    mixing the orbitals",
                "0.0 1.0     coef  Occ    mixing the orbitals",
                "1.0 0.0 0.0 0.0    coef  Virt   mixing the orbitals",
                "0.0 1.0 0.0 0.0    coef  Virt   mixing the orbitals",
                "0.0 0.0 1.0 0.0    coef  Virt   mixing the orbitals",
                "0.0 0.0 0.0 1.0    coef  Virt   mixing the orbitals",
                "1                no components first state",
                "1                no components second state",
                "1                no components third state",
                "1                no components fourth state",
                "1                no components fifth state",
                "1                no components sixth state",
                "1                no components seventh state",
                "1                no components eighth state",
                "1    1   1.d0",
                "1    2   1.d0",
                "1    3   1.d0",
                "1    4   1.d0",
                "2    1   1.d0",
                "2    2   1.d0",
                "2    3   1.d0",
                "2    4   1.d0",
            ]

            with open(inp_name, "w") as inp:
                inp.write("\n".join(lines))
                inp.write("\n")

            runfile.write(
                f"{OVERDIA_EXE} < {inp_name} > {out_name}\n"
            )

print(
    f"Generated {NFRAGS*(NFRAGS-1)//2} input files "
    f"and run_diab.sh"
)