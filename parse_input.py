#!/usr/bin/env python3

## ToDo parse CT states from start_define_states block
import re
from pathlib import Path
from typing import List
from models import State



# def parse_frdex_input(inpfile):
def parse_frdex_input(inpfile):

    with open(inpfile) as f:

        lines = f.readlines()

    states = []

    states.extend(
        parse_le_states(lines)
    )

    states.extend(
        parse_ct_states(lines)
    )

    return states

# EXAMPLE

def parse_le_states(lines):

    states = []

    inside = False

    for line in lines:

        line = line.strip()

        if line.startswith("start_read_states"):
            inside = True
            continue

        if line.startswith("end_read_states"):
            break

        if not inside:
            continue

        m = re.match(r"d(\d+)\s+(\d+)\s+s(\d+)", line)

        if m:

            diab, frag, root = m.groups()

            states.append(

                State(

                    label=f"LE_F{frag}_S{root}",

                    state_type="LE",

                    fragment=int(frag),

                    root=int(root),

                    local_index=int(diab)

                )

            )

    return states

def parse_ct_states(lines):

    states = []

    inside = False

    i = 0

    while i < len(lines):

        line = lines[i].strip()

        if line.startswith("start_define_states"):

            inside = True
            i += 1
            continue

        if line.startswith("end_define_states"):
            break

        if not inside:
            i += 1
            continue

        if "_start" in line:

            diab = int(
                re.search(r"d(\d+)_start", line).group(1)
            )

            transition = lines[i+2].strip()

            m = re.match(
                r"(\d+)H\d+\s+(\d+)L\d+",
                transition
            )

            if m:

                donor, acceptor = map(int, m.groups())

                states.append(

                    State(

                        label=f"CT_F{donor}_F{acceptor}",

                        state_type="CT",

                        donor=donor,

                        acceptor=acceptor,

                        local_index=diab

                    )

                )

        i += 1

    return states


if __name__ == "__main__":

    files = sorted(
    Path(".").glob("*.inp")
)


    print(f"Found {len(files)} output files.\n")

    for f in files:

        try:
            states = parse_frdex_input(f)

            for s in states:
                print(s)

        except Exception as exc:
            print(f"Failed to parse {f}: {exc}")
