#!/usr/bin/env python3

import re
from dataclasses import dataclass, field


@dataclass
class ExcitedState:

    index: int
    energy: float
    wavelength: float
    oscillator: float

    transitions: list = field(default_factory=list)

    dipole = None

state_re = re.compile(
    r"Excited State\s+(\d+).*?([\d.]+)\s+eV.*?([\d.]+)\s+nm.*?f=([\d.]+)"
)

trans_re = re.compile(
    r"(\d+)\s*->\s*(\d+)\s+([-\d\.]+)"
)

def parse_gaussian_log(logfile):

    with open(logfile) as f:
        lines = f.readlines()

    states = []

    current = None
 
    state_re = re.compile(
    r"Excited State\s+(\d+).*?([\d.]+)\s+eV.*?([\d.]+)\s+nm.*?f=([\d.]+)"
)
    for line in lines:

        m = state_re.search(line)


        if m:

            current = ExcitedState(

                index=int(m.group(1)),
                energy=float(m.group(2)),
                wavelength=float(m.group(3)),
                oscillator=float(m.group(4))

            )

            states.append(current)

            continue

        if current is None:
            continue
        
        trans_re = re.compile(
    r"(\d+)\s*->\s*(\d+)\s+([-\d\.]+)"
)


        m = trans_re.search(line)

        if m:

            occ = int(m.group(1))
            virt = int(m.group(2))
            coeff = float(m.group(3))

            current.transitions.append(

                (occ, virt, coeff)

            )

    return states

def dominant_transition(state):

    if len(state.transitions) == 0:
        return None

    return max(

        state.transitions,

        key=lambda x: abs(x[2])

    )

def classify_state(state):

    if state.oscillator < 1e-3:

        return "dark"

    dom = dominant_transition(state)

    if dom is None:

        return "unknown"

    occ, virt, coeff = dom

    if abs(coeff) > 0.6:

        return "LE"

    return "mixed"

def suggest_le_states(states, nstates=2):

    candidates = []

    for state in states:

        label = classify_state(state)

        if label != "LE":
            continue

        score = state.oscillator

        dom = dominant_transition(state)

        if dom:

            score *= abs(dom[2])

        candidates.append(

            (score, state)

        )

    candidates.sort(

        reverse=True,

        key=lambda x: x[0]

    )

    return [

        f"s{s.index}"

        for _, s in candidates[:nstates]

    ]

def write_report(states, outfile):

    with open(outfile, "w") as f:

        f.write(

            "State  Energy(eV)     f     Type\n"

        )

        for s in states:

            f.write(

                f"s{s.index:2d} "

                f"{s.energy:8.3f} "

                f"{s.oscillator:8.4f} "

                f"{classify_state(s)}\n"

            )