#!/usr/bin/env python3

import re
from dataclasses import dataclass, field


# ============================================================
# Data structure
# ============================================================

@dataclass
class ExcitedState:

    index: int
    energy: float          # eV
    wavelength: float      # nm
    oscillator: float

    transitions: list = field(default_factory=list)

    dipole = None


# ============================================================
# Regular expressions
# ============================================================

state_re = re.compile(
    r"Excited State\s+(\d+).*?([\d.]+)\s+eV.*?([\d.]+)\s+nm.*?f=([\d.]+)"
)

trans_re = re.compile(
    r"(\d+)\s*->\s*(\d+)\s+([-\d\.]+)"
)


# ============================================================
# Gaussian parser
# ============================================================

def parse_gaussian_log(logfile):

    with open(logfile) as f:
        lines = f.readlines()

    states = []
    current = None

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

        m = trans_re.search(line)

        if m:

            occ = int(m.group(1))
            virt = int(m.group(2))
            coeff = float(m.group(3))

            current.transitions.append(
                (occ, virt, coeff)
            )

    return states


# ============================================================
# Utilities
# ============================================================

def dominant_transition(state):

    if not state.transitions:
        return None

    return max(
        state.transitions,
        key=lambda x: abs(x[2])
    )


# ============================================================
# Simple classification
# ============================================================

def classify_state(state):

    if state.oscillator < 1.0e-3:
        return "dark"

    dom = dominant_transition(state)

    if dom is None:
        return "unknown"

    occ, virt, coeff = dom

    if abs(coeff) >= 0.50:
        return "LE"

    if abs(coeff) >= 0.30:
        return "mixed"

    return "multiconfigurational"


# ============================================================
# Score states for LE selection
# ============================================================

def score_state(state):

    label = classify_state(state)

    if label != "LE":
        return -1.0

    dom = dominant_transition(state)

    coeff = abs(dom[2]) if dom else 0.0

    # Prefer:
    #
    # 1) bright states
    # 2) single-configurational states
    # 3) lower-energy states
    #

    score = (
        10.0 * state.oscillator +
        coeff -
        0.02 * state.energy
    )

    return score


# ============================================================
# Suggest LE states
# ============================================================

# def suggest_le_states(states, nstates=2):

#     ranked = sorted(
#         states,
#         key=score_state,
#         reverse=True
#     )

#     selected = []

#     for state in ranked:

#         if classify_state(state) != "LE":
#             continue

#         selected.append(f"s{state.index}")

#         if len(selected) == nstates:
#             break

#     return selected


def suggest_le_states(
        states,
        nstates=2,
        emin=4.5,
        emax=6.90,
    ):
    
    candidates = []

    for state in states:
    
        if not (emin <= state.energy <= emax):
            continue
    
        if classify_state(state) == "dark":
            continue
    
        dom = dominant_transition(state)
    
        coeff = abs(dom[2]) if dom else 0.0
    
        score = (
            state.oscillator
            * coeff
        )
    
        candidates.append((score, state))
    
    candidates.sort(
    reverse=True,
    key=lambda x: x[0]
    )
    
    return [
        f"s{s.index}"
        for _, s in candidates[:nstates]
    ]

# ============================================================
# Report
# ============================================================

def write_report(states, outfile):

    with open(outfile, "w") as f:

        f.write(
            "{:>5s} {:>8s} {:>8s} {:>8s} {:>10s} {:>12s}\n".format(
                "State",
                "E(eV)",
                "f",
                "Coeff",
                "Type",
                "Score"
            )
        )

        f.write("-"*60 + "\n")

        for state in states:

            dom = dominant_transition(state)

            coeff = abs(dom[2]) if dom else 0.0

            f.write(
                "{:>5s} {:8.3f} {:8.4f} {:8.3f} {:>10s} {:12.3f}\n".format(
                    f"s{state.index}",
                    state.energy,
                    state.oscillator,
                    coeff,
                    classify_state(state),
                    score_state(state)
                )
            )

        f.write("\n")

        best = suggest_le_states(states)

        f.write("Suggested LE states:\n")

        for s in best:
            f.write(f"   {s}\n")