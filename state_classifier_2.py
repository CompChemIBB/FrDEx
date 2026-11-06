# ============================================================
# Score states for LE selection
# ============================================================

def score_state(
    state,
    emin=4.5,
    emax=6.5,
):
    """
    Score a state for inclusion in the diabatic basis.

    Preference is given to

        1) LE states
        2) states inside the desired energy window
        3) bright transitions
        4) dominant single excitation
        5) lower excitation energies
    """

    # reject non-LE states
    if classify_state(state) != "LE":
        return -1e9

    # reject outside energy window
    if not (emin <= state.energy <= emax):
        return -1e9

    dom = dominant_transition(state)

    coeff = abs(dom[2]) if dom else 0.0

    score = (
        8.0 * state.oscillator
        + 2.0 * coeff
        - 0.30 * (state.energy - emin)
    )

    return score


# ============================================================
# Suggest LE states
# ============================================================

def suggest_le_states(
    states,
    nstates=3,
    emin=4.5,
    emax=6.5,
):
    """
    Return the n best LE states.
    """

    ranked = sorted(
        states,
        key=lambda s: score_state(s, emin, emax),
        reverse=True,
    )

    selected = []

    for state in ranked:

        if score_state(state, emin, emax) < 0:
            continue

        selected.append(f"s{state.index}")

        if len(selected) == nstates:
            break

    return selected


# ============================================================
# Report
# ============================================================

def write_report(
    states,
    outfile,
    nstates=3,
    emin=4.5,
    emax=6.5,
):

    with open(outfile, "w") as f:

        f.write(
            "{:>5s} {:>8s} {:>8s} {:>8s} {:>10s} {:>12s}\n".format(
                "State",
                "E(eV)",
                "f",
                "Coeff",
                "Type",
                "Score",
            )
        )

        f.write("-" * 65 + "\n")

        ranked = sorted(
            states,
            key=lambda s: score_state(s, emin, emax),
            reverse=True,
        )

        for state in ranked:

            dom = dominant_transition(state)

            coeff = abs(dom[2]) if dom else 0.0

            f.write(
                "{:>5s} {:8.3f} {:8.4f} {:8.3f} {:>10s} {:12.3f}\n".format(
                    f"s{state.index}",
                    state.energy,
                    state.oscillator,
                    coeff,
                    classify_state(state),
                    score_state(state, emin, emax),
                )
            )

        f.write("\n")

        best = suggest_le_states(
            states,
            nstates,
            emin,
            emax,
        )

        f.write("Suggested LE states:\n")

        for s in best:
            f.write(f"   {s}\n")