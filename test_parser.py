#!/usr/bin/env python3

from pathlib import Path
import numpy as np

from parse_output import parse_frdex_output


files = sorted(
    Path(".").glob("*.out")
)

print(f"Found {len(files)} output files.\n")

for f in files:

    try:

        H, mu, vel, mag = parse_frdex_output(f)

        assert mu.shape[0] == H.shape[0]
        assert vel.shape[0] == H.shape[0]
        assert mag.shape[0] == H.shape[0]

        print(
            f"{f.name:45s}"
            f"H={H.shape}"
            f"diag={np.round(np.diag(H),4)}"
            f"symmetric={np.allclose(H, H.T)}"
            f" mu={mu.shape}"
            f" vel={vel.shape}"
            f" mag={mag.shape}"
        )

    except Exception as e:

        print(
            f"{f.name:45s}"
            f"FAILED : {e}"
        )