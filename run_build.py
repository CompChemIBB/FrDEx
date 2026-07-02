# run_build.py

from parse_output import *
from build_exat import *

H,mu,vel,mag = parse_frdex_output(
    "diab_fratd_resp_4ct_imot_frag1_frag2.out"
)

write_site(H)

write_coup(H)

write_vectors(mu,"dipo.in")

write_vectors(
    vel,
    "dipovel.in"
)

write_vectors(
    mag,
    "dipomag.in"
)