#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np

def get_h_coords(anchor, nitro, cme, phi, theta, mod_c):
    """Calculates Methyl H coordinates based on anchor-Nitro-CMe geometry."""
    a = nitro - anchor
    b = cme - nitro
    mod_b = np.linalg.norm(b)
    d_vec = np.cross(a, b)
    mod_d = np.linalg.norm(d_vec)
    
    v_unit_b = b / mod_b
    v_unit_d = d_vec / mod_d
    v_unit_perp = np.cross(v_unit_d, v_unit_b)
    
    # Vector c from CMe to H
    c_vec = mod_c * (np.sin(theta)*np.cos(phi)*v_unit_perp + 
                     np.sin(theta)*np.sin(phi)*v_unit_d - 
                     np.cos(theta)*v_unit_b)
    return c_vec + cme

def parse_pdb(pdb_file):
    """Extracts base atoms and C1' for G, DC, and DCP residues."""
    residues = []
    current_res_id = None
    current_atoms = {}
    current_res_name = None

    # Atoms to keep for Guanine (G), Cytosine (DC/DCP), and Potassium (K)
    keep_list = ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "O6", "N7", "C8", "N9", "N2", 
                 "H1", "H8", "H21", "H22", "H3", "H5", "H6", "H41", "H42", "C1'", "K", "K+"]

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                atom_name = line[12:16].strip()
                
                if res_name in ["G", "K", "DC", "DCP", "DC5", "DC3"]:
                    if res_seq != current_res_id:
                        if current_atoms:
                            residues.append({'name': current_res_name, 'seq': current_res_id, 'atoms': current_atoms})
                        current_atoms = {}
                        current_res_id = res_seq
                        current_res_name = res_name
                    
                    if atom_name in keep_list:
                        coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        current_atoms[atom_name] = coords
        if current_atoms:
            residues.append({'name': current_res_name, 'seq': current_res_id, 'atoms': current_atoms})
    return residues

def write_fragment_xyz(res, f_id, phiC, thetaC, blC, base_name):
    """Generates an individual standard XYZ coordinate file for the given fragment."""
    xyz_lines = []
    atoms = res['atoms']
    
    if res['name'] == "K":
        for atom_name, xyz in atoms.items():
            xyz_lines.append(f"K   {xyz[0]:12.8f} {xyz[1]:12.8f} {xyz[2]:12.8f}")
    else:
        if res['name'] == "G":
            nitro, anchor = atoms["N9"], atoms["C8"]
        elif res['name'] in ["DC", "DCP", "DC5", "DC3"]:
            nitro, anchor = atoms["N1"], atoms["C2"]
            
        cme = atoms["C1'"]
        
        # 1. Base atoms (excluding C1')
        for name, xyz in atoms.items():
            if name == "C1'": continue
            xyz_lines.append(f"{name[:1]:3} {xyz[0]:12.8f} {xyz[1]:12.8f} {xyz[2]:12.8f}")
            
        # 2. Methyl Carbon
        xyz_lines.append(f"C   {cme[0]:12.8f} {cme[1]:12.8f} {cme[2]:12.8f}")
        
        # 3. Methyl Hydrogens
        for j in range(3):
            h_pos = get_h_coords(anchor, nitro, cme, phiC[j], thetaC[j], blC[j])
            xyz_lines.append(f"H   {h_pos[0]:12.8f} {h_pos[1]:12.8f} {h_pos[2]:12.8f}")
            
    # Write the formatted XYZ file
    # xyz_filename = f"{base_name}_frag_{f_id}_{res['name']}_{res['seq']}.xyz"
    xyz_filename = f"frag{f_id}.xyz"
    with open(xyz_filename, 'w') as xyz_f:
        xyz_f.write(f"{len(xyz_lines)}\n")
        xyz_f.write(f"Fragment {f_id} - Residue {res['name']}{res['seq']} extracted from PDB\n")
        for line in xyz_lines:
            xyz_f.write(line + "\n")

def main():
    # Use argparse to handle command-line flags cleanly
    parser = argparse.ArgumentParser(description="Convert PDB structures to Gaussian .com fragment files with optional XYZ file output.")
    parser.add_argument("pdb_file", help="Path to the input PDB file (e.g., g_quartet.pdb)")
    parser.add_argument("-x", "--xyz", action="store_true", help="Produce individual .xyz files for each fragment block.")
    
    args = parser.parse_args()

    pdb_name = args.pdb_file
    base_name = os.path.splitext(pdb_name)[0]
    residues = parse_pdb(pdb_name)
    
    if not residues:
        print(f"Error: No matching residues (G, DC, DCP, K) found in {pdb_name}.")
        sys.exit(1)

    # Methyl parameters
    phiC = [2.0944, -2.0944, 0.0]
    thetaC = [1.91602245, 1.91602245, 1.89577663]
    blC = [1.09, 1.09, 1.09]

    # Calculate Total System Charge (All organic methyl bases = 0 1, K+ ions = 1 1)
    fragment_charges = []
    for res in residues:
        if res['name'] == "K":
            charge, mult = 1, 1
        else:
            charge, mult = 0, 1
        fragment_charges.append((charge, mult))

    sys_charge = sum(c for c, m in fragment_charges)
    sys_mult = 1 if sys_charge % 2 == 0 else 2

    # --- Part 1: Writing the Gaussian Input File ---
    with open(f"{base_name}.com", "w") as f:
        f.write(f"%mem=96GB\n%nprocshared=48\n%chk={base_name}.chk\n")
        # Updated to standard Ground State Single-point/NPA configuration for VIP
        f.write("# M052X/6-31G+(d) Pop=NPA SCRF=(Solvent=Water)\n\n")
        f.write(f"Neutral Reference ground state calculation for {pdb_name}\n\n")
        
        # System Charge line dynamically formatted for the exact number of parsed fragments
        charge_line = f"{sys_charge} {sys_mult} " + " ".join([f"{c} {m}" for c, m in fragment_charges])
        f.write(charge_line + "\n")

        for i, res in enumerate(residues):
            f_id = i + 1
            atoms = res['atoms']

            if res['name'] == "K":
                for atom_name, xyz in atoms.items():
                    f.write(f"K(fragment={f_id}) {xyz[0]:12.8f} {xyz[1]:12.8f} {xyz[2]:12.8f}\n")
                continue

            if res['name'] == "G":
                nitro, anchor = atoms["N9"], atoms["C8"]
            elif res['name'] in ["DC", "DCP", "DC5", "DC3"]: 
                nitro, anchor = atoms["N1"], atoms["C2"]

            cme = atoms["C1'"]

            # Write existing base atoms (excluding C1')
            for name, xyz in atoms.items():
                if name == "C1'": continue
                f.write(f"{name[:1]}(fragment={f_id}) {xyz[0]:12.8f} {xyz[1]:12.8f} {xyz[2]:12.8f}\n")
            
            # Write Capping Methyl Carbon
            f.write(f"C(fragment={f_id}) {cme[0]:12.8f} {cme[1]:12.8f} {cme[2]:12.8f}\n")

            # Calculate and Write Methyl Hydrogens
            for j in range(3):
                h_pos = get_h_coords(anchor, nitro, cme, phiC[j], thetaC[j], blC[j])
                f.write(f"H(fragment={f_id}) {h_pos[0]:12.8f} {h_pos[1]:12.8f} {h_pos[2]:12.8f}\n")
        f.write("\n")

    print(f"Success: Main Gaussian input file generated: {base_name}.com")

    # --- Part 2: Optional Writing of Individual XYZ Files ---
    if args.xyz:
        print("Generating individual XYZ files for visual inspection...")
        for i, res in enumerate(residues):
            f_id = i + 1
            write_fragment_xyz(res, f_id, phiC, thetaC, blC, base_name)
            print(f" -> Generated XYZ for Fragment {f_id} ({res['name']}{res['seq']})")
    else:
        print("Skipping XYZ generation. (Use '--xyz' or '-x' flag to enable it next time)")

if __name__ == "__main__":
    main()