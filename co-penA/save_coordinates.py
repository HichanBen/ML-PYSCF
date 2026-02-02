#!/usr/bin/env python3
"""
Coordinate File Generator for Cu-Penicillamine Complex
Saves coordinates in multiple formats for different software
"""

import os

# Your coordinates from the original input
COORDS_DATA = """
C -6.48761 -4.97242 1.48230
C -7.93727 -4.98602 0.94199
H -8.63742 -5.45562 1.66343
H -7.98292 -5.57699 -0.00369
H -8.30189 -3.96377 0.71400
C -6.08494 -6.43841 1.75677
H -6.79668 -6.90212 2.47714
H -5.07090 -6.50674 2.20062
H -6.10603 -7.04331 0.82404
C -5.57115 -4.30739 0.40727
C -4.09466 -4.36089 0.79896
O -3.21209 -4.93698 -0.06507
H -2.26810 -4.96945 0.14337
O -3.69987 -3.91890 1.84571
N -6.01231 -2.91451 0.17743
H -5.67484 -4.85544 -0.55208
S -6.34630 -4.02050 3.07487
C -5.52051 0.92555 2.46932
C -5.36110 2.35886 1.89397
H -5.77444 2.43361 0.86449
H -4.29831 2.67674 1.86222
H -5.90928 3.09786 2.52328
C -4.92834 0.91160 3.88738
H -3.85195 1.19773 3.83542
H -4.99420 -0.10280 4.32871
H -5.41240 1.64820 4.55757
C -7.05510 0.54715 2.36663
C -7.93117 1.41086 3.25636
O -7.82656 1.34945 4.47366
N -7.33659 -0.86298 2.65689
H -7.35549 0.71830 1.30860
O -8.84776 2.23713 2.71054
H -9.41545 2.79163 3.26665
Cu -5.94420 -2.02580 1.91951
S -4.47973 -0.26743 1.51396
H -8.30389 -1.11923 2.35558
H -5.35116 -2.45939 -0.49326
"""

OUTPUT_DIR = "/mnt/user-data/outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_coords():
    """Parse coordinate string into atoms and positions"""
    atoms = []
    coords = []
    for line in COORDS_DATA.strip().split('\n'):
        parts = line.split()
        if len(parts) == 4:
            atoms.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
    return atoms, coords

def save_xyz(atoms, coords, filename="cu_penicillamine.xyz"):
    """Save in XYZ format (for VMD, Avogadro, ChimeraX)"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write("Cu-Penicillamine Complex (Initial Geometry)\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"{atom:3s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")
    return filepath

def save_pdb(atoms, coords, filename="cu_penicillamine.pdb"):
    """Save in PDB format (for PyMOL, VMD)"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, 'w') as f:
        f.write("REMARK   Cu-Penicillamine Complex\n")
        f.write("REMARK   Generated for DFT calculations\n")
        for i, (atom, coord) in enumerate(zip(atoms, coords), 1):
            # PDB format: ATOM, serial, name, resName, chainID, resSeq, x, y, z, occupancy, tempFactor, element
            f.write(f"ATOM  {i:5d}  {atom:3s} PEN A   1    "
                   f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                   f"  1.00  0.00          {atom:2s}\n")
        f.write("END\n")
    return filepath

def save_mol2(atoms, coords, filename="cu_penicillamine.mol2"):
    """Save in MOL2 format (for some docking programs)"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("Cu-Penicillamine\n")
        f.write(f"{len(atoms)} 0 0 0 0\n")
        f.write("SMALL\n")
        f.write("NO_CHARGES\n\n")
        
        f.write("@<TRIPOS>ATOM\n")
        for i, (atom, coord) in enumerate(zip(atoms, coords), 1):
            f.write(f"{i:7d} {atom:3s}{i:<4d} {coord[0]:10.4f} {coord[1]:10.4f} {coord[2]:10.4f} "
                   f"{atom:5s} 1 PEN     0.0000\n")
        
        f.write("@<TRIPOS>BOND\n")
        # Note: bonds would need to be defined here for complete MOL2
        # Simplified version without bond information
    return filepath

def save_gaussian(atoms, coords, filename="cu_penicillamine.gjf"):
    """Save in Gaussian input format"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, 'w') as f:
        f.write("%NProcShared=4\n")
        f.write("%Mem=4GB\n")
        f.write("%Chk=cu_penicillamine.chk\n")
        f.write("# B3LYP/def2SVP Opt Freq\n\n")
        f.write("Cu-Penicillamine Complex\n\n")
        f.write("0 2\n")  # Charge 0, Multiplicity 2 (doublet)
        for atom, coord in zip(atoms, coords):
            f.write(f"{atom:3s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")
        f.write("\n")
    return filepath

def save_turbomole(atoms, coords, filename="coord"):
    """Save in TURBOMOLE coord format"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    # Convert Angstrom to Bohr (TURBOMOLE uses atomic units)
    bohr = 1.889726125
    with open(filepath, 'w') as f:
        f.write("$coord\n")
        for atom, coord in zip(atoms, coords):
            x, y, z = [c * bohr for c in coord]
            f.write(f"{x:20.14f} {y:20.14f} {z:20.14f}  {atom.lower()}\n")
        f.write("$end\n")
    return filepath

def count_atoms(atoms):
    """Count atoms by element"""
    from collections import Counter
    counts = Counter(atoms)
    return counts

def main():
    print("\n" + "="*70)
    print("  Cu-Penicillamine Coordinate File Generator")
    print("="*70 + "\n")
    
    # Parse coordinates
    atoms, coords = parse_coords()
    
    print(f"Total atoms: {len(atoms)}\n")
    
    # Count atoms
    counts = count_atoms(atoms)
    print("Atomic composition:")
    for element, count in sorted(counts.items()):
        print(f"  {element}: {count}")
    
    # Derive molecular formula
    formula_parts = []
    for element in ['C', 'H', 'N', 'O', 'S', 'Cu']:
        if element in counts:
            count = counts[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
    formula = ''.join(formula_parts)
    print(f"\nMolecular formula: {formula}")
    
    print("\n" + "-"*70)
    print("Generating coordinate files...\n")
    
    # Save in different formats
    files_created = []
    
    xyz_file = save_xyz(atoms, coords)
    files_created.append(("XYZ", xyz_file, "VMD, Avogadro, ChimeraX, Jmol"))
    print(f"✓ XYZ format saved")
    
    pdb_file = save_pdb(atoms, coords)
    files_created.append(("PDB", pdb_file, "PyMOL, VMD, ChimeraX"))
    print(f"✓ PDB format saved")
    
    mol2_file = save_mol2(atoms, coords)
    files_created.append(("MOL2", mol2_file, "UCSF Chimera, some docking programs"))
    print(f"✓ MOL2 format saved")
    
    gjf_file = save_gaussian(atoms, coords)
    files_created.append(("Gaussian", gjf_file, "Gaussian, GaussView"))
    print(f"✓ Gaussian input saved")
    
    turbo_file = save_turbomole(atoms, coords)
    files_created.append(("TURBOMOLE", turbo_file, "TURBOMOLE"))
    print(f"✓ TURBOMOLE coord saved")
    
    print("\n" + "="*70)
    print("Summary of Generated Files:")
    print("="*70 + "\n")
    
    for format_name, filepath, software in files_created:
        print(f"{format_name:12s}: {filepath}")
        print(f"{'':12s}  Compatible with: {software}")
        print()
    
    print("All files saved to:", OUTPUT_DIR)
    print("\n" + "="*70)
    print("Files are ready for use in DFT calculations and visualization!")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
