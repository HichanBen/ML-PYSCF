#!/usr/bin/env python3
"""
Quick Validation Script for Cu-Penicillamine DFT Results
Checks geometry, energies, and convergence
"""

import os
import re
import numpy as np

OUTPUT_DIR = "/mnt/user-data/outputs"

def print_header(title):
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70)

def check_file_exists(filename):
    """Check if output file exists"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    exists = os.path.exists(filepath)
    status = "✓" if exists else "✗"
    print(f"  {status} {filename}")
    return exists

def read_xyz(filename):
    """Read XYZ file and return coordinates"""
    filepath = os.path.join(OUTPUT_DIR, filename)
    if not os.path.exists(filepath):
        return None, None
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    natoms = int(lines[0].strip())
    atoms = []
    coords = []
    
    for line in lines[2:2+natoms]:
        parts = line.split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    
    return atoms, np.array(coords)

def calculate_bond_length(coords, i, j):
    """Calculate distance between atoms i and j"""
    return np.linalg.norm(coords[i] - coords[j])

def validate_geometry(atoms, coords):
    """Validate key geometric parameters"""
    print_header("Geometry Validation")
    
    # Find Cu, S, N, O atoms
    cu_idx = [i for i, a in enumerate(atoms) if a == 'Cu']
    s_idx = [i for i, a in enumerate(atoms) if a == 'S']
    n_idx = [i for i, a in enumerate(atoms) if a == 'N']
    o_idx = [i for i, a in enumerate(atoms) if a == 'O']
    
    print(f"\nAtom counts:")
    print(f"  Cu: {len(cu_idx)} (expected: 1)")
    print(f"  S:  {len(s_idx)} (expected: 2)")
    print(f"  N:  {len(n_idx)} (expected: 2)")
    print(f"  O:  {len(o_idx)} (expected: 5)")
    
    if len(cu_idx) != 1:
        print("  ⚠ Warning: Expected exactly 1 Cu atom!")
        return
    
    cu = cu_idx[0]
    
    print(f"\nKey bond lengths (Å):")
    
    # Check Cu-S bonds
    for s in s_idx:
        dist = calculate_bond_length(coords, cu, s)
        status = "✓" if 2.0 < dist < 2.5 else "⚠"
        print(f"  {status} Cu{cu}-S{s}: {dist:.3f} Å (expected: 2.10-2.30)")
    
    # Check Cu-N bonds
    for n in n_idx:
        dist = calculate_bond_length(coords, cu, n)
        if dist < 3.0:  # Only report if reasonable bonding distance
            status = "✓" if 1.8 < dist < 2.2 else "⚠"
            print(f"  {status} Cu{cu}-N{n}: {dist:.3f} Å (expected: 1.90-2.10)")
    
    # Calculate coordination number
    coord_atoms = []
    for i, atom in enumerate(atoms):
        if i != cu:
            dist = calculate_bond_length(coords, cu, i)
            if dist < 2.5:  # Bonding distance cutoff
                coord_atoms.append((atom, i, dist))
    
    print(f"\nCu coordination environment:")
    print(f"  Coordination number: {len(coord_atoms)}")
    for atom, idx, dist in sorted(coord_atoms, key=lambda x: x[2]):
        print(f"    {atom}{idx}: {dist:.3f} Å")

def compare_geometries():
    """Compare initial and optimized geometries"""
    print_header("Geometry Comparison")
    
    atoms_i, coords_i = read_xyz("initial_geometry.xyz")
    atoms_o, coords_o = read_xyz("optimized_geometry.xyz")
    
    if atoms_i is None or atoms_o is None:
        print("  ✗ Could not read geometry files")
        return
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean((coords_i - coords_o)**2))
    print(f"\nRMSD between initial and optimized: {rmsd:.4f} Å")
    
    if rmsd < 0.01:
        print("  ⚠ Warning: Very small RMSD - geometry barely changed")
    elif rmsd < 0.5:
        print("  ✓ Good: Moderate geometry relaxation")
    elif rmsd < 2.0:
        print("  ✓ Significant geometry change (normal for first optimization)")
    else:
        print("  ⚠ Large geometry change - verify structure is reasonable")
    
    # Find maximum displacement
    displacements = np.linalg.norm(coords_i - coords_o, axis=1)
    max_disp_idx = np.argmax(displacements)
    max_disp = displacements[max_disp_idx]
    
    print(f"\nMaximum atomic displacement:")
    print(f"  Atom {max_disp_idx} ({atoms_i[max_disp_idx]}): {max_disp:.4f} Å")

def check_fukui_functions():
    """Analyze Fukui function results"""
    print_header("Fukui Function Analysis")
    
    fukui_file = os.path.join(OUTPUT_DIR, "fukui_analysis.txt")
    if not os.path.exists(fukui_file):
        print("  ✗ Fukui analysis file not found")
        return
    
    with open(fukui_file, 'r') as f:
        lines = f.readlines()
    
    # Parse Fukui data
    data = []
    for line in lines:
        if line.strip() and not line.startswith(('=', 'Atom', '-', 'INTER')):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    idx = int(parts[0])
                    symbol = parts[1]
                    f_plus = float(parts[2])
                    f_minus = float(parts[3])
                    f_zero = float(parts[4])
                    dual = float(parts[5])
                    data.append((idx, symbol, f_plus, f_minus, f_zero, dual))
                except:
                    pass
    
    if not data:
        print("  ✗ Could not parse Fukui data")
        return
    
    # Find most reactive sites
    data_sorted_plus = sorted(data, key=lambda x: abs(x[2]), reverse=True)
    data_sorted_minus = sorted(data, key=lambda x: abs(x[3]), reverse=True)
    data_sorted_dual = sorted(data, key=lambda x: abs(x[5]), reverse=True)
    
    print("\nTop 5 electrophilic sites (highest f-):")
    for i, (idx, symbol, fp, fm, fz, dual) in enumerate(data_sorted_minus[:5], 1):
        print(f"  {i}. Atom {idx:3d} ({symbol:2s}): f- = {fm:8.5f}")
    
    print("\nTop 5 nucleophilic sites (highest f+):")
    for i, (idx, symbol, fp, fm, fz, dual) in enumerate(data_sorted_plus[:5], 1):
        print(f"  {i}. Atom {idx:3d} ({symbol:2s}): f+ = {fp:8.5f}")
    
    print("\nTop 5 sites by dual descriptor |Δf|:")
    for i, (idx, symbol, fp, fm, fz, dual) in enumerate(data_sorted_dual[:5], 1):
        character = "electrophilic" if dual > 0 else "nucleophilic"
        print(f"  {i}. Atom {idx:3d} ({symbol:2s}): Δf = {dual:8.5f} ({character})")

def check_summary():
    """Read and display calculation summary"""
    print_header("Calculation Summary")
    
    summary_file = os.path.join(OUTPUT_DIR, "calculation_summary.txt")
    if not os.path.exists(summary_file):
        print("  ✗ Summary file not found")
        return
    
    with open(summary_file, 'r') as f:
        content = f.read()
    
    # Extract key information
    energy_match = re.search(r'Electronic energy:\s+([-\d.]+)\s+Hartree', content)
    energy_ev_match = re.search(r'Electronic energy:\s+([-\d.]+)\s+eV', content)
    converged_match = re.search(r'SCF converged:\s+(\w+)', content)
    gap_match = re.search(r'HOMO-LUMO gap:\s+([\d.]+)\s+eV', content)
    
    if energy_match:
        print(f"\nElectronic Energy: {energy_match.group(1)} Hartree")
    if energy_ev_match:
        print(f"                   {energy_ev_match.group(1)} eV")
    
    if converged_match:
        status = converged_match.group(1)
        symbol = "✓" if status.lower() == "yes" else "✗"
        print(f"\nSCF Convergence: {symbol} {status}")
    
    if gap_match:
        gap = float(gap_match.group(1))
        print(f"\nHOMO-LUMO Gap: {gap:.4f} eV")
        if gap < 1.0:
            print("  → Very reactive system")
        elif gap < 3.0:
            print("  → Moderately reactive (typical for Cu complexes)")
        elif gap < 5.0:
            print("  → Relatively stable")
        else:
            print("  → Very stable system")

def main():
    """Main validation workflow"""
    
    print("\n" + "="*70)
    print("  Cu-PENICILLAMINE DFT CALCULATION VALIDATION")
    print("="*70)
    
    # Check output files
    print_header("Output Files Check")
    
    files_to_check = [
        "initial_geometry.xyz",
        "optimized_geometry.xyz",
        "density.cube",
        "mep.cube",
        "homo_alpha.cube",
        "fukui_analysis.txt",
        "calculation_summary.txt"
    ]
    
    all_exist = True
    for filename in files_to_check:
        if not check_file_exists(filename):
            all_exist = False
    
    if not all_exist:
        print("\n⚠ Some output files are missing!")
        print("  Make sure the calculation completed successfully.")
        return
    
    # Read and validate optimized geometry
    atoms, coords = read_xyz("optimized_geometry.xyz")
    if atoms is not None:
        validate_geometry(atoms, coords)
    
    # Compare geometries
    compare_geometries()
    
    # Check Fukui functions
    check_fukui_functions()
    
    # Display summary
    check_summary()
    
    # Final recommendations
    print_header("Recommendations")
    
    print("""
Next Steps:

1. VISUALIZE RESULTS:
   - Load .cube files in VMD, Avogadro, or ChimeraX
   - Examine MEP for reactive sites
   - View HOMO/LUMO orbitals
   
2. VERIFY CHEMISTRY:
   - Check if Cu coordination makes sense
   - Verify bond lengths are reasonable
   - Confirm charge distribution is logical
   
3. PRODUCTION CALCULATIONS (if needed):
   - Use ORCA for full D3(BJ) dispersion
   - Calculate vibrational frequencies
   - Add solvent effects (COSMO/PCM)
   - Perform NBO or QTAIM analysis
   
4. COMPARE WITH LITERATURE:
   - Similar Cu-thiol complexes
   - Known antimalarial copper complexes
   - Fe(III)PPIX binding studies

For detailed interpretation, see: GUIDE_CuPenicillamine_DFT.md
""")
    
    print("="*70)
    print("  Validation Complete!")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
