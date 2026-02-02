# Cu-Penicillamine DFT Analysis Package

Complete toolkit for DFT analysis of copper-penicillamine complexes, based on methodologies from Sun et al. (2017) and Nkungli & Ghogomu (2017).

## 📦 Package Contents

### Main Scripts
1. **cu_penicillamine_analysis.py** - Main PySCF calculation script
2. **validate_results.py** - Results validation and checking
3. **cu_penicillamine_orca.inp** - ORCA input for production calculations

### Documentation
4. **GUIDE_CuPenicillamine_DFT.md** - Complete user guide
5. **README.md** - This file

---

## 🚀 Quick Start

### Option 1: PySCF (Recommended for Learning)

```bash
# Install dependencies
pip install pyscf numpy scipy

# Run calculation
python cu_penicillamine_analysis.py

# Validate results
python validate_results.py
```

**Time:** ~30-60 minutes  
**Outputs:** See `/mnt/user-data/outputs/`

### Option 2: ORCA (Recommended for Production)

```bash
# Run ORCA
orca cu_penicillamine_orca.inp > cu_penicillamine.out &

# Monitor
tail -f cu_penicillamine.out

# Validate
grep "FINAL SINGLE POINT ENERGY" cu_penicillamine.out
```

**Time:** ~10-30 minutes (4 cores)  
**Outputs:** Current directory + .cube files

---

## 📊 What Gets Calculated

| Property | Method | Output File |
|----------|--------|-------------|
| **Optimized Geometry** | BP86-D3/def2-SVP | optimized_geometry.xyz |
| **Electronic Energy** | B3LYP-D3/def2-TZVP | calculation_summary.txt |
| **Electron Density** | B3LYP | density.cube |
| **MEP** | B3LYP | mep.cube |
| **HOMO/LUMO** | B3LYP | homo_alpha.cube, lumo_alpha.cube |
| **Fukui Functions** | B3LYP | fukui_analysis.txt |

---

## 📁 Expected Outputs

After running successfully:

```
/mnt/user-data/outputs/
├── initial_geometry.xyz          # Input structure
├── optimized_geometry.xyz        # Relaxed structure
├── density.cube                  # Electron density
├── mep.cube                      # Molecular electrostatic potential
├── homo_alpha.cube               # HOMO orbital
├── lumo_alpha.cube               # LUMO orbital
├── fukui_analysis.txt            # Local reactivity descriptors
└── calculation_summary.txt       # Complete summary
```

---

## ✅ Validation Checklist

Run `python validate_results.py` to check:

- [x] All output files generated
- [x] SCF converged
- [x] Geometry optimization completed
- [x] Reasonable bond lengths (Cu-S: 2.1-2.3 Å, Cu-N: 1.9-2.1 Å)
- [x] Fukui functions calculated
- [x] HOMO-LUMO gap reasonable (2-4 eV for Cu complexes)

---

## 🔍 Interpreting Results

### 1. Geometry
- Check Cu coordination (square planar/tetrahedral?)
- Verify bond lengths against literature
- RMSD between initial/optimized should be 0.1-2.0 Å

### 2. MEP (Molecular Electrostatic Potential)
- **Blue regions:** Positive potential → Nucleophilic attack sites
- **Red regions:** Negative potential → Electrophilic attack sites
- Expected: Cu (blue), S/O (red)

### 3. Fukui Functions
- **f+ (nucleophilic attack):** High on sites that donate electrons
- **f- (electrophilic attack):** High on sites that accept electrons
- **Δf > 0:** Electrophilic character
- **Δf < 0:** Nucleophilic character

### 4. HOMO-LUMO
- **Gap < 3 eV:** Reactive complex
- **Gap 3-5 eV:** Moderately stable
- **Gap > 5 eV:** Very stable
- Cu d-orbitals usually in HOMO

---

## 🎨 Visualization

### Using VMD
```bash
vmd density.cube
vmd mep.cube
vmd homo_alpha.cube
```

### Using Avogadro
```bash
avogadro optimized_geometry.xyz
# Then: Extensions → Surfaces → Isosurface
```

### Using ChimeraX
```bash
chimerax density.cube
# Then: Volume → Surface
```

---

## 🔧 Troubleshooting

### Problem: "SCF did not converge"
**Solution:**
```python
# In cu_penicillamine_analysis.py, adjust:
mf.max_cycle = 200
mf.conv_tol = 1e-5
mf.diis_space = 12
```

### Problem: "Geometry optimization failed"
**Solution:**
- Check initial geometry (bond lengths reasonable?)
- Try smaller basis: `basis='def2-sv(p)'`
- Reduce convergence criteria
- Use ORCA instead (more robust)

### Problem: "Memory error"
**Solution:**
- Reduce grid density: `nx=40, ny=40, nz=40`
- Use smaller basis set
- Increase system memory

### Problem: "Fukui functions all near zero"
**Solution:**
- Verify N+1 and N-1 calculations converged
- Check spin multiplicities correct
- Try Hirshfeld charges instead of Mulliken

---

## 📚 References

### Methodology Papers

1. **Sun, Q., et al. (2017)**  
   "PySCF: the Python-based simulations of chemistry framework"  
   *WIREs Comput Mol Sci*, e1340  
   https://doi.org/10.1002/wcms.1340

2. **Nkungli, N. K., & Ghogomu, J. N. (2017)**  
   "Theoretical analysis of the binding of iron(III) protoporphyrin IX to 
   4-methoxyacetophenone thiosemicarbazone via DFT-D3, MEP, QTAIM, NCI, 
   ELF, and LOL studies"  
   *J Mol Model*, 23:200  
   https://doi.org/10.1007/s00894-017-3370-4

### Software

- **PySCF:** https://pyscf.org/
- **ORCA:** https://www.faccts.de/orca/
- **VMD:** https://www.ks.uiuc.edu/Research/vmd/

---

## 📧 Support

For questions about:
- **PySCF usage:** https://github.com/pyscf/pyscf/issues
- **ORCA usage:** https://orcaforum.kofo.mpg.de
- **This package:** Check GUIDE_CuPenicillamine_DFT.md

---

## 🎯 Next Steps

1. **Run the calculation** (choose PySCF or ORCA)
2. **Validate results** using validate_results.py
3. **Visualize outputs** (VMD/Avogadro/ChimeraX)
4. **Interpret chemistry** (see Guide)
5. **Compare with literature**

For antimalarial activity analysis:
- Compare with Fe(III)PPIX binding
- Identify coordination sites
- Assess reactivity patterns
- Correlate with biological activity

---

## 📄 License

Educational and research use. Based on open-source software (PySCF) and 
published methodologies.

---

## 🔄 Version

**Version:** 1.0  
**Created:** February 2026  
**Based on:** PySCF 2.x, ORCA 5.x

---

**Happy Computing! 🧪⚛️**

For detailed methodology and interpretation, see **GUIDE_CuPenicillamine_DFT.md**
