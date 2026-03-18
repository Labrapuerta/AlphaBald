# Problem 3 Solutions Guide

This document outlines how to solve each question in Problem_3 using the AlphaBald pipeline functions.

## Overview of Problem 3

You have a modeled protein structure (given as CA atoms in the notebook) that needs to be assessed and potentially fixed.

---

## Question-by-Question Solutions

### 1. What's the protein fold and family?

**Available Function:** `src.Analysis.identify_protein_family()`

```python
from src.Analysis import identify_protein_family, ModelAssessor

# Option 1: Direct function call
result = identify_protein_family(pdb_file="model.pdb")
print(f"Family: {result['family_name']}")
print(f"Fold: {result['fold']}")

# Option 2: Using ModelAssessor
assessor = ModelAssessor("model.pdb")
family = assessor.identify_family()
```

**Manual Steps:**
1. Extract sequence from PDB
2. Search Pfam: `hmmscan databases/hmm/Pfam/Pfam-A.hmm sequence.fa`
3. Check SCOP/CATH databases for fold classification

---

### 2. Obtain HMM profile and alignment (p3b.hmm, p3b.aln)

**Available Functions:**
- `src.Analysis.create_hmm_profile()`
- `src.Analysis.align_hmm_to_sequence()`

```python
from src.Analysis import create_hmm_profile, align_hmm_to_sequence

# First, create MSA with homologs
# (Use the pipeline to find homologs first)

# Create HMM from alignment
hmm_path = create_hmm_profile(
    alignment_file="Alignments/aligned_homologs.fa",
    output_file="Problems/Problem_3/p3b.hmm",
    name="protein_family"
)

# Align HMM to closest homolog
alignment = align_hmm_to_sequence(
    hmm_file="Problems/Problem_3/p3b.hmm",
    sequence_file="closest_homolog.fa",
    output_file="Problems/Problem_3/p3b.aln"
)
```

**Commands:**
```bash
# Build HMM from alignment
hmmbuild p3b.hmm alignment.fa

# Align sequence to HMM
hmmalign -o p3b.aln p3b.hmm homolog.fa
```

---

### 3. Add a cation (p3c.pdb)

**Available Function:** `src.Analysis.add_cation_to_structure()`

```python
from src.Analysis import add_cation_to_structure

# Add calcium (common for many enzymes)
output = add_cation_to_structure(
    pdb_file="model.pdb",
    cation="CA",  # Calcium
    output_file="Problems/Problem_3/p3c.pdb"
)

# For zinc-binding proteins:
# cation="ZN"

# For magnesium-dependent enzymes:
# cation="MG"
```

**Note:** The function tries to identify the binding site automatically based on coordinating residues (Asp, Glu, His, Cys).

---

### 4. Mark important functional residues with @ (in p3b.aln)

**Available Function:** `src.Analysis.identify_functional_residues()`

```python
from src.Analysis import identify_functional_residues

# Identify and mark conserved residues
important = identify_functional_residues(
    alignment_file="Problems/Problem_3/p3b.aln",
    output_file="Problems/Problem_3/p3b_marked.aln",
    conservation_threshold=0.9
)

for res in important:
    print(f"Position {res['position']}: {res['residue']} (conservation: {res['conservation']:.0%})")
```

---

### 5. Active site analysis (p3e_fig1)

**Available Function:** `src.Analysis.analyze_active_site()`

```python
from src.Analysis import analyze_active_site

# Analyze active site
analysis = analyze_active_site(
    pdb_file="model.pdb",
    active_site_residues=[123, 156, 189],  # From functional residue analysis
    reference_pdb="reference_structure.pdb"  # Optional
)

print(f"Active: {analysis['is_active_prediction']}")
print(f"Preserved: {analysis['is_preserved']}")
```

**For visualization (p3e_fig1):**
- Use PyMOL or Chimera to visualize active site
- Color conserved residues green, non-conserved red
- Export as PNG

---

### 6. Validate model regions (p3f_fig1, p3f_fig2, p3f_fig3)

**Available Function:** `src.Analysis.validate_model_regions()`

```python
from src.Analysis import validate_model_regions

# Run validation
validation = validate_model_regions(
    pdb_file="model.pdb",
    method="all"  # Uses ProSA and DOPE
)

print("Problematic regions:")
for region in validation['problematic_regions']:
    print(f"  Residue {region['residue']}: {region['issue']}")
```

**Manual Steps:**
1. Submit to ProSA web server: https://prosa.services.came.sbg.ac.at/prosa.php
2. Use multiple window sizes (10, 40) for detailed analysis
3. Save energy profile images

---

### 7. Calculate secondary structure with DSSP (p3b.dssp)

**Available Function:** `src.Analysis.run_dssp_analysis()`

```python
from src.Analysis import run_dssp_analysis

# Run DSSP
dssp_result = run_dssp_analysis(
    pdb_file="model.pdb",
    output_file="Problems/Problem_3/p3b.dssp"
)

# Summary of secondary structure
print("Secondary structure composition:")
for ss, count in dssp_result['ss_summary'].items():
    print(f"  {ss}: {count} residues")
```

**Command:**
```bash
dssp -i model.pdb -o p3b.dssp
```

---

### 8. Fix the problem (p3h.pdb with cation)

**Available Function:** `src.Analysis.fix_model_problems()`

```python
from src.Analysis import ModelAssessor

# Full workflow
assessor = ModelAssessor("model.pdb", output_dir="Problems/Problem_3")

# Run DSSP
assessor.run_dssp("p3b.dssp")

# Find problems
problems = assessor.find_problematic_regions()

# Fix problems (generates MODELLER script)
fixed_path = assessor.fix_problems(
    problems['problematic_regions'],
    output_file="p3h.pdb"
)

# Add cation to fixed model
assessor.pdb_file = fixed_path
assessor.add_cation("CA", "p3h.pdb")
```

**MODELLER Loop Refinement:**
```python
from src.modeller import generate_loop_refinement_script

script = generate_loop_refinement_script(
    pdb_file="model.pdb",
    target_id="target",
    loop_residues=[(120, 130), (200, 210)],  # Problem regions
    num_models=5
)

# Save and run with MODELLER
with open("fix_loops.py", "w") as f:
    f.write(script)
```

---

## Complete Workflow Example

```python
from src.Analysis import ModelAssessor
from src.modeller import create_modeller_scripts_for_pipeline

# Step 1: Create assessor
assessor = ModelAssessor("model.pdb", output_dir="Problems/Problem_3")

# Step 2: Identify family
family = assessor.identify_family()
print(f"Protein family: {family.get('family_name', 'Unknown')}")

# Step 3: Run DSSP
dssp = assessor.run_dssp("p3b.dssp")

# Step 4: Validate model
validation = assessor.find_problematic_regions()

# Step 5: Add cation
assessor.add_cation("CA", "p3c.pdb")

# Step 6: Fix problems
if validation['problematic_regions']:
    assessor.fix_problems(output_file="p3h.pdb")
```

---

## What's Still Needed (External Tools)

1. **ProSA Analysis:** Must be done via web server or local installation
2. ~~**Visualization:** Use PyMOL/Chimera for figure generation~~ **Now automated with `src.Analysis.visualization`**
3. ~~**MODELLER:** Required for loop refinement~~ **Scripts auto-generated, just run them**
4. **Manual Curation:** Some steps require expert judgment (e.g., choosing the right cation)

---

## PyMOL Visualization (NEW)

Since PyMOL is installed in your conda environment, you can now auto-generate all figures:

```python
from src.Analysis import StructureVisualizer

# Create visualizer
viz = StructureVisualizer("model.pdb", output_dir="Problems/Problem_3")

# Generate active site figure (p3e_fig1)
viz.active_site(
    residues=[123, 156, 189],
    output_name="p3e_fig1.png",
    conserved=[123, 156],      # Green
    non_conserved=[189]        # Red
)

# Generate problematic regions figures (p3f_fig1-3)
viz.problematic_regions(
    regions=[(120, 130), (200, 210)],
    output_name="p3f_fig1.png"
)

# Generate secondary structure figure (p3g_fig1)
viz.secondary_structure("p3g_fig1.png")

# Generate cation binding figure
viz.cation_site("p3h_fig1.png", cation="CA", binding_residues=[123, 156])
```

**Individual functions also available:**
```python
from src.Analysis import (
    visualize_active_site,
    visualize_problematic_regions,
    visualize_secondary_structure,
    visualize_cation_binding,
    compare_structures
)

# Compare original and fixed structures
compare_structures("model.pdb", "p3h.pdb", "p3h_fig2.png", align=True)
```

---

## Summary of Implemented Functions

| Question | Function | Status |
|----------|----------|--------|
| Fold/Family | `identify_protein_family()` | Implemented |
| HMM Profile | `create_hmm_profile()` | Implemented |
| HMM Alignment | `align_hmm_to_sequence()` | Implemented |
| Add Cation | `add_cation_to_structure()` | Implemented |
| Functional Residues | `identify_functional_residues()` | Implemented |
| Active Site | `analyze_active_site()` | Implemented |
| Validation | `validate_model_regions()` | Implemented |
| DSSP | `run_dssp_analysis()` | Implemented |
| Fix Problems | `fix_model_problems()` | Implemented |
