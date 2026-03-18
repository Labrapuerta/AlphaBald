"""
Analysis module for protein structure assessment.

Provides tools for:
- Protein fold/family identification
- HMM profile creation and searching
- Cation placement
- Active site analysis
- Secondary structure analysis
- Model quality assessment
- PyMOL visualization
"""

from .assessment import (
    identify_protein_family,
    create_hmm_profile,
    align_hmm_to_sequence,
    add_cation_to_structure,
    identify_functional_residues,
    analyze_active_site,
    validate_model_regions,
    run_dssp_analysis,
    fix_model_problems,
    ModelAssessor
)

from .visualization import (
    visualize_active_site,
    visualize_problematic_regions,
    visualize_secondary_structure,
    visualize_cation_binding,
    compare_structures,
    StructureVisualizer,
    run_pymol_script
)

__all__ = [
    # Assessment functions
    'identify_protein_family',
    'create_hmm_profile',
    'align_hmm_to_sequence',
    'add_cation_to_structure',
    'identify_functional_residues',
    'analyze_active_site',
    'validate_model_regions',
    'run_dssp_analysis',
    'fix_model_problems',
    'ModelAssessor',
    # Visualization functions
    'visualize_active_site',
    'visualize_problematic_regions',
    'visualize_secondary_structure',
    'visualize_cation_binding',
    'compare_structures',
    'StructureVisualizer',
    'run_pymol_script'
]
