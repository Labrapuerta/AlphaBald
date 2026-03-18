"""
Structure superimposition utilities for homology modeling.

Provides functions to align and superimpose protein structures based on
sequence alignments.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio import SeqIO
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom


def get_ca_atoms(structure: Structure, chain_id: str) -> List[Atom]:
    """
    Get all C-alpha atoms from a specific chain.

    Parameters
    ----------
    structure : Structure
        BioPython structure object
    chain_id : str
        Chain identifier (e.g., 'A', 'B')

    Returns
    -------
    list
        List of C-alpha Atom objects
    """
    chain = structure[0][chain_id]
    return [
        res['CA'] for res in chain.get_residues()
        if res.has_id('CA') and res.id[0] == ' '  # Standard residues only
    ]


def parse_alignment(alignment_file: Path) -> Dict[str, str]:
    """
    Parse a FASTA alignment file.

    Parameters
    ----------
    alignment_file : Path
        Path to aligned FASTA file

    Returns
    -------
    dict
        Mapping of sequence ID to aligned sequence string
    """
    alignment = {}
    for record in SeqIO.parse(str(alignment_file), "fasta"):
        alignment[record.id] = str(record.seq)
    return alignment


def get_aligned_positions(
    target_aligned: str,
    template_aligned: str
) -> List[Tuple[int, int]]:
    """
    Get corresponding residue positions from aligned sequences.

    Only returns positions where both target and template have residues
    (neither is a gap).

    Parameters
    ----------
    target_aligned : str
        Aligned target sequence (with gaps as '-')
    template_aligned : str
        Aligned template sequence (with gaps as '-')

    Returns
    -------
    list
        List of (target_pos, template_pos) tuples (0-indexed)
    """
    positions = []
    target_pos = 0
    template_pos = 0

    for t_char, p_char in zip(target_aligned, template_aligned):
        if t_char != '-' and p_char != '-':
            positions.append((target_pos, template_pos))

        if t_char != '-':
            target_pos += 1
        if p_char != '-':
            template_pos += 1

    return positions


def get_ca_mapping(
    pdb_id: str,
    alignment_file: Path = Path("Alignments/aligned_templates.fasta"),
    templates_dir: Path = Path("Templates"),
    target_id: Optional[str] = None
) -> Tuple[Dict[int, Atom], Structure]:
    """
    Create mapping of target residue positions to template C-alpha atoms.

    Parameters
    ----------
    pdb_id : str
        Template PDB ID with chain (e.g., "1ABC_A")
    alignment_file : Path
        Path to aligned FASTA file containing target and template
    templates_dir : Path
        Directory containing PDB files
    target_id : str, optional
        Target sequence ID in alignment. If None, uses first sequence.

    Returns
    -------
    tuple
        (mapping dict of target_pos -> CA atom, Structure object)
    """
    # Parse PDB ID
    pdb_4 = pdb_id[:4].lower()
    chain_id = pdb_id[5] if len(pdb_id) > 5 else "A"

    # Load structure
    pdb_file = templates_dir / f"pdb{pdb_4}.ent"
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_4, str(pdb_file))

    # Get CA atoms from template
    ca_atoms = get_ca_atoms(structure, chain_id)

    # Parse alignment
    alignment = parse_alignment(alignment_file)

    # Find target and template in alignment
    if target_id is None:
        target_id = list(alignment.keys())[0]

    if target_id not in alignment:
        raise ValueError(f"Target '{target_id}' not found in alignment")

    if pdb_id not in alignment:
        raise ValueError(f"Template '{pdb_id}' not found in alignment")

    aligned_target = alignment[target_id]
    aligned_template = alignment[pdb_id]

    # Get aligned positions
    positions = get_aligned_positions(aligned_target, aligned_template)

    # Create mapping
    mapping = {}
    for target_pos, template_pos in positions:
        if template_pos < len(ca_atoms):
            mapping[target_pos] = ca_atoms[template_pos]

    return mapping, structure


def superimpose_structures(
    fixed_structure: Structure,
    moving_structure: Structure,
    fixed_atoms: List[Atom],
    moving_atoms: List[Atom]
) -> float:
    """
    Superimpose moving structure onto fixed structure using C-alpha atoms.

    Parameters
    ----------
    fixed_structure : Structure
        Reference structure (will not be modified)
    moving_structure : Structure
        Structure to be superimposed (will be rotated/translated)
    fixed_atoms : list
        List of C-alpha atoms from fixed structure
    moving_atoms : list
        List of corresponding C-alpha atoms from moving structure

    Returns
    -------
    float
        RMSD value after superimposition
    """
    if len(fixed_atoms) != len(moving_atoms):
        raise ValueError(
            f"Atom lists must have same length: {len(fixed_atoms)} vs {len(moving_atoms)}"
        )

    if len(fixed_atoms) < 3:
        raise ValueError("Need at least 3 atoms for superimposition")

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)

    # Apply rotation/translation to all atoms in moving structure
    sup.apply(moving_structure.get_atoms())

    return sup.rms


def superimpose_template_to_target(
    template_pdb_id: str,
    target_structure: Structure,
    alignment_file: Path,
    templates_dir: Path = Path("Templates"),
    target_id: Optional[str] = None
) -> Tuple[Structure, float]:
    """
    Superimpose a template structure onto a target structure.

    Parameters
    ----------
    template_pdb_id : str
        Template PDB ID with chain (e.g., "1ABC_A")
    target_structure : Structure
        Target structure to superimpose onto
    alignment_file : Path
        Path to sequence alignment file
    templates_dir : Path
        Directory containing template PDB files
    target_id : str, optional
        Target sequence ID in alignment

    Returns
    -------
    tuple
        (superimposed template Structure, RMSD)
    """
    # Get CA mapping for template
    template_mapping, template_structure = get_ca_mapping(
        template_pdb_id, alignment_file, templates_dir, target_id
    )

    # Get target CA atoms at aligned positions
    target_chain_id = list(target_structure[0].get_chains())[0].id
    target_ca = get_ca_atoms(target_structure, target_chain_id)

    # Build atom lists
    fixed_atoms = []
    moving_atoms = []

    for target_pos, template_atom in template_mapping.items():
        if target_pos < len(target_ca):
            fixed_atoms.append(target_ca[target_pos])
            moving_atoms.append(template_atom)

    # Superimpose
    rmsd = superimpose_structures(
        target_structure, template_structure,
        fixed_atoms, moving_atoms
    )

    return template_structure, rmsd


def write_structure(structure: Structure, output_path: Path):
    """
    Write a structure to a PDB file.

    Parameters
    ----------
    structure : Structure
        BioPython structure object
    output_path : Path
        Output file path
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))


def get_match_mismatch_positions(
    target_aligned: str,
    template_aligned: str
) -> Tuple[List[int], List[int], List[int]]:
    """
    Get positions of matches, mismatches, and gaps from alignment.

    Parameters
    ----------
    target_aligned : str
        Aligned target sequence
    template_aligned : str
        Aligned template sequence

    Returns
    -------
    tuple
        (match_positions, mismatch_positions, gap_positions) - all 1-indexed for template
    """
    matches = []
    mismatches = []
    gaps = []

    template_pos = 0
    for t_char, p_char in zip(target_aligned, template_aligned):
        if p_char != '-':
            template_pos += 1
            if t_char == '-':
                gaps.append(template_pos)
            elif t_char == p_char:
                matches.append(template_pos)
            else:
                mismatches.append(template_pos)

    return matches, mismatches, gaps


class SuperimpositionVisualizer:
    """
    Visualize superimposed structures using py3Dmol with MSA-based coloring.

    Parameters
    ----------
    alignment_file : Path or str
        Path to aligned FASTA file
    templates_dir : Path or str, optional
        Directory containing PDB files. Default "Templates"
    output_dir : Path or str, optional
        Directory to save superimposed structures. Default "Superimposed"
    """

    # Color scheme for visualization
    COLORS = {
        'match': '0x00FF00',       # Green for matches
        'mismatch': '0xFF0000',    # Red for mismatches
        'gap': '0x808080',         # Gray for gaps
        'target': '0x0000FF',      # Blue for target
    }

    def __init__(
        self,
        alignment_file: str = "Alignments/aligned_templates.fasta",
        templates_dir: str = "Templates",
        output_dir: str = "Superimposed"
    ):
        self.alignment_file = Path(alignment_file)
        self.templates_dir = Path(templates_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Parse alignment
        self.alignment = parse_alignment(self.alignment_file)
        self.target_id = list(self.alignment.keys())[0]
        self.target_aligned = self.alignment[self.target_id]

        # Store superimposed structures
        self.superimposed = {}

    def superimpose_and_save(
        self,
        template_ids: List[str],
        reference_template: Optional[str] = None
    ) -> Dict[str, Tuple[Path, float]]:
        """
        Superimpose templates and save to output directory.

        Parameters
        ----------
        template_ids : list
            List of template PDB IDs (e.g., ["1ABC_A", "2DEF_B"])
        reference_template : str, optional
            Template to use as reference. If None, uses first template.

        Returns
        -------
        dict
            Mapping of template_id to (output_path, rmsd)
        """
        parser = PDBParser(QUIET=True)
        results = {}

        # Load reference structure
        if reference_template is None:
            reference_template = template_ids[0]

        ref_pdb_4 = reference_template[:4].lower()
        ref_chain = reference_template[5] if len(reference_template) > 5 else "A"
        ref_file = self.templates_dir / f"pdb{ref_pdb_4}.ent"

        if not ref_file.exists():
            print(f"Reference PDB not found: {ref_file}")
            return results

        ref_structure = parser.get_structure(ref_pdb_4, str(ref_file))
        ref_ca = get_ca_atoms(ref_structure, ref_chain)

        # Save reference
        ref_output = self.output_dir / f"{reference_template}_ref.pdb"
        write_structure(ref_structure, ref_output)
        results[reference_template] = (ref_output, 0.0)

        # Superimpose each template to reference
        for template_id in template_ids:
            if template_id == reference_template:
                continue

            if template_id not in self.alignment:
                print(f"Template {template_id} not in alignment, skipping")
                continue

            try:
                pdb_4 = template_id[:4].lower()
                chain_id = template_id[5] if len(template_id) > 5 else "A"
                pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"

                if not pdb_file.exists():
                    continue

                structure = parser.get_structure(pdb_4, str(pdb_file))
                ca_atoms = get_ca_atoms(structure, chain_id)

                # Get aligned positions
                template_aligned = self.alignment[template_id]
                ref_aligned = self.alignment[reference_template]
                positions = get_aligned_positions(ref_aligned, template_aligned)

                # Build atom pairs
                fixed_atoms = []
                moving_atoms = []
                for ref_pos, temp_pos in positions:
                    if ref_pos < len(ref_ca) and temp_pos < len(ca_atoms):
                        fixed_atoms.append(ref_ca[ref_pos])
                        moving_atoms.append(ca_atoms[temp_pos])

                if len(fixed_atoms) >= 3:
                    rmsd = superimpose_structures(
                        ref_structure, structure, fixed_atoms, moving_atoms
                    )

                    # Save superimposed structure
                    output_path = self.output_dir / f"{template_id}_superimposed.pdb"
                    write_structure(structure, output_path)
                    results[template_id] = (output_path, rmsd)
                    self.superimposed[template_id] = structure

                    print(f"  {template_id}: RMSD = {rmsd:.2f} Å")

            except Exception as e:
                print(f"  Error superimposing {template_id}: {e}")

        return results

    def visualize_py3dmol(
        self,
        template_ids: List[str],
        width: int = 800,
        height: int = 600
    ):
        """
        Create interactive py3Dmol visualization with match/mismatch coloring.

        Parameters
        ----------
        template_ids : list
            List of template PDB IDs to visualize
        width : int
            Viewer width in pixels
        height : int
            Viewer height in pixels

        Returns
        -------
        py3Dmol.view
            The viewer object for display in Jupyter
        """
        try:
            import py3Dmol
        except ImportError:
            print("py3Dmol not installed. Run: pip install py3Dmol")
            return None

        view = py3Dmol.view(width=width, height=height)

        # Color cycle for different templates
        template_base_colors = [
            '0x1f77b4', '0xff7f0e', '0x2ca02c', '0xd62728',
            '0x9467bd', '0x8c564b', '0xe377c2', '0x7f7f7f'
        ]

        for i, template_id in enumerate(template_ids):
            if template_id not in self.alignment:
                continue

            pdb_4 = template_id[:4].lower()
            pdb_file = self.output_dir / f"{template_id}_superimposed.pdb"

            if not pdb_file.exists():
                pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"

            if not pdb_file.exists():
                continue

            # Read PDB content
            with open(pdb_file) as f:
                pdb_content = f.read()

            # Add model
            view.addModel(pdb_content, 'pdb')

            # Get match/mismatch positions
            template_aligned = self.alignment[template_id]
            matches, mismatches, gaps = get_match_mismatch_positions(
                self.target_aligned, template_aligned
            )

            # Set base style - cartoon with low opacity
            base_color = template_base_colors[i % len(template_base_colors)]
            view.setStyle(
                {'model': i},
                {'cartoon': {'color': base_color, 'opacity': 0.5}}
            )

            # Color matches in green
            if matches:
                view.setStyle(
                    {'model': i, 'resi': matches},
                    {'cartoon': {'color': self.COLORS['match'], 'opacity': 0.9}}
                )

            # Color mismatches in red
            if mismatches:
                view.setStyle(
                    {'model': i, 'resi': mismatches},
                    {'cartoon': {'color': self.COLORS['mismatch'], 'opacity': 0.9}}
                )

            # Color gaps in gray
            if gaps:
                view.setStyle(
                    {'model': i, 'resi': gaps},
                    {'cartoon': {'color': self.COLORS['gap'], 'opacity': 0.5}}
                )

        # Center and zoom
        view.zoomTo()

        return view

    def visualize_with_target(
        self,
        target_pdb: str,
        template_ids: List[str],
        width: int = 800,
        height: int = 600
    ):
        """
        Visualize templates superimposed onto target structure.

        Parameters
        ----------
        target_pdb : str
            Path to target PDB file
        template_ids : list
            Template IDs to show
        width : int
            Viewer width
        height : int
            Viewer height

        Returns
        -------
        py3Dmol.view
            The viewer object
        """
        try:
            import py3Dmol
        except ImportError:
            print("py3Dmol not installed")
            return None

        view = py3Dmol.view(width=width, height=height)

        # Add target structure in blue
        if Path(target_pdb).exists():
            with open(target_pdb) as f:
                view.addModel(f.read(), 'pdb')
            view.setStyle({'model': 0}, {'cartoon': {'color': self.COLORS['target']}})

        # Add templates
        for i, template_id in enumerate(template_ids, start=1):
            pdb_file = self.output_dir / f"{template_id}_superimposed.pdb"
            if not pdb_file.exists():
                continue

            with open(pdb_file) as f:
                view.addModel(f.read(), 'pdb')

            # Get coloring from alignment
            if template_id in self.alignment:
                template_aligned = self.alignment[template_id]
                matches, mismatches, gaps = get_match_mismatch_positions(
                    self.target_aligned, template_aligned
                )

                view.setStyle({'model': i}, {'cartoon': {'opacity': 0.6}})

                if matches:
                    view.setStyle(
                        {'model': i, 'resi': matches},
                        {'cartoon': {'color': self.COLORS['match']}}
                    )
                if mismatches:
                    view.setStyle(
                        {'model': i, 'resi': mismatches},
                        {'cartoon': {'color': self.COLORS['mismatch']}}
                    )

        view.zoomTo()
        return view

    def get_alignment_stats(self, template_id: str) -> Dict[str, int]:
        """
        Get alignment statistics for a template.

        Parameters
        ----------
        template_id : str
            Template PDB ID

        Returns
        -------
        dict
            Statistics with match/mismatch/gap counts
        """
        if template_id not in self.alignment:
            return {}

        template_aligned = self.alignment[template_id]
        matches, mismatches, gaps = get_match_mismatch_positions(
            self.target_aligned, template_aligned
        )

        total = len(matches) + len(mismatches)
        return {
            'matches': len(matches),
            'mismatches': len(mismatches),
            'gaps': len(gaps),
            'identity': (len(matches) / total * 100) if total > 0 else 0,
            'coverage': total
        }

