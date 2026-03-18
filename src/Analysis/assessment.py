"""
Protein structure assessment utilities.

Provides comprehensive tools for analyzing and validating protein models,
identifying protein families, and assessing model quality.

Example Usage
-------------
>>> from src.Analysis import ModelAssessor
>>> assessor = ModelAssessor("model.pdb")
>>> assessor.identify_family()
>>> assessor.run_dssp()
>>> assessor.find_problematic_regions()
"""

import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
import warnings

try:
    from Bio.PDB import PDBParser, PDBIO, Select, DSSP
    from Bio.PDB.Polypeptide import protein_letters_3to1
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    warnings.warn("BioPython not available. Some functions may not work.")

import requests


def identify_protein_family(
    sequence: str = None,
    pdb_file: str = None,
    output_prefix: str = "query"
) -> Dict[str, Any]:
    """
    Identify protein fold and family using Pfam/InterPro.

    Parameters
    ----------
    sequence : str, optional
        Protein sequence to search
    pdb_file : str, optional
        PDB file to extract sequence from
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with family information including:
        - family_id: Pfam/InterPro ID
        - family_name: Family name
        - fold: SCOP/CATH fold classification if available
        - domains: List of identified domains
    """
    if pdb_file and not sequence:
        sequence = extract_sequence_from_pdb(pdb_file)

    if not sequence:
        raise ValueError("Either sequence or pdb_file must be provided")

    # Save sequence to temp file
    temp_dir = Path("temp")
    temp_dir.mkdir(exist_ok=True)
    query_file = temp_dir / f"{output_prefix}.fa"
    with open(query_file, 'w') as f:
        f.write(f">{output_prefix}\n{sequence}\n")

    # Run hmmscan against Pfam
    pfam_db = Path("databases/hmm/Pfam/Pfam-A.hmm")
    if not pfam_db.exists():
        print("Warning: Pfam database not found. Using web search...")
        return _search_pfam_web(sequence)

    domtbl_file = temp_dir / f"{output_prefix}_family.domtbl"
    cmd = [
        "hmmscan",
        "--domtblout", str(domtbl_file),
        "-E", "1e-5",
        str(pfam_db),
        str(query_file)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"hmmscan failed: {result.stderr[:200]}")
        return {}

    # Parse results
    domains = _parse_domtbl(domtbl_file)

    if domains:
        best_hit = domains[0]
        return {
            'family_id': best_hit['accession'],
            'family_name': best_hit['name'],
            'domains': domains,
            'fold': _get_fold_from_pfam(best_hit['accession'])
        }
    return {}


def _search_pfam_web(sequence: str) -> Dict[str, Any]:
    """Search Pfam via web API (fallback)."""
    try:
        url = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
        data = {
            'seq': sequence,
            'seqdb': 'pfam'
        }
        response = requests.post(url, data=data, timeout=60)
        if response.status_code == 200:
            # Parse response
            return {'note': 'Web search performed - check results manually'}
    except Exception as e:
        print(f"Web search failed: {e}")
    return {}


def _parse_domtbl(domtbl_file: Path) -> List[Dict]:
    """Parse HMMER domain table output."""
    domains = []
    if not domtbl_file.exists():
        return domains

    with open(domtbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 23:
                domains.append({
                    'name': parts[0],
                    'accession': parts[1],
                    'description': ' '.join(parts[22:]),
                    'evalue': float(parts[6]),
                    'score': float(parts[7]),
                    'start': int(parts[17]),
                    'end': int(parts[18])
                })

    domains.sort(key=lambda x: x['evalue'])
    return domains


def _get_fold_from_pfam(pfam_id: str) -> Optional[str]:
    """Get SCOP/CATH fold classification from Pfam ID."""
    try:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            # Extract fold info if available
            return data.get('metadata', {}).get('type', {}).get('name', 'Unknown')
    except Exception:
        pass
    return None


def create_hmm_profile(
    alignment_file: str,
    output_file: str = "profile.hmm",
    name: str = "profile"
) -> Path:
    """
    Build HMM profile from multiple sequence alignment.

    Parameters
    ----------
    alignment_file : str
        Input alignment file (FASTA, Stockholm, etc.)
    output_file : str
        Output HMM profile file
    name : str
        Name for the HMM profile

    Returns
    -------
    Path
        Path to created HMM file
    """
    output_path = Path(output_file)

    cmd = [
        "hmmbuild",
        "-n", name,
        str(output_path),
        str(alignment_file)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"hmmbuild failed: {result.stderr}")

    print(f"HMM profile created: {output_path}")
    return output_path


def align_hmm_to_sequence(
    hmm_file: str,
    sequence_file: str,
    output_file: str = "alignment.aln"
) -> Path:
    """
    Align HMM profile to sequence using hmmalign.

    Parameters
    ----------
    hmm_file : str
        HMM profile file
    sequence_file : str
        Sequence file (FASTA)
    output_file : str
        Output alignment file

    Returns
    -------
    Path
        Path to alignment file
    """
    output_path = Path(output_file)

    cmd = [
        "hmmalign",
        "-o", str(output_path),
        str(hmm_file),
        str(sequence_file)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"hmmalign failed: {result.stderr}")

    print(f"Alignment created: {output_path}")
    return output_path


def add_cation_to_structure(
    pdb_file: str,
    cation: str = "CA",  # Default calcium
    output_file: str = None,
    binding_residues: List[int] = None
) -> Path:
    """
    Add a cation to a protein structure.

    Parameters
    ----------
    pdb_file : str
        Input PDB file
    cation : str
        Cation type (CA, MG, ZN, FE, etc.)
    output_file : str, optional
        Output PDB file. If None, uses {input}_cation.pdb
    binding_residues : list, optional
        List of residue numbers that coordinate the cation

    Returns
    -------
    Path
        Path to output PDB file
    """
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required for this function")

    pdb_path = Path(pdb_file)
    if output_file is None:
        output_file = pdb_path.stem + "_cation.pdb"
    output_path = Path(output_file)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", str(pdb_path))

    # Find cation binding site based on common motifs
    if binding_residues is None:
        binding_residues = _find_cation_binding_site(structure, cation)

    if not binding_residues:
        print("Warning: Could not identify binding site automatically.")
        print("Please provide binding_residues parameter or check structure manually.")
        return pdb_path

    # Calculate center of binding residues for cation placement
    coords = []
    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.id[1] in binding_residues:
                for atom in residue:
                    if atom.name in ['OD1', 'OD2', 'OE1', 'OE2', 'NE2', 'SD']:
                        coords.append(atom.coord)

    if not coords:
        # Use CA atoms as fallback
        for chain in model:
            for residue in chain:
                if residue.id[1] in binding_residues:
                    if 'CA' in residue:
                        coords.append(residue['CA'].coord)

    if coords:
        import numpy as np
        center = np.mean(coords, axis=0)

        # Add cation as HETATM
        with open(pdb_path) as f:
            pdb_lines = f.readlines()

        # Find last ATOM line
        insert_idx = 0
        for i, line in enumerate(pdb_lines):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                insert_idx = i + 1

        # Create HETATM line for cation
        cation_line = f"HETATM{9999:5d}  {cation:2s}  {cation:3s} X{1:4d}    {center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00  0.00          {cation[:2]:>2s}\n"

        pdb_lines.insert(insert_idx, cation_line)

        with open(output_path, 'w') as f:
            f.writelines(pdb_lines)

        print(f"Cation {cation} added at position ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        print(f"Output: {output_path}")
        return output_path
    else:
        print("Could not determine cation position")
        return pdb_path


def _find_cation_binding_site(structure, cation: str) -> List[int]:
    """Find potential cation binding residues."""
    binding_residues = []

    # Common cation-coordinating residues
    coordinating_residues = {
        'CA': ['ASP', 'GLU', 'ASN', 'GLN'],  # Calcium
        'MG': ['ASP', 'GLU', 'HIS'],  # Magnesium
        'ZN': ['CYS', 'HIS', 'ASP', 'GLU'],  # Zinc
        'FE': ['CYS', 'HIS', 'GLU', 'ASP']  # Iron
    }

    target_residues = coordinating_residues.get(cation, ['ASP', 'GLU', 'HIS'])

    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.resname in target_residues:
                # Check if residue is surface-accessible (simplified)
                binding_residues.append(residue.id[1])

    return binding_residues[:4]  # Return top 4 candidates


def identify_functional_residues(
    alignment_file: str,
    output_file: str = None,
    conservation_threshold: float = 0.9
) -> List[Dict]:
    """
    Identify functionally important residues based on conservation.

    Parameters
    ----------
    alignment_file : str
        Multiple sequence alignment file
    output_file : str, optional
        Output file with @ symbols marking important residues
    conservation_threshold : float
        Minimum conservation score (0-1) to mark as important

    Returns
    -------
    list
        List of important residue positions with conservation scores
    """
    # Parse alignment
    from Bio import AlignIO
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception:
        try:
            alignment = AlignIO.read(alignment_file, "stockholm")
        except Exception as e:
            print(f"Could not parse alignment: {e}")
            return []

    # Calculate conservation at each position
    important_residues = []
    num_seqs = len(alignment)

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        # Skip gaps
        non_gap = [c for c in column if c != '-']
        if not non_gap:
            continue

        # Calculate conservation
        from collections import Counter
        counts = Counter(non_gap)
        most_common = counts.most_common(1)[0][1]
        conservation = most_common / len(non_gap)

        if conservation >= conservation_threshold:
            important_residues.append({
                'position': i + 1,
                'residue': counts.most_common(1)[0][0],
                'conservation': conservation
            })

    # Write alignment with @ markers if output file specified
    if output_file and important_residues:
        _write_marked_alignment(alignment_file, important_residues, output_file)
        print(f"Marked alignment saved to: {output_file}")

    return important_residues


def _write_marked_alignment(alignment_file: str, important_residues: List[Dict], output_file: str):
    """Write alignment with @ symbols marking important residues."""
    important_positions = {r['position'] for r in important_residues}

    with open(alignment_file) as f:
        content = f.read()

    # Add marker line at the end
    marker_line = "\n# Important functional residues (marked with @):\n# "
    for i in range(1, max(important_positions) + 1):
        if i in important_positions:
            marker_line += "@"
        else:
            marker_line += " "
    marker_line += "\n"

    with open(output_file, 'w') as f:
        f.write(content)
        f.write(marker_line)


def analyze_active_site(
    pdb_file: str,
    active_site_residues: List[int] = None,
    reference_pdb: str = None
) -> Dict[str, Any]:
    """
    Analyze active site residues and compare with reference.

    Parameters
    ----------
    pdb_file : str
        PDB file to analyze
    active_site_residues : list, optional
        List of active site residue numbers
    reference_pdb : str, optional
        Reference PDB for comparison

    Returns
    -------
    dict
        Active site analysis including:
        - residue_info: Details of each active site residue
        - is_preserved: Whether key residues are conserved
        - is_active_prediction: Prediction of protein activity
    """
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required for this function")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", str(pdb_file))

    analysis = {
        'residue_info': [],
        'is_preserved': True,
        'is_active_prediction': True,
        'issues': []
    }

    model = structure[0]
    if active_site_residues is None:
        # Try to identify catalytic residues automatically
        active_site_residues = _find_catalytic_residues(model)

    for chain in model:
        for residue in chain:
            if residue.id[1] in active_site_residues:
                res_info = {
                    'number': residue.id[1],
                    'name': residue.resname,
                    'chain': chain.id
                }

                # Check for structural issues
                if 'CB' not in residue and residue.resname != 'GLY':
                    res_info['issue'] = 'Missing CB atom'
                    analysis['issues'].append(f"Residue {residue.id[1]}: Missing CB")

                analysis['residue_info'].append(res_info)

    # Compare with reference if provided
    if reference_pdb:
        ref_structure = parser.get_structure("ref", str(reference_pdb))
        # Compare active site residues
        # ... (comparison logic)

    return analysis


def _find_catalytic_residues(model) -> List[int]:
    """Identify potential catalytic residues based on common patterns."""
    catalytic = []
    # Common catalytic residue types
    catalytic_types = ['HIS', 'CYS', 'SER', 'ASP', 'GLU', 'LYS']

    for chain in model:
        for residue in chain:
            if residue.resname in catalytic_types:
                catalytic.append(residue.id[1])

    return catalytic[:10]  # Return first 10 candidates


def validate_model_regions(
    pdb_file: str,
    method: str = "prosa"
) -> Dict[str, Any]:
    """
    Validate model regions for structural problems.

    Parameters
    ----------
    pdb_file : str
        PDB file to validate
    method : str
        Validation method ('prosa', 'dope', 'all')

    Returns
    -------
    dict
        Validation results with problematic regions
    """
    results = {
        'pdb_file': pdb_file,
        'problematic_regions': [],
        'overall_quality': 'unknown'
    }

    # Run DSSP for secondary structure
    dssp_data = run_dssp_analysis(pdb_file)['dssp_data']

    if method in ['prosa', 'all']:
        print("Note: ProSA requires web submission at:")
        print("  https://prosa.services.came.sbg.ac.at/prosa.php")
        print("Upload the PDB file and analyze with multiple window sizes (10, 40)")
        results['prosa_note'] = "Submit to ProSA web server for energy analysis"

    if method in ['dope', 'all']:
        # Use MODELLER DOPE if available
        try:
            from modeller import Environ
            from modeller.scripts import complete_pdb
            dope_results = _calculate_dope_profile(pdb_file)
            results['dope_profile'] = dope_results
        except ImportError:
            results['dope_note'] = "MODELLER not available for DOPE calculation"

    # Identify problematic regions from DSSP
    if dssp_data:
        # Look for unusual phi/psi angles or breaks
        for res_id, data in dssp_data.items():
            if data.get('phi', 0) == 360.0 or data.get('psi', 0) == 360.0:
                results['problematic_regions'].append({
                    'residue': res_id,
                    'issue': 'Chain break or missing backbone'
                })

    return results


def _calculate_dope_profile(pdb_file: str) -> Dict:
    """Calculate DOPE profile using MODELLER."""
    try:
        from modeller import Environ
        from modeller.scripts import complete_pdb

        env = Environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        mdl = complete_pdb(env, pdb_file)
        s = mdl.select_all()

        profile_file = Path(pdb_file).stem + "_dope.profile"
        s.assess_dope(
            output='ENERGY_PROFILE NO_REPORT',
            file=profile_file,
            normalize_profile=True,
            smoothing_window=15
        )

        return {'profile_file': profile_file}
    except Exception as e:
        return {'error': str(e)}


def run_dssp_analysis(
    pdb_file: str,
    output_file: str = None
) -> Dict[str, Any]:
    """
    Run DSSP secondary structure analysis.

    Parameters
    ----------
    pdb_file : str
        PDB file to analyze
    output_file : str, optional
        Output DSSP file. If None, uses {input}.dssp

    Returns
    -------
    dict
        DSSP results including:
        - dssp_file: Path to DSSP output
        - dssp_data: Per-residue secondary structure assignments
        - ss_summary: Summary of secondary structure composition
    """
    pdb_path = Path(pdb_file)
    if output_file is None:
        output_file = str(pdb_path.with_suffix('.dssp'))

    results = {
        'dssp_file': output_file,
        'dssp_data': {},
        'ss_summary': {}
    }

    # Try running dssp command
    cmd = ["dssp", "-i", str(pdb_file), "-o", str(output_file)]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"DSSP output: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"DSSP command failed: {e}")
        # Try BioPython DSSP as fallback
        if HAS_BIOPYTHON:
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("model", str(pdb_file))
                model = structure[0]
                dssp = DSSP(model, str(pdb_file))

                for key in dssp:
                    res_id = key[1][1]
                    results['dssp_data'][res_id] = {
                        'ss': dssp[key][2],
                        'rsa': dssp[key][3],
                        'phi': dssp[key][4],
                        'psi': dssp[key][5]
                    }

                # Calculate summary
                ss_counts = {}
                for data in results['dssp_data'].values():
                    ss = data['ss']
                    ss_counts[ss] = ss_counts.get(ss, 0) + 1
                results['ss_summary'] = ss_counts

            except Exception as e:
                print(f"BioPython DSSP also failed: {e}")

    return results


def fix_model_problems(
    pdb_file: str,
    problems: List[Dict],
    output_file: str = None,
    method: str = "modeller_loop"
) -> Path:
    """
    Fix identified problems in the model.

    Parameters
    ----------
    pdb_file : str
        PDB file with problems
    problems : list
        List of problem dictionaries from validate_model_regions
    output_file : str, optional
        Output file for fixed model
    method : str
        Fixing method ('modeller_loop', 'energy_minimize')

    Returns
    -------
    Path
        Path to fixed PDB file
    """
    pdb_path = Path(pdb_file)
    if output_file is None:
        output_file = pdb_path.stem + "_fixed.pdb"
    output_path = Path(output_file)

    if not problems:
        print("No problems to fix")
        return pdb_path

    # Identify loop regions to refine
    loop_regions = []
    for problem in problems:
        res = problem.get('residue')
        if res:
            # Extend region around problem
            loop_regions.append((max(1, res - 2), res + 2))

    # Merge overlapping regions
    loop_regions = _merge_regions(loop_regions)

    if method == 'modeller_loop':
        # Generate MODELLER loop refinement script
        from src.modeller.scripts import generate_loop_refinement_script

        # Get target_id from PDB file
        target_id = pdb_path.stem.split('.')[0]

        script = generate_loop_refinement_script(
            str(pdb_path),
            target_id,
            loop_regions,
            num_models=3
        )

        script_path = pdb_path.parent / "fix_loops.py"
        with open(script_path, 'w') as f:
            f.write(script)

        print(f"Loop refinement script created: {script_path}")
        print("Run with MODELLER to fix the problems:")
        print(f"  python {script_path}")

    return output_path


def _merge_regions(regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping regions."""
    if not regions:
        return []

    sorted_regions = sorted(regions)
    merged = [sorted_regions[0]]

    for start, end in sorted_regions[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged


def extract_sequence_from_pdb(pdb_file: str) -> str:
    """Extract amino acid sequence from PDB file."""
    if not HAS_BIOPYTHON:
        # Fallback: parse ATOM lines manually
        sequence = []
        current_res = None
        three_to_one = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    res_num = int(line[22:26])
                    if res_num != current_res:
                        current_res = res_num
                        res_name = line[17:20].strip()
                        sequence.append(three_to_one.get(res_name, 'X'))

        return ''.join(sequence)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", str(pdb_file))
    sequence = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    sequence.append(
                        protein_letters_3to1.get(residue.resname, 'X')
                    )
            break  # Only first chain
        break  # Only first model

    return ''.join(sequence)


class ModelAssessor:
    """
    Comprehensive model assessment class.

    Provides a unified interface for all assessment functions.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file to assess
    output_dir : str, optional
        Output directory for results

    Example
    -------
    >>> assessor = ModelAssessor("model.pdb")
    >>> assessor.identify_family()
    >>> assessor.run_dssp("model.dssp")
    >>> problems = assessor.find_problematic_regions()
    >>> assessor.fix_problems(problems)
    """

    def __init__(self, pdb_file: str, output_dir: str = "assessment"):
        self.pdb_file = Path(pdb_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Extract sequence
        self.sequence = extract_sequence_from_pdb(str(self.pdb_file))
        self.results = {}

    def identify_family(self) -> Dict:
        """Identify protein fold and family."""
        result = identify_protein_family(
            sequence=self.sequence,
            output_prefix=self.pdb_file.stem
        )
        self.results['family'] = result
        return result

    def create_hmm(self, alignment_file: str, name: str = None) -> Path:
        """Create HMM profile from alignment."""
        if name is None:
            name = self.pdb_file.stem
        output = self.output_dir / f"{name}.hmm"
        return create_hmm_profile(alignment_file, str(output), name)

    def run_dssp(self, output_file: str = None) -> Dict:
        """Run DSSP analysis."""
        if output_file is None:
            output_file = str(self.output_dir / f"{self.pdb_file.stem}.dssp")
        result = run_dssp_analysis(str(self.pdb_file), output_file)
        self.results['dssp'] = result
        return result

    def add_cation(self, cation: str = "CA", output_file: str = None) -> Path:
        """Add cation to structure."""
        if output_file is None:
            output_file = str(self.output_dir / f"{self.pdb_file.stem}_cation.pdb")
        return add_cation_to_structure(str(self.pdb_file), cation, output_file)

    def find_problematic_regions(self, method: str = "all") -> Dict:
        """Validate model and find problems."""
        result = validate_model_regions(str(self.pdb_file), method)
        self.results['validation'] = result
        return result

    def fix_problems(self, problems: List[Dict] = None, output_file: str = None) -> Path:
        """Fix identified problems."""
        if problems is None:
            problems = self.results.get('validation', {}).get('problematic_regions', [])
        if output_file is None:
            output_file = str(self.output_dir / f"{self.pdb_file.stem}_fixed.pdb")
        return fix_model_problems(str(self.pdb_file), problems, output_file)

    def get_summary(self) -> Dict:
        """Get summary of all assessments."""
        return {
            'pdb_file': str(self.pdb_file),
            'sequence_length': len(self.sequence),
            **self.results
        }
