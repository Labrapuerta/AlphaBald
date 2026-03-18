"""
Template domain analysis module.

Downloads template structures and retrieves domain annotations from PDBe/InterPro.
"""

import os
import subprocess
import gzip
import shutil
import urllib.request
from pathlib import Path
from typing import List, Dict, Optional, Any, Tuple

import requests
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBList


class TemplateProcessor:
    """
    Download and analyze template structures for homology modeling.

    Parameters
    ----------
    homologs : pd.DataFrame
        DataFrame with homolog hits (must have 'PDB_ID' column)
    num_templates : int, optional
        Number of top templates to process. Default 5
    target_path : str or Path, optional
        Path to target FASTA file. Default "target/target.fa"
    verbose : bool, optional
        Print progress messages. Default True

    Attributes
    ----------
    templates : pd.DataFrame
        Top template entries from homologs
    pdb_sequences : dict
        Mapping of PDB_ID to sequence string
    domain_data : list
        List of [pdb_id, sequence, domains] for each template
    target_sequence : str
        Target protein sequence
    """

    def __init__(
        self,
        homologs: pd.DataFrame,
        num_templates: int = 5,
        target_path: str = "target/target.fa",
        verbose: bool = True
    ):
        self.homologs = homologs
        self.num_templates = num_templates
        self.verbose = verbose

        # Paths
        self.target_path = Path(target_path)
        self.templates_dir = Path("Templates")
        self.temp_dir = Path("temp")
        self.alignments_dir = Path("Alignments")
        self.pfam_dir = Path("databases/hmm/Pfam")

        # Ensure directories exist
        self.templates_dir.mkdir(exist_ok=True)
        self.temp_dir.mkdir(exist_ok=True)
        self.alignments_dir.mkdir(exist_ok=True)

        # Load target
        self.target_fasta = SeqIO.read(str(self.target_path), "fasta")
        self.target_sequence = str(self.target_fasta.seq)
        self.target_id = self.target_fasta.id

        # Get top templates
        self.templates = homologs.head(num_templates).copy()

        # Download and process
        self._log("Downloading template structures...")
        self.pdb_files = self._download_templates()

        self._log("Extracting template sequences...")
        self.pdb_sequences = self._extract_sequences()

        self._log("Creating alignment FASTA...")
        self.fasta_for_alignment = self.temp_dir / "templates.fasta"
        self._create_fasta_file()

        self._log("Fetching domain annotations...")
        self.domain_data = self._fetch_all_domains()

        # Search for target domains using HMM
        self._log("Searching target domains with HMM...")
        self.target_domains = self._search_target_domains()

    def _log(self, message: str):
        """Print message if verbose mode is on."""
        if self.verbose:
            print(message)

    def _download_templates(self) -> Dict[str, Path]:
        """Download PDB files for selected templates."""
        pdb_files = {}
        pdbl = PDBList(verbose=False)

        # Get unique 4-letter PDB IDs
        pdb_ids = self.templates["PDB_ID"].str[:4].unique().tolist()

        for pdb_id in pdb_ids:
            pdb_id_lower = pdb_id.lower()
            pdb_file = self.templates_dir / f"pdb{pdb_id_lower}.ent"

            if not pdb_file.exists():
                try:
                    pdbl.retrieve_pdb_file(
                        pdb_id,
                        pdir=str(self.templates_dir),
                        file_format="pdb",
                        overwrite=True
                    )
                except Exception as e:
                    self._log(f"  Warning: Could not download {pdb_id}: {e}")

            if pdb_file.exists():
                pdb_files[pdb_id.upper()] = pdb_file

        return pdb_files

    def _extract_sequences(self) -> Dict[str, str]:
        """Extract sequences from downloaded PDB files."""
        sequences = {}

        for _, row in self.templates.iterrows():
            full_id = row["PDB_ID"]
            pdb_4 = full_id[:4].lower()
            chain = full_id[5] if len(full_id) > 5 else "A"

            pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
            if not pdb_file.exists():
                continue

            seq = self._parse_chain_sequence(pdb_file, chain)
            if seq:
                sequences[full_id] = seq

        return sequences

    @staticmethod
    def _parse_chain_sequence(pdb_file: Path, chain_id: str) -> Optional[str]:
        """Parse sequence for a specific chain from PDB file."""
        try:
            for record in SeqIO.parse(str(pdb_file), "pdb-atom"):
                chain_letter = record.id.split(':')[-1]
                if chain_letter == chain_id:
                    return str(record.seq)
        except Exception:
            pass
        return None

    def _create_fasta_file(self):
        """Create combined FASTA file with target and template sequences."""
        with open(self.fasta_for_alignment, 'w') as f:
            # Write target first
            f.write(f">{self.target_id}\n{self.target_sequence}\n")
            # Write templates
            for full_id, sequence in self.pdb_sequences.items():
                f.write(f">{full_id}\n{sequence}\n")

    def _fetch_all_domains(self) -> List[List[Any]]:
        """Fetch domain annotations for all templates."""
        domain_data = []

        for full_id, sequence in self.pdb_sequences.items():
            pdb_id = full_id[:4].lower()
            chain = full_id[5] if len(full_id) > 5 else "A"

            domains = self._get_pdb_domains(pdb_id, chain)
            domain_data.append([full_id, sequence, domains])

            if self.verbose and domains:
                self._log(f"  {full_id}: {len(domains)} domain(s)")

        return domain_data

    @staticmethod
    def _get_pdb_domains(pdb_id: str, chain: str = "A") -> List[Dict[str, Any]]:
        """
        Fetch InterPro domain mappings for a PDB structure from PDBe API.

        Parameters
        ----------
        pdb_id : str
            4-letter PDB ID
        chain : str
            Chain identifier

        Returns
        -------
        list
            List of domain dictionaries with id, name, start, end
        """
        pdb_id = pdb_id.lower()[:4]
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/interpro/{pdb_id}"

        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return []

            data = response.json()
            domains = []

            if pdb_id not in data:
                return []

            interpro_mappings = data[pdb_id].get('InterPro', {})

            for domain_id, domain_info in interpro_mappings.items():
                domain_name = domain_info.get('name', domain_id)
                mappings = domain_info.get('mappings', [])

                for mapping in mappings:
                    if mapping.get('chain_id') == chain:
                        domains.append({
                            'id': domain_id,
                            'name': domain_name,
                            'start': mapping.get('start', {}).get('residue_number'),
                            'end': mapping.get('end', {}).get('residue_number')
                        })

            return domains

        except Exception:
            return []

    def align_templates(self) -> Path:
        """
        Run ClustalW alignment on target + templates.

        Returns
        -------
        Path
            Path to aligned FASTA file
        """
        aligned_file = self.alignments_dir / "aligned_templates.fasta"

        cmd = [
            "clustalw",
            f"-INFILE={self.fasta_for_alignment}",
            f"-OUTFILE={aligned_file}",
            "-OUTPUT=FASTA",
            "-OUTORDER=INPUT"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self._log(f"Warning: ClustalW alignment may have failed: {result.stderr[:200]}")

        return aligned_file

    def run_hmmscan(self, sequence: str, output_prefix: str = "query") -> Optional[Path]:
        """
        Run hmmscan against Pfam database for a sequence.

        Parameters
        ----------
        sequence : str
            Protein sequence to search
        output_prefix : str
            Prefix for output files

        Returns
        -------
        Path or None
            Path to domain table output, or None if failed
        """
        pfam_hmm = self.pfam_dir / "Pfam-A.hmm"
        if not pfam_hmm.exists():
            self._log("Pfam database not found. Run Setup() first.")
            return None

        # Write sequence to temp file
        query_file = self.temp_dir / f"{output_prefix}.fa"
        with open(query_file, 'w') as f:
            f.write(f">{output_prefix}\n{sequence}\n")

        tblout_file = self.temp_dir / f"{output_prefix}_pfam.tbl"
        domtbl_file = self.temp_dir / f"{output_prefix}_pfam.domtbl"

        cmd = [
            "hmmscan",
            "--tblout", str(tblout_file),
            "--domtblout", str(domtbl_file),
            "-E", "1e-5",
            str(pfam_hmm),
            str(query_file)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self._log(f"hmmscan failed: {result.stderr[:200]}")
            return None

        return domtbl_file

    def _parse_hmmscan_domtbl(self, domtbl_file: Path) -> List[Dict[str, Any]]:
        """
        Parse hmmscan domain table output.

        Parameters
        ----------
        domtbl_file : Path
            Path to domtblout file

        Returns
        -------
        list
            List of domain dictionaries
        """
        if not domtbl_file or not domtbl_file.exists():
            return []

        domains = []
        with open(domtbl_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 23:
                    domains.append({
                        'id': parts[1],  # Pfam accession
                        'name': parts[0],  # Domain name
                        'description': ' '.join(parts[22:]),
                        'e_value': float(parts[6]),
                        'score': float(parts[7]),
                        'start': int(parts[17]),  # ali_from (alignment start in query)
                        'end': int(parts[18]),    # ali_to (alignment end in query)
                    })

        # Sort by start position and remove overlapping duplicates
        domains.sort(key=lambda x: (x['start'], -x['score']))
        return domains

    def _search_target_domains(self) -> List[Dict[str, Any]]:
        """
        Search for domains in the target sequence using Pfam HMM database.

        Returns
        -------
        list
            List of domain dictionaries for target
        """
        domtbl_file = self.run_hmmscan(self.target_sequence, output_prefix="target")
        domains = self._parse_hmmscan_domtbl(domtbl_file)

        if domains and self.verbose:
            self._log(f"  Target: {len(domains)} domain(s) found")
            for d in domains:
                self._log(f"    - {d['name']} ({d['id']}): {d['start']}-{d['end']}")

        return domains

    def search_sequence_domains(self, sequence: str, name: str = "query") -> List[Dict[str, Any]]:
        """
        Search for domains in any sequence using Pfam.

        Parameters
        ----------
        sequence : str
            Protein sequence to analyze
        name : str
            Name for the sequence

        Returns
        -------
        list
            List of domain dictionaries
        """
        domtbl_file = self.run_hmmscan(sequence, output_prefix=name)
        return self._parse_hmmscan_domtbl(domtbl_file)

    def get_target_data(self) -> Tuple[str, str, List[Dict]]:
        """
        Get target information for visualization.

        Returns
        -------
        tuple
            (target_id, target_sequence, target_domains)
        """
        return (self.target_id, self.target_sequence, self.target_domains)

    def get_template_data_list(self) -> List[List[Any]]:
        """
        Get template data formatted for DomainVisualizer.

        Returns
        -------
        list
            List of [pdb_id, sequence, domains] entries
        """
        return self.domain_data

    def __repr__(self) -> str:
        return f"TemplateProcessor(templates={len(self.templates)}, sequences={len(self.pdb_sequences)})"
