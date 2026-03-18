"""
Preprocessing module for loading and validating target sequences.
"""

import os
from pathlib import Path
from typing import Optional, Dict, Any
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.PDB import Selection


class TargetPreprocessor:
    """
    Load and preprocess target sequence from FASTA or PDB file.

    Parameters
    ----------
    target_dir : str or Path, optional
        Directory containing target file. Default "target".
    fasta_file : str, optional
        Specific FASTA file name. If None, auto-detects from target_dir.
    pdb_file : str, optional
        Specific PDB file name. If None and fasta_file is None, auto-detects.

    Attributes
    ----------
    sequence : str
        Target protein sequence
    name : str
        Target sequence identifier
    file_path : Path
        Path to the source file
    source_type : str
        Either "fasta" or "pdb"
    resolution : float or None
        Resolution if loaded from PDB, otherwise None
    """

    def __init__(
        self,
        target_dir: str = "target",
        fasta_file: Optional[str] = None,
        pdb_file: Optional[str] = None
    ):
        self.target_dir = Path(target_dir)
        self.sequence: str = ""
        self.name: str = ""
        self.file_path: Optional[Path] = None
        self.source_type: str = ""
        self.resolution: Optional[float] = None

        if pdb_file:
            self._load_from_pdb(self.target_dir / pdb_file)
        elif fasta_file:
            self._load_from_fasta(self.target_dir / fasta_file)
        else:
            self._auto_detect_and_load()

    def _auto_detect_and_load(self):
        """Auto-detect target file in the target directory."""
        if not self.target_dir.exists():
            raise FileNotFoundError(f"Target directory not found: {self.target_dir}")

        files = list(self.target_dir.iterdir())

        # Prefer FASTA files
        for f in files:
            if f.suffix.lower() in ['.fa', '.fasta', '.faa']:
                self._load_from_fasta(f)
                return

        # Fall back to PDB
        for f in files:
            if f.suffix.lower() in ['.pdb', '.ent']:
                self._load_from_pdb(f)
                return

        raise FileNotFoundError(
            f"No FASTA or PDB file found in {self.target_dir}. "
            "Please add a .fa, .fasta, or .pdb file."
        )

    def _load_from_fasta(self, file_path: Path):
        """Load sequence from FASTA file."""
        if not file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

        record = SeqIO.read(str(file_path), "fasta")
        self.sequence = str(record.seq)
        self.name = record.id
        self.file_path = file_path
        self.source_type = "fasta"

    def _load_from_pdb(self, file_path: Path):
        """Load sequence from PDB file."""
        if not file_path.exists():
            raise FileNotFoundError(f"PDB file not found: {file_path}")

        parser = PDBParser(QUIET=True, PERMISSIVE=True)
        structure = parser.get_structure("target", str(file_path))

        # Get resolution if available
        self.resolution = structure.header.get("resolution")
        self.name = structure.header.get("name", file_path.stem)

        # Extract sequence from ATOM records
        residues = [
            res for res in Selection.unfold_entities(structure, "R")
            if res.get_id()[0] == " "  # Standard residues only
        ]
        self.sequence = "".join([
            protein_letters_3to1.get(res.get_resname(), "X")
            for res in residues
        ])

        self.file_path = file_path
        self.source_type = "pdb"

    def to_dict(self) -> Dict[str, Any]:
        """Return target information as dictionary."""
        return {
            "name": self.name,
            "sequence": self.sequence,
            "length": len(self.sequence),
            "source_type": self.source_type,
            "file_path": str(self.file_path),
            "resolution": self.resolution
        }

    def write_fasta(self, output_path: Optional[Path] = None) -> Path:
        """
        Write the target sequence to a FASTA file.

        Parameters
        ----------
        output_path : Path, optional
            Output file path. If None, writes to target/target.fa

        Returns
        -------
        Path
            Path to the written file
        """
        if output_path is None:
            output_path = self.target_dir / "target.fa"

        with open(output_path, "w") as f:
            f.write(f">{self.name}\n{self.sequence}\n")

        return output_path

    def __repr__(self) -> str:
        return f"TargetPreprocessor(name='{self.name}', length={len(self.sequence)}, source='{self.source_type}')"
