"""
Setup module for initializing directories and databases for the homology modeling pipeline.
"""

import os
import shutil
import gzip
import urllib.request
import subprocess
from pathlib import Path


class Setup:
    """
    Setup of all directories and databases for the homology modeling pipeline.

    This class creates necessary directories and downloads/prepares databases:
    - SwissProt BLAST database for homology searches
    - PDB sequence database (pdbaa) for structural template searches
    - Pfam HMM database for domain analysis

    Parameters
    ----------
    skip_databases : bool, optional
        If True, skip database downloads (useful if already set up). Default False.
    clear_temp : bool, optional
        If True, clear temporary directories on init. Default True.

    Attributes
    ----------
    base_dir : Path
        Base project directory
    databases_dir : Path
        Path to databases folder
    """

    def __init__(self, skip_databases: bool = False, clear_temp: bool = True):
        self.base_dir = Path.cwd()
        self.databases_dir = self.base_dir / "databases"
        self.swissprot_dir = self.databases_dir / "swissprot"
        self.pdb_seq_dir = self.databases_dir / "pdb_seq"
        self.pfam_dir = self.databases_dir / "hmm" / "Pfam"

        self.create_directories(clear_temp=clear_temp)

        if not skip_databases:
            print("Setting up databases (this may take a while on first run)...")
            self.setup_swissprot_database()
            self.setup_pdb_database()
            self.setup_pfam_database()
            print("Database setup complete.")

    def create_directories(self, clear_temp: bool = True):
        """
        Create necessary directories for the pipeline.

        Parameters
        ----------
        clear_temp : bool
            If True, clear and recreate temporary directories.
        """
        # Directories to clear on each run
        temp_dirs = ["Alignments", "Templates", "temp"]

        # Directories to create if they don't exist
        persistent_dirs = [
            "databases",
            "databases/hmm/Pfam",
            "databases/pdb_seq",
            "databases/swissprot",
            "target"
        ]

        if clear_temp:
            for directory in temp_dirs:
                dir_path = self.base_dir / directory
                if dir_path.exists():
                    shutil.rmtree(dir_path)
                dir_path.mkdir(exist_ok=True)
        else:
            for directory in temp_dirs:
                (self.base_dir / directory).mkdir(exist_ok=True)

        for directory in persistent_dirs:
            (self.base_dir / directory).mkdir(parents=True, exist_ok=True)

    def setup_swissprot_database(self) -> bool:
        """
        Download and setup SwissProt BLAST database.

        Returns
        -------
        bool
            True if database already exists, False if newly downloaded.
        """
        pin_file = self.swissprot_dir / "swissprot.pin"

        if pin_file.exists():
            print("  SwissProt database already exists.")
            return True

        print("  Downloading SwissProt database...")
        result = subprocess.run(
            ["update_blastdb.pl", "--decompress", "swissprot"],
            cwd=str(self.swissprot_dir),
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  Warning: SwissProt download may have failed: {result.stderr}")

        return False

    def setup_pdb_database(self) -> bool:
        """
        Download and setup PDB sequence BLAST database (pdbaa).

        Returns
        -------
        bool
            True if database already exists, False if newly downloaded.
        """
        pin_file = self.pdb_seq_dir / "pdbaa.pin"

        if pin_file.exists():
            print("  PDB sequence database already exists.")
            return True

        print("  Downloading PDB sequence database...")
        result = subprocess.run(
            ["update_blastdb.pl", "--decompress", "pdbaa"],
            cwd=str(self.pdb_seq_dir),
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  Warning: PDB database download may have failed: {result.stderr}")

        return False

    def setup_pfam_database(self) -> str:
        """
        Download and setup Pfam HMM database.

        Returns
        -------
        str
            Path to the Pfam HMM file.
        """
        hmm_file = self.pfam_dir / "Pfam-A.hmm"
        gz_file = self.pfam_dir / "Pfam-A.hmm.gz"
        pressed_check = self.pfam_dir / "Pfam-A.hmm.h3m"

        if pressed_check.exists():
            print("  Pfam database already exists and is indexed.")
            return str(hmm_file)

        self.pfam_dir.mkdir(parents=True, exist_ok=True)

        if not hmm_file.exists():
            if not gz_file.exists():
                print("  Downloading Pfam database...")
                url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
                urllib.request.urlretrieve(url, str(gz_file))

            print("  Extracting Pfam database...")
            with gzip.open(gz_file, 'rb') as f_in:
                with open(hmm_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            gz_file.unlink()

        print("  Indexing Pfam database with hmmpress...")
        result = subprocess.run(
            ["hmmpress", str(hmm_file)],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"hmmpress failed: {result.stderr}")

        return str(hmm_file)

    def create_pdb_fasta(self) -> Path:
        """
        Create FASTA file from PDB BLAST database for HMMER searches.

        Returns
        -------
        Path
            Path to the generated FASTA file.
        """
        fasta_file = self.pdb_seq_dir / "pdbaa.fasta"

        if fasta_file.exists():
            return fasta_file

        print("  Converting PDB database to FASTA...")
        result = subprocess.run(
            ["blastdbcmd", "-db", str(self.pdb_seq_dir / "pdbaa"),
             "-entry", "all", "-out", str(fasta_file)],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"blastdbcmd failed: {result.stderr}")

        return fasta_file
