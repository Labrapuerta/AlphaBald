"""
Template retrieval module for homology modeling pipeline.

This module handles searching for homologous structures using PSI-BLAST and HMMER.
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, List, Tuple

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


class TemplateRetriever:
    """
    Search for homologous template structures using PSI-BLAST and HMMER.

    The search strategy:
    1. Build PSSM from PSI-BLAST against SwissProt
    2. Search PDB database using the PSSM
    3. Create MSA from top SwissProt hits
    4. Search PDB using the MSA profile
    5. Optionally run jackhmmer for iterative HMM search
    6. Rank and combine all hits

    Parameters
    ----------
    target_path : str or Path
        Path to target FASTA file
    e_value_threshold : float, optional
        E-value cutoff for filtering hits. Default 1e-5
    num_iterations : int, optional
        Number of PSI-BLAST iterations. Default 5
    num_hits : int, optional
        Number of top hits to use for MSA. Default 10
    run_hmmer : bool, optional
        Whether to run jackhmmer search. Default True
    verbose : bool, optional
        Print progress messages. Default True

    Attributes
    ----------
    homologs : pd.DataFrame
        Combined and ranked template hits
    top_swissprot_hits : pd.DataFrame
        Top hits from SwissProt search
    top_psi_hits : pd.DataFrame
        Hits from PSSM-based PDB search
    top_msa_hits : pd.DataFrame
        Hits from MSA-based PDB search
    """

    # BLAST output columns
    PDB_BLAST_COLUMNS = [
        "PDB_ID", "Bitscore", "E-value", "Identity", "Query_Coverage",
        "Query_Start", "Query_End", "Subject_Start", "Subject_End",
        "Query_Length", "Subject_Title"
    ]

    def __init__(
        self,
        target_path: str = "target/target.fa",
        e_value_threshold: float = 1e-5,
        num_iterations: int = 5,
        num_hits: int = 10,
        run_hmmer: bool = True,
        verbose: bool = True
    ):
        self.target_path = Path(target_path)
        self.e_value_threshold = e_value_threshold
        self.num_iterations = num_iterations
        self.num_hits = num_hits
        self.verbose = verbose

        # Paths
        self.temp_dir = Path("temp")
        self.alignments_dir = Path("Alignments")
        self.pssm_file = self.temp_dir / "target.pssm"
        self.swissprot_db = Path("databases/swissprot/swissprot")
        self.pdb_db = Path("databases/pdb_seq/pdbaa")
        self.pdb_fasta = Path("databases/pdb_seq/pdbaa.fasta")
        self.unaligned_fasta = self.temp_dir / "unaligned_homologs.fasta"
        self.aligned_fasta = self.alignments_dir / "aligned_homologs.fa"

        # Ensure directories exist
        self.temp_dir.mkdir(exist_ok=True)
        self.alignments_dir.mkdir(exist_ok=True)

        # Run pipeline
        self._log("Starting template search pipeline...")
        self._generate_pssm()
        self._search_pdb_with_pssm()
        self._create_msa_and_search()

        if run_hmmer:
            self._run_jackhmmer()
            self._build_hmm_profile()

        self._rank_homologs()
        self._log(f"Found {len(self.homologs)} unique templates.")

    def _log(self, message: str):
        """Print message if verbose mode is on."""
        if self.verbose:
            print(message)

    def _run_command(self, cmd: List[str], description: str = "") -> subprocess.CompletedProcess:
        """Run a shell command and handle errors."""
        self._log(f"  {description}..." if description else f"  Running: {' '.join(cmd[:3])}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 and self.verbose:
            print(f"  Warning: {description} may have issues: {result.stderr[:200]}")
        return result

    def _generate_pssm(self):
        """Generate Position-Specific Scoring Matrix using PSI-BLAST against SwissProt."""
        self._log("Step 1: Building PSSM from SwissProt search...")

        cmd = [
            "psiblast",
            "-query", str(self.target_path),
            "-db", str(self.swissprot_db),
            "-num_iterations", str(self.num_iterations),
            "-out_pssm", str(self.pssm_file),
            "-outfmt", "6 sacc bitscore evalue",
            "-out", str(self.temp_dir / "target_swissprot.out")
        ]
        self._run_command(cmd, "Running PSI-BLAST against SwissProt")

        # Parse results
        output_file = self.temp_dir / "target_swissprot.out"
        if output_file.exists() and output_file.stat().st_size > 0:
            self.top_swissprot_hits = pd.read_csv(
                output_file, sep="\t", header=None,
                names=["Accession", "Bitscore", "E-value"]
            )
            self.top_swissprot_hits = self.top_swissprot_hits[
                self.top_swissprot_hits["E-value"] < self.e_value_threshold
            ].drop_duplicates(subset=["Accession"])
            self.top_swissprot_hits.sort_values(
                by=["Bitscore", "E-value"],
                ascending=[False, True],
                inplace=True
            )
        else:
            self.top_swissprot_hits = pd.DataFrame(columns=["Accession", "Bitscore", "E-value"])

    def _search_pdb_with_pssm(self):
        """Search PDB database using the generated PSSM."""
        self._log("Step 2: Searching PDB with PSSM...")

        cmd = [
            "psiblast",
            "-in_pssm", str(self.pssm_file),
            "-db", str(self.pdb_db),
            "-out", str(self.temp_dir / "target_pdbaa.out"),
            "-outfmt", "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
        ]
        self._run_command(cmd, "Running PSI-BLAST against PDB")

        self.top_psi_hits = self._parse_pdb_blast_output(
            self.temp_dir / "target_pdbaa.out",
            origin="PSSM-based Search"
        )

    def _create_msa_and_search(self):
        """Create MSA from top SwissProt hits and search PDB."""
        self._log("Step 3: Creating MSA and searching PDB...")

        # Prepare unaligned FASTA with target + top SwissProt hits
        self._prepare_unaligned_fasta()

        # Run ClustalW alignment
        if self.unaligned_fasta.exists():
            cmd = [
                "clustalw",
                f"-INFILE={self.unaligned_fasta}",
                f"-OUTFILE={self.aligned_fasta}",
                "-OUTPUT=FASTA",
                "-OUTORDER=INPUT"
            ]
            self._run_command(cmd, "Running ClustalW alignment")

        # Search PDB with MSA
        if self.aligned_fasta.exists():
            cmd = [
                "psiblast",
                "-in_msa", str(self.aligned_fasta),
                "-db", str(self.pdb_db),
                "-out", str(self.temp_dir / "msa_top_hits.out"),
                "-outfmt", "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
            ]
            self._run_command(cmd, "Running MSA-based PSI-BLAST")

            self.top_msa_hits = self._parse_pdb_blast_output(
                self.temp_dir / "msa_top_hits.out",
                origin="MSA-based Search"
            )
        else:
            self.top_msa_hits = pd.DataFrame(columns=self.PDB_BLAST_COLUMNS + ["Origin"])

    def _prepare_unaligned_fasta(self):
        """Prepare FASTA file with target and top SwissProt sequences."""
        # Remove old file
        if self.unaligned_fasta.exists():
            self.unaligned_fasta.unlink()

        # Write target sequence first
        with open(self.unaligned_fasta, "w") as fout:
            with open(self.target_path) as fin:
                content = fin.read()
                fout.write(content)
                if not content.endswith("\n"):
                    fout.write("\n")

        # Get top SwissProt sequences
        if len(self.top_swissprot_hits) == 0:
            return

        accessions = self.top_swissprot_hits.head(self.num_hits)["Accession"].tolist()
        if not accessions:
            return

        cmd = [
            "blastdbcmd",
            "-db", str(self.swissprot_db),
            "-entry", ",".join(accessions)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        with open(self.unaligned_fasta, "a") as fout:
            fout.write(result.stdout)

    def _run_jackhmmer(self):
        """Run jackhmmer iterative HMM search against PDB sequences."""
        self._log("Step 4: Running jackhmmer search...")

        # Create FASTA from BLAST db if needed
        if not self.pdb_fasta.exists():
            cmd = [
                "blastdbcmd",
                "-db", str(self.pdb_db),
                "-entry", "all",
                "-out", str(self.pdb_fasta)
            ]
            self._run_command(cmd, "Converting PDB database to FASTA")

        cmd = [
            "jackhmmer",
            "-N", str(self.num_iterations),
            "--tblout", str(self.temp_dir / "jackhmmer_hits.out"),
            str(self.target_path),
            str(self.pdb_fasta)
        ]
        self._run_command(cmd, "Running jackhmmer")

    def _build_hmm_profile(self):
        """Build HMM profile from the aligned sequences."""
        self._log("Step 5: Building HMM profile...")

        if not self.aligned_fasta.exists():
            return

        self.hmm_profile = self.temp_dir / "target.hmm"
        cmd = [
            "hmmbuild",
            str(self.hmm_profile),
            str(self.aligned_fasta)
        ]
        self._run_command(cmd, "Running hmmbuild")

    def search_pfam_domains(self, sequence_file: Optional[str] = None) -> pd.DataFrame:
        """
        Search Pfam database for domains in a sequence using hmmscan.

        Parameters
        ----------
        sequence_file : str, optional
            Path to FASTA file. If None, uses target sequence.

        Returns
        -------
        pd.DataFrame
            DataFrame with domain hits (name, accession, start, end, e-value)
        """
        if sequence_file is None:
            sequence_file = str(self.target_path)

        pfam_db = Path("databases/hmm/Pfam/Pfam-A.hmm")
        if not pfam_db.exists():
            self._log("Warning: Pfam database not found. Run Setup() first.")
            return pd.DataFrame()

        output_file = self.temp_dir / "pfam_domains.out"
        domtbl_file = self.temp_dir / "pfam_domains.domtbl"

        cmd = [
            "hmmscan",
            "--tblout", str(output_file),
            "--domtblout", str(domtbl_file),
            "-E", str(self.e_value_threshold),
            str(pfam_db),
            sequence_file
        ]
        self._run_command(cmd, "Running hmmscan against Pfam")

        return self._parse_hmmscan_output(domtbl_file)

    def _parse_hmmscan_output(self, domtbl_file: Path) -> pd.DataFrame:
        """Parse hmmscan domain table output."""
        if not domtbl_file.exists():
            return pd.DataFrame()

        domains = []
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
                        'e_value': float(parts[6]),
                        'score': float(parts[7]),
                        'start': int(parts[17]),
                        'end': int(parts[18]),
                        'ali_start': int(parts[19]),
                        'ali_end': int(parts[20])
                    })

        return pd.DataFrame(domains)

    def hmmsearch_with_profile(self, database_fasta: Optional[str] = None) -> pd.DataFrame:
        """
        Search a sequence database using the built HMM profile.

        Parameters
        ----------
        database_fasta : str, optional
            Path to FASTA database. If None, uses PDB FASTA.

        Returns
        -------
        pd.DataFrame
            DataFrame with homolog hits from HMM search
        """
        if not hasattr(self, 'hmm_profile') or not self.hmm_profile.exists():
            self._log("Warning: HMM profile not built. Run pipeline with run_hmmer=True.")
            return pd.DataFrame()

        if database_fasta is None:
            database_fasta = str(self.pdb_fasta)

        output_file = self.temp_dir / "hmmsearch_hits.out"
        tblout_file = self.temp_dir / "hmmsearch_hits.tbl"

        cmd = [
            "hmmsearch",
            "--tblout", str(tblout_file),
            "-E", str(self.e_value_threshold),
            str(self.hmm_profile),
            database_fasta
        ]
        self._run_command(cmd, "Running hmmsearch with profile")

        return self._parse_hmmsearch_output(tblout_file)

    def _parse_hmmsearch_output(self, tblout_file: Path) -> pd.DataFrame:
        """Parse hmmsearch table output."""
        if not tblout_file.exists():
            return pd.DataFrame()

        hits = []
        with open(tblout_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 9:
                    target_name = parts[0]
                    # Extract PDB ID from name (format: pdb|XXXX|chain)
                    if '|' in target_name:
                        pdb_parts = target_name.split('|')
                        pdb_id = f"{pdb_parts[1].upper()}_{pdb_parts[2]}" if len(pdb_parts) > 2 else pdb_parts[1]
                    else:
                        pdb_id = target_name[:6] if len(target_name) >= 6 else target_name

                    hits.append({
                        'PDB_ID': pdb_id,
                        'E-value': float(parts[4]),
                        'score': float(parts[5]),
                        'bias': float(parts[6]),
                        'description': ' '.join(parts[18:]) if len(parts) > 18 else ''
                    })

        df = pd.DataFrame(hits)
        if len(df) > 0:
            df = df.drop_duplicates(subset=['PDB_ID'])
            df = df.sort_values('E-value')
        return df

    def _parse_pdb_blast_output(self, output_file: Path, origin: str) -> pd.DataFrame:
        """Parse BLAST tabular output from PDB search."""
        if not output_file.exists() or output_file.stat().st_size == 0:
            return pd.DataFrame(columns=self.PDB_BLAST_COLUMNS + ["Protein Name", "Organism", "Origin"])

        df = pd.read_csv(output_file, sep="\t", header=None, names=self.PDB_BLAST_COLUMNS)

        # Parse protein name and organism from Subject_Title
        df[["Protein Name", "Organism"]] = df["Subject_Title"].apply(self._parse_pdb_title)
        df.drop(columns=["Subject_Title"], inplace=True)

        # Filter and deduplicate
        df = df[df["E-value"] < self.e_value_threshold].drop_duplicates(subset=["PDB_ID"])
        df.sort_values(by=["Bitscore", "E-value"], ascending=[False, True], inplace=True)
        df["Origin"] = origin

        return df

    @staticmethod
    def _parse_pdb_title(title: str) -> pd.Series:
        """Parse protein name and organism from PDB title string."""
        # Extract organism from brackets
        if '[' in title:
            parts = title.rsplit('[', 1)
            organism = parts[1].replace(']', '').strip()
            name_part = parts[0].strip()
        else:
            organism = "Unknown"
            name_part = title

        # Clean up chain prefix
        if name_part.startswith("Chain"):
            name_split = name_part.split(',', 1)
            clean_name = name_split[1].strip() if len(name_split) > 1 else name_part
        else:
            clean_name = name_part

        return pd.Series([clean_name.upper(), organism])

    def _rank_homologs(self):
        """Combine and rank all template hits."""
        self._log("Ranking all template hits...")

        # Combine PSI-BLAST and MSA hits
        dfs = [df for df in [self.top_psi_hits, self.top_msa_hits] if len(df) > 0]

        if not dfs:
            self.homologs = pd.DataFrame()
            return

        self.homologs = pd.concat(dfs, ignore_index=True)

        # Convert to numeric
        for col in ["E-value", "Query_Coverage", "Identity"]:
            if col in self.homologs.columns:
                self.homologs[col] = pd.to_numeric(self.homologs[col], errors='coerce')

        # Remove duplicates keeping best hit
        self.homologs = self.homologs.drop_duplicates(subset=["PDB_ID"], keep="first")

        # Sort by E-value, then Identity, then Coverage
        self.homologs.sort_values(
            by=["E-value", "Identity", "Query_Coverage"],
            ascending=[True, False, False],
            inplace=True
        )
        self.homologs.reset_index(drop=True, inplace=True)

    def plot_homologs(self, top_n: int = 10):
        """
        Plot distribution of top protein names and organisms.

        Parameters
        ----------
        top_n : int
            Number of top entries to show in each plot
        """
        if len(self.homologs) == 0:
            print("No homologs to plot.")
            return

        # Create copy with unique PDB IDs (ignoring chain)
        df = self.homologs.copy()
        df["PDB_4"] = df["PDB_ID"].str[:4]
        df = df.drop_duplicates(subset=["PDB_4"])

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Protein names
        name_counts = df["Protein Name"].value_counts().head(top_n)
        axes[0].barh(range(len(name_counts)), name_counts.values)
        axes[0].set_yticks(range(len(name_counts)))
        axes[0].set_yticklabels(name_counts.index, fontsize=9)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Count")
        axes[0].set_title(f"Top {top_n} Protein Types")

        # Organisms
        org_counts = df["Organism"].value_counts().head(top_n)
        axes[1].barh(range(len(org_counts)), org_counts.values)
        axes[1].set_yticks(range(len(org_counts)))
        axes[1].set_yticklabels(org_counts.index, fontsize=9)
        axes[1].invert_yaxis()
        axes[1].set_xlabel("Count")
        axes[1].set_title(f"Top {top_n} Organisms")

        plt.tight_layout()
        plt.show()

    def get_top_templates(self, n: int = 5) -> pd.DataFrame:
        """
        Get top N template candidates.

        Parameters
        ----------
        n : int
            Number of templates to return

        Returns
        -------
        pd.DataFrame
            Top template hits
        """
        return self.homologs.head(n).copy()

    def __repr__(self) -> str:
        return f"TemplateRetriever(found={len(self.homologs)} templates, e_value<{self.e_value_threshold})"
