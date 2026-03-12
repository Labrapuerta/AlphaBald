import os
import subprocess
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBList
import matplotlib.pyplot as plt

class Retrieve:
    def __init__(self, target_path = "target", e_value_threshold=1e-5, num_iterations=5, num_hits=10):
        self.target_file = os.listdir(target_path)[0]
        self.target_path = Path(target_path) / self.target_file 
        self.psmm = "temp/target.pssm"
        self.temp_dir = Path("temp")        ### Temporary files will be stored here
        self.templates_dir = Path("Templates") ### PDB files will be downloaded here
        self.target_alignment_fasta = Path("temp/unaligned_homologs.fasta") ## Target on top, followed by the swissprot hits used for the MSA-based search
        self.aligned_homologs = Path("alignments/aligned_homologs.fa") ## MSA of the target and the swissprot hits used for the MSA-based search
        self.swissprot_db = Path("databases/swissprot")
        self.pdb_db = Path("databases/pdb_seq")
        self.pdb_fasta_db = Path("databases/pdb_seq/pdbaa.fasta")

        self.database_exists(self.swissprot_db, "swissprot.pin", "swissprot")
        self.database_exists(self.pdb_db, "pdbaa.pin", "pdbaa")

        self.generate_PSMM(num_iterations=num_iterations)
        self.retrieve_homologs(e_value_threshold=e_value_threshold)
        self.create_alignment(num_hits=num_hits)
        self.run_jackhmmer(num_hits=num_hits)
        self.manual_hhmer(num_hits=num_hits)
        self.rank_homologs()


    def setup_directories(self):
        os.makedirs("temp", exist_ok=True)
        os.makedirs("Templates", exist_ok=True)
        os.makedirs("Alignments", exist_ok=True)
        os.makedirs("databases/swissprot", exist_ok=True)
        os.makedirs("databases/pdb_seq", exist_ok=True)
        os.makedirs("databases/hmm/Pfam", exist_ok=True)

    def database_exists(self, db_path, db_file, db):
        if (db_path / db_file).exists():
            return True
        else:
            os.makedirs(db_path, exist_ok=True)
            subprocess.run(["update_blastdb.pl", "--decompress", db], cwd=db_path)
            return False
        

    def generate_PSMM(self, num_iterations=5):
        if not self.temp_dir.exists():
            self.temp_dir.mkdir()
        # Generate PSMM for the target sequence
        psiblast_cmd = [
            "psiblast",
            "-query", str(self.target_path),
            "-db", str(self.swissprot_db / "swissprot"),
            "-num_iterations", str(num_iterations),
            "-out_pssm", str(self.psmm),
            "-outfmt", "6 sacc bitscore evalue",
            "-out", str(self.temp_dir / "target_swissprot.out")
        ]
        subprocess.run(psiblast_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE )
        self.top_swissprot_hits = pd.read_csv(self.temp_dir / "target_swissprot.out", sep="\t", header=None, names=["Accession", "Bitscore", "E-value"])
        self.top_swissprot_hits = self.top_swissprot_hits[self.top_swissprot_hits["E-value"] < 1e-5]
        self.top_swissprot_hits = self.top_swissprot_hits.drop_duplicates(subset=["Accession"])
        self.top_swissprot_hits.sort_values(by=["Bitscore", "E-value"], inplace=True, ascending= [False, True])

    def retrieve_homologs(self, e_value_threshold=1e-5):
        # Search for homologous sequences in the PDB database using PSI-BLAST
        psiblast_cmd = [
            "psiblast",
            "-in_pssm", str(self.psmm),
            "-db", str(self.pdb_db / "pdbaa"),
            "-out", str(self.temp_dir / "target_pdbaa.out"),
            "-outfmt", "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"]
        
        subprocess.run(psiblast_cmd, stderr =subprocess.PIPE, stdout=subprocess.PIPE)
        self.top_psi_hits = pd.read_csv(self.temp_dir / "target_pdbaa.out", sep="\t", header=None, names=["PDB_ID", "Bitscore", "E-value", "Identity", "Query_Coverage", "Query_Start", "Query_End", "Subject_Start", "Subject_End", "Query_Length", "Subject_Identity"])
        self.top_psi_hits[["Protein Name", "Organism"]] = self.top_psi_hits["Subject_Identity"].apply(self.parse_pdb_title)
        self.top_psi_hits.drop(columns=["Subject_Identity"], inplace=True)
        self.top_psi_hits = self.top_psi_hits[self.top_psi_hits["E-value"] < e_value_threshold]
        self.top_psi_hits = self.top_psi_hits.drop_duplicates(subset=["PDB_ID"])
        self.top_psi_hits.sort_values(by=["Bitscore", "E-value"], inplace=True, ascending=[False, True])
        self.top_psi_hits["Origin"] = "Automatic PSI-BLAST Search"

    def create_alignment(self, num_hits=10):

        os.makedirs("Alignments", exist_ok=True)

        if self.target_alignment_fasta.exists():
            os.remove(self.target_alignment_fasta)

        with open(self.target_alignment_fasta, "w") as fasta_out:
            with open(Path("target") / self.target_file, "r") as target_fasta:
                target_content = target_fasta.read()
                fasta_out.write(target_content)
                if not target_content.endswith("\n"):
                    fasta_out.write("\n")

        selected_hits = self.top_swissprot_hits[:num_hits]["Accession"].to_list()
        selected_hits = ",".join(selected_hits)

        cmd = [
            "blastdbcmd",
            "-db", str(self.swissprot_db / "swissprot"),
            "-entry", selected_hits
        ]
        x = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with open(self.target_alignment_fasta, "a") as fasta_out:
            fasta_out.write(x.stdout.decode())

        cmd = [
            "clustalw",
            "-INFILE=" + str(self.target_alignment_fasta),
            "-OUTFILE=" + str(self.aligned_homologs),
            "-OUTPUT=FASTA", "-OUTORDER=INPUT"
        ]
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        cmd = [
            "psiblast",
            "-in_msa", str(self.aligned_homologs),
            "-db", str(self.pdb_db / "pdbaa"),
            "-out", str(self.temp_dir / "msa_top_hits.out"),
            "-outfmt", "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
        ]

        subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        self.top_msa_hits = pd.read_csv(self.temp_dir / "msa_top_hits.out", sep="\t", header=None, names=["PDB_ID", "Bitscore", "E-value", "Identity", "Query_Coverage", "Query_Start", "Query_End", "Subject_Start", "Subject_End",  "Query_Length", "Subject_Identity"])
        self.top_msa_hits[["Protein Name", "Organism"]] = self.top_msa_hits["Subject_Identity"].apply(self.parse_pdb_title)
        self.top_msa_hits.drop(columns=["Subject_Identity"], inplace=True)
        self.top_msa_hits = self.top_msa_hits[self.top_msa_hits["E-value"] < 1e-5]
        self.top_msa_hits = self.top_msa_hits.drop_duplicates(subset=["PDB_ID"])
        self.top_msa_hits.sort_values(by=["Bitscore", "E-value"], inplace=True, ascending=[False, True])
        self.top_msa_hits["Origin"] = "MSA-based Search"


    def run_jackhmmer(self, num_hits=10):
        ### Create the FASTA for the database search
        if not self.pdb_fasta_db.exists():
            cmd = ["blastdbcmd", "-db", str(self.pdb_db / "pdbaa"), "-entry", "all", "-out", str(self.pdb_db / "pdbaa.fasta")]
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        jackhammer_cmd = ["jackhmmer", "-N", str(num_hits), "--tblout", str(self.temp_dir / "jackhmmer_hits.out"), str(self.target_path), str(self.pdb_fasta_db)]
        subprocess.run(jackhammer_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    def manual_hhmer(self, num_hits=10):
        hmmbuild = ["hmmbuild", str(self.temp_dir / "target.hmm"), str(self.aligned_homologs)]
        subprocess.run(hmmbuild, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



    def rank_homologs(self):
        self.homologs = pd.concat([self.top_psi_hits, self.top_msa_hits], ignore_index=True)

        self.homologs["E-value"] = pd.to_numeric(self.homologs["E-value"])
        self.homologs["Query_Coverage"] = pd.to_numeric(self.homologs["Query_Coverage"])
        self.homologs["Identity"] = pd.to_numeric(self.homologs["Identity"])
        self.homologs.sort_values(by=["E-value", "Identity", "Query_Coverage"], inplace=True, ascending=[True, False, False])
        self.homologs.reset_index(drop=True, inplace=True)
    def parse_pdb_title(self, title):

        if '[' in title:
            parts = title.rsplit('[', 1) # Split from the right side
            organism = parts[1].replace(']', '').strip()
            name_part = parts[0].strip() 
        else:
            organism = "Unknown"
            name_part = title

        if name_part.startswith("Chain"):
            name_split = name_part.split(',', 1)
            if len(name_split) > 1:
                clean_name = name_split[1].strip()
            else:
                clean_name = name_part
        else:
            clean_name = name_part

        return pd.Series([clean_name.upper(), organism])
    
    def plot_homologs(self):
        self.homologs_unique = self.homologs.copy()
        self.homologs_unique["PDB_ID"] = self.homologs["PDB_ID"].str.split("_").str[0]
        self.homologs_unique.drop_duplicates(subset=["PDB_ID"], inplace=True)
        fig, axes = plt.subplots(1, 2, figsize=(17, 5))
        counts = self.homologs_unique["Protein Name"].value_counts().head(10)
        axes[0].bar(counts.index, counts.values)
        axes[0].tick_params(axis="x", rotation=45)
        axes[0].set_title("Top 10 Proteins")
        axes[0].set_ylabel("Count")

        organism_counts = self.homologs_unique["Organism"].value_counts().head(10)
        axes[1].bar(organism_counts.index, organism_counts.values)
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].set_title("Top 10 Organisms")
        axes[1].set_ylabel("Count")
        plt.show()
