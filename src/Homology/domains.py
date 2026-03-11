import requests
import os
import subprocess
import urllib.request
import gzip
import shutil
import pathlib
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBList

class TemplateProcessor:
    #### This is to see the domains of the templates and superimpose them. Later
    def __init__(self, templates, number_of_templates=2):
        self.data = templates
        self.templates_dir = Path("Templates")
        self.target_path = Path("target") / os.listdir("target")[0]
        self.top_templates = templates.head(number_of_templates)
        self.hmm_database_path = Path("Databases") / "hmm" / "Pfam"
        self.setup_pfam_database() ## Ensure the Pfam database is set up and ready to use
        #self.download_templates() ## Download the PDB files for the selected number of templates
        self.fasta_for_alignment = Path("temp") / "templates.fasta" # Fasta file for the target and templates to be aligned)
        #self.aligned_sequences = Path("alignments") / "aligned_templates.fasta" Fasta file for the aligned sequences of the target and templates
        

        self.pdb_sequences = {}
        for _, template in self.top_templates.iterrows():
            full_id = template['PDB_ID']
            pdb_4 = full_id[:4].lower()
            chain = full_id[5]
            # Find the downloaded file
            pdb_file = self.templates_dir / f"pdb{pdb_4}.ent"
            if pdb_file.exists():
                seq = self.parse_template_sequence(pdb_file, chain)
                if seq:
                    self.pdb_sequences[full_id] = seq # Key is now "9H0M_B"
        self.target_fasta = SeqIO.read(self.target_path, "fasta") # Read the target fasta file

        self.target_sequence = self.target_fasta.seq # Get the target sequence from the fasta file


        self.create_fasta_file()
        self.domain_data = []
        for full_id, sequence in self.pdb_sequences.items():
            pdb_id = full_id[:4].lower()
            domains = self.get_pdb_domains(pdb_id)
            temp = [full_id, sequence, domains]
            self.domain_data.append(temp)
    

    def get_pdb_domains(self, pdb_id):
        chain = pdb_id.split("_")[-1] if "_" in pdb_id else "A"  # Default to chain A if not specified
        pdb_id = pdb_id.lower()[0:4]
        url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/interpro/{pdb_id}"


        try:
            response = requests.get(url)
            if response.status_code != 200:
                print(f"Warning: No domain data found for {pdb_id}")
                return []
                
            data = response.json()
            domains = []
            
            # Parse the PDBe JSON structure
            if pdb_id.lower() in data:
                interpro_mappings = data[pdb_id.lower()].get('InterPro', {})
                
                for domain_id, domain_info in interpro_mappings.items():
                    domain_name = domain_info.get('name', domain_id)
                    mappings = domain_info.get('mappings', [])
                    
                    for mapping in mappings:
                        # Only grab domains that belong to the specific chain we are using!
                        if mapping.get('chain_id') == chain:
                            domains.append({
                                'id': domain_id,
                                'name': domain_name,
                                'start': mapping.get('start', {}).get('residue_number'),
                                'end': mapping.get('end', {}).get('residue_number')
                            })
            return domains
            
        except Exception as e:
            print(f"Error fetching PDB domains: {e}")
            return []
        
    def download_templates(self):
        if not self.templates_dir.exists():
            self.templates_dir.mkdir()
        top_templates = self.top_templates["PDB_ID"].str[:4].tolist()  
        pdbl = PDBList(verbose=False)
        for pdb_id in top_templates:
            pdbl.retrieve_pdb_file(pdb_id, pdir=self.templates_dir, file_format="pdb", overwrite=True)  
        return {pdb_id: self.templates_dir / f"pdb{pdb_id.lower()}.ent" for pdb_id in top_templates}

    def parse_template_sequence(self, pdb_file, chain_id):
        for record in SeqIO.parse(pdb_file, "pdb-atom"):
            chain_letter = record.id.split(':')[-1]
            if chain_letter == chain_id:
                return str(record.seq)
        return None
    
    def create_fasta_file(self):
        with open(self.fasta_for_alignment, 'w') as fasta_file:
            fasta_file.write(f">{self.target_fasta.id}\n{self.target_sequence}\n")
            # FIX 3: Write the full ID to the fasta file
            for full_id, sequence in self.pdb_sequences.items():
                fasta_file.write(f">{full_id}\n{sequence}\n")

    def setup_pfam_database(self):
        hmm_file = os.path.join(self.hmm_database_path, "Pfam-A.hmm")
        gz_file = hmm_file + ".gz"
        pressed_check = hmm_file + ".h3m"
        if os.path.exists(pressed_check):
            return hmm_file
                
        os.makedirs(self.hmm_database_path, exist_ok=True)
        
        if not os.path.exists(hmm_file):
            if not os.path.exists(gz_file):
                url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
                
                urllib.request.urlretrieve(url, gz_file)
            with gzip.open(gz_file, 'rb') as f_in:
                with open(hmm_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            os.remove(gz_file) 
            
        try:
            subprocess.run(["hmmpress", hmm_file], check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            raise
        except FileNotFoundError:
            raise