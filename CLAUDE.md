# AlphaBald - Homology Modeling Pipeline

## Project Overview
This is a homology modeling pipeline for protein structure prediction. The pipeline searches for template structures, performs alignments, and prepares data for MODELLER-based model building.

## Quick Start
```bash
conda activate AlphaBald
# Then run The_Fate_of_Baldo.ipynb
```

## Project Structure
```
AlphaBald/
├── src/
│   ├── pipeline.py          # Main HomologyPipeline class
│   ├── Setup/
│   │   ├── setup.py          # Database and directory setup
│   │   └── preprocessing.py  # Target file preprocessing
│   ├── Homology/
│   │   ├── retrieve.py       # Template search (PSI-BLAST, HMMER)
│   │   ├── domains.py        # Template domain analysis
│   │   └── superimpose.py    # Structure superimposition
│   ├── UI/
│   │   └── app.py            # Interactive visualization widgets
│   ├── modeller/             # NEW: MODELLER script generators
│   │   └── scripts.py        # Model building, refinement, evaluation
│   └── Analysis/             # NEW: Structure assessment tools
│       ├── assessment.py     # Family ID, DSSP, cation placement
│       └── visualization.py  # PyMOL figure generation
├── target/                   # Place your target sequence here (target.fa)
├── databases/                # BLAST and HMM databases
├── Templates/                # Downloaded PDB template files
├── Modeller_Templates/       # Cropped PDBs and PIR files for MODELLER
├── Alignments/               # Output alignments
├── Problems/                 # Practice problems (e.g., Problem_3)
└── temp/                     # Temporary files
```

## Pipeline Commands Reference

### 1. Database Setup
```bash
# SwissProt BLAST database
update_blastdb.pl --decompress swissprot

# PDB sequence database
update_blastdb.pl --decompress pdbaa

# Pfam HMM database
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### 2. Template Search - PSI-BLAST
```bash
# Build PSSM from SwissProt search (5 iterations)
psiblast -query target/target.fa \
         -db databases/swissprot/swissprot \
         -num_iterations 5 \
         -out_pssm temp/target.pssm \
         -outfmt "6 sacc bitscore evalue" \
         -out temp/target_swissprot.out

# Search PDB using PSSM
psiblast -in_pssm temp/target.pssm \
         -db databases/pdb_seq/pdbaa \
         -out temp/target_pdbaa.out \
         -outfmt "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
```

### 3. Multiple Sequence Alignment
```bash
# Extract sequences from BLAST database
blastdbcmd -db databases/swissprot/swissprot -entry acc1,acc2,acc3

# Align with ClustalW
clustalw -INFILE=temp/unaligned_homologs.fasta \
         -OUTFILE=Alignments/aligned_homologs.fa \
         -OUTPUT=FASTA \
         -OUTORDER=INPUT

# MSA-based PSI-BLAST search
psiblast -in_msa Alignments/aligned_homologs.fa \
         -db databases/pdb_seq/pdbaa \
         -out temp/msa_top_hits.out \
         -outfmt "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
```

### 4. HMMER Searches
```bash
# Convert BLAST DB to FASTA for jackhmmer
blastdbcmd -db databases/pdb_seq/pdbaa -entry all -out databases/pdb_seq/pdbaa.fasta

# Jackhmmer iterative search
jackhmmer -N 10 \
          --tblout temp/jackhmmer_hits.out \
          target/target.fa \
          databases/pdb_seq/pdbaa.fasta

# Build HMM profile from alignment
hmmbuild temp/target.hmm Alignments/aligned_homologs.fa

# Search with HMM profile
hmmsearch temp/target.hmm databases/pdb_seq/pdbaa.fasta
```

### 5. Structure Analysis
```bash
# DSSP for secondary structure
dssp -i structure.pdb -o structure.dssp

# ProSA for energy analysis (via web or local)
# Submit at: https://prosa.services.came.sbg.ac.at/prosa.php
```

## Python Pipeline Usage
```python
from src.pipeline import HomologyPipeline

# Full pipeline
pipeline = HomologyPipeline(
    target_path="target/target.fa",
    e_value_threshold=1e-5,
    num_iterations=5,
    num_templates=5
)

# Access results
pipeline.homologs           # DataFrame of ranked templates
pipeline.templates          # TemplateProcessor with domain info

# Visualization
from src.UI.app import CoverageVisualizer, DomainVisualizer
CoverageVisualizer(pipeline.homologs, number_of_templates=3).show()
```

## Key Dependencies
- BLAST+ (psiblast, blastdbcmd, update_blastdb.pl)
- HMMER (hmmbuild, hmmsearch, jackhmmer, hmmpress)
- ClustalW
- MODELLER
- BioPython
- PyMOL (optional, for visualization)

## Common Issues

### "Identity" column missing error
The BLAST output format must include `pident` for identity percentage.

### Database not found
Run `Setup()` first or manually download databases to `databases/` folder.

### ClustalW alignment fails
Ensure sequences in FASTA have unique IDs and no special characters.

## Conventions
- Target sequence goes in `target/target.fa`
- PDB IDs use format: `XXXX_Y` (4-letter code + chain)
- E-value threshold default: 1e-5
- Identity and coverage are in percentages (0-100)
