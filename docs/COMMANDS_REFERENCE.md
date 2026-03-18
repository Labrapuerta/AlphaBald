# Bioinformatics Commands Reference

Complete reference for all commands used in the AlphaBald homology modeling pipeline.
All paths are relative to the project root (`AlphaBald/`).

## Database Paths

| Database | Path | Description |
|----------|------|-------------|
| SwissProt | `databases/swissprot/swissprot` | Curated protein sequences |
| PDB sequences | `databases/pdb_seq/pdbaa` | PDB protein sequences |
| PDB FASTA | `databases/pdb_seq/pdbaa.fasta` | FASTA version for HMMER |
| Pfam HMM | `databases/hmm/Pfam/Pfam-A.hmm` | Protein family HMM profiles |

---

## 1. Database Setup

### Download SwissProt BLAST Database
```bash
cd databases/swissprot/
update_blastdb.pl --decompress swissprot
```

### Download PDB Sequence Database
```bash
cd databases/pdb_seq/
update_blastdb.pl --decompress pdbaa
```

### Download and Index Pfam Database
```bash
cd databases/hmm/Pfam/
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### Convert BLAST Database to FASTA (for HMMER)
```bash
blastdbcmd -db databases/pdb_seq/pdbaa \
           -entry all \
           -out databases/pdb_seq/pdbaa.fasta
```

---

## 2. BLAST Commands

### Basic BLASTP Search
```bash
blastp -query target/target.fa \
       -db databases/swissprot/swissprot \
       -out temp/blastp_results.out \
       -outfmt "6 sacc bitscore evalue pident"
```

### PSI-BLAST: Build PSSM from SwissProt
```bash
psiblast -query target/target.fa \
         -db databases/swissprot/swissprot \
         -num_iterations 5 \
         -out_pssm temp/target.pssm \
         -outfmt "6 sacc bitscore evalue" \
         -out temp/target_swissprot.out
```

### PSI-BLAST: Search PDB with PSSM
```bash
psiblast -in_pssm temp/target.pssm \
         -db databases/pdb_seq/pdbaa \
         -out temp/target_pdbaa.out \
         -outfmt "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
```

### PSI-BLAST: Search with MSA as Input
```bash
psiblast -in_msa Alignments/aligned_homologs.fa \
         -db databases/pdb_seq/pdbaa \
         -out temp/msa_top_hits.out \
         -outfmt "6 sacc bitscore evalue pident qcovs qstart qend sstart send qlen stitle"
```

### Extract Sequences from BLAST Database
```bash
# Single sequence
blastdbcmd -db databases/swissprot/swissprot \
           -entry P12345

# Multiple sequences
blastdbcmd -db databases/swissprot/swissprot \
           -entry "P12345,Q67890,R11111"

# All sequences to FASTA
blastdbcmd -db databases/pdb_seq/pdbaa \
           -entry all \
           -out output.fasta
```

### BLAST Output Format Codes
| Code | Field |
|------|-------|
| sacc | Subject accession |
| bitscore | Bit score |
| evalue | E-value |
| pident | Percentage identity |
| qcovs | Query coverage |
| qstart | Query start |
| qend | Query end |
| sstart | Subject start |
| send | Subject end |
| qlen | Query length |
| stitle | Subject title |

---

## 3. HMMER Commands

### Jackhmmer: Iterative Search
```bash
jackhmmer -N 5 \
          --tblout temp/jackhmmer_hits.out \
          target/target.fa \
          databases/pdb_seq/pdbaa.fasta
```

Options:
- `-N 5`: Number of iterations
- `--tblout`: Table output file
- `--domtblout`: Domain table output
- `-E 0.001`: E-value threshold

### Build HMM Profile from Alignment
```bash
hmmbuild temp/target.hmm Alignments/aligned_homologs.fa
```

Options:
- `--amino`: Force amino acid alphabet
- `-n NAME`: Name the HMM profile

### Search with HMM Profile
```bash
hmmsearch --tblout temp/hmmsearch_hits.out \
          temp/target.hmm \
          databases/pdb_seq/pdbaa.fasta
```

### Scan Sequence Against HMM Database (Pfam)
```bash
hmmscan --tblout temp/pfam_hits.out \
        databases/hmm/Pfam/Pfam-A.hmm \
        target/target.fa
```

### Index HMM Database
```bash
hmmpress databases/hmm/Pfam/Pfam-A.hmm
```

---

## 4. Multiple Sequence Alignment

### ClustalW Alignment
```bash
clustalw -INFILE=temp/unaligned_homologs.fasta \
         -OUTFILE=Alignments/aligned_homologs.fa \
         -OUTPUT=FASTA \
         -OUTORDER=INPUT
```

Options:
- `-OUTPUT=FASTA`: Output format (FASTA, CLUSTAL, PIR)
- `-OUTORDER=INPUT`: Keep input order
- `-GAPOPEN=10`: Gap opening penalty
- `-GAPEXT=0.2`: Gap extension penalty

### MUSCLE Alignment (alternative)
```bash
muscle -in temp/unaligned.fasta \
       -out Alignments/aligned.fa
```

### T-Coffee Alignment (alternative)
```bash
t_coffee temp/unaligned.fasta \
         -output fasta \
         -outfile Alignments/aligned.fa
```

---

## 5. Structure Analysis

### DSSP: Secondary Structure Assignment
```bash
# From PDB file
dssp -i Templates/pdb1abc.ent -o temp/1abc.dssp

# Short output format
dssp -i structure.pdb --output-format dssp
```

DSSP Output Codes:
| Code | Structure |
|------|-----------|
| H | Alpha helix |
| B | Beta bridge |
| E | Extended strand |
| G | 3-10 helix |
| I | Pi helix |
| T | Turn |
| S | Bend |
| - | Coil |

### ProSA Web Server
Submit structures at: https://prosa.services.came.sbg.ac.at/prosa.php

---

## 6. MODELLER Commands

### Basic Model Building Script
```python
from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

# Directories for input atom files
env.io.atom_files_directory = ['Templates/', '.']

# Build model
a = AutoModel(env,
              alnfile='Alignments/target_template.ali',
              knowns='1ABC_A',           # Template PDB ID
              sequence='target',          # Target sequence ID
              assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5
a.make()
```

### PIR Alignment Format for MODELLER
```
>P1;template
structureX:1ABC_A:1:A:217:A::::
MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIE...*

>P1;target
sequence:target:::::::0.00:0.00
MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIE...*
```

---

## 7. Visualization Commands

### PyMOL
```bash
# Open structure
pymol Templates/pdb1abc.ent

# Superimpose structures
pymol -c -d "load struct1.pdb; load struct2.pdb; align struct1, struct2; save aligned.pdb"
```

### PyMOL Python Commands
```python
import pymol
from pymol import cmd

cmd.load("structure.pdb")
cmd.show("cartoon")
cmd.color("blue", "ss h")  # Color helices
cmd.color("red", "ss s")   # Color sheets
cmd.save("output.png")
```

---

## 8. Common Workflows

### Complete Template Search Pipeline
```bash
# 1. Build PSSM
psiblast -query target/target.fa \
         -db databases/swissprot/swissprot \
         -num_iterations 5 \
         -out_pssm temp/target.pssm \
         -out temp/swissprot_hits.out

# 2. Search PDB
psiblast -in_pssm temp/target.pssm \
         -db databases/pdb_seq/pdbaa \
         -out temp/pdb_hits.out \
         -outfmt "6 sacc bitscore evalue pident qcovs"

# 3. Extract sequences for MSA
blastdbcmd -db databases/swissprot/swissprot \
           -entry "ACC1,ACC2,ACC3" > temp/homologs.fasta

# 4. Create MSA
clustalw -INFILE=temp/homologs.fasta \
         -OUTFILE=Alignments/aligned.fa \
         -OUTPUT=FASTA

# 5. Build HMM profile
hmmbuild temp/profile.hmm Alignments/aligned.fa

# 6. Search with HMM
hmmsearch temp/profile.hmm databases/pdb_seq/pdbaa.fasta
```

### Domain Analysis Pipeline
```bash
# Scan against Pfam
hmmscan --tblout temp/domains.out \
        databases/hmm/Pfam/Pfam-A.hmm \
        target/target.fa

# Parse results
grep -v "^#" temp/domains.out | awk '{print $1, $3, $4, $7}'
```

---

## 9. File Format References

### FASTA Format
```
>sequence_id description
MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIE
MKPHPWFFGKIPRAKAEEMLSKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHF
```

### Stockholm Format (for HMMER)
```
# STOCKHOLM 1.0
seq1    ACDEFGHIKLMNPQRSTVWY
seq2    ACDEFGHIKLMNPQRSTVWY
//
```

### PIR Format (for MODELLER)
```
>P1;sequence_name
type:code:start:chain:end:chain:name:source:resolution:r-factor
SEQUENCE*
```

---

## 10. Useful One-Liners

```bash
# Count sequences in FASTA
grep -c "^>" file.fasta

# Extract sequence IDs
grep "^>" file.fasta | sed 's/>//'

# Get sequence length
awk '/^>/ {if (seq) print name, length(seq); name=$0; seq=""} /^[^>]/ {seq=seq$0} END {print name, length(seq)}' file.fasta

# Convert multiline FASTA to single line
awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);N++;next;} {printf("%s",$0)} END {printf("\n")}' file.fasta
```
