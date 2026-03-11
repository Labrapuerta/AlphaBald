from Bio.PDB import Superimposer, PDBIO
from Bio import SeqIO
from Bio.PDB import PDBParser
import Path


def get_ca_mapping(pdb_id, msa_seq):
    pdb_id, chain_id = pdb_id.split("_")[0], pdb_id.split("_")[1]
    pdb_file = Path("Templates") / f"pdb{pdb_id.lower()}.ent"
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    chain = structure[0][chain_id]



    ca_atoms = [res['CA'] for res in chain.get_residues() if res.has_id('CA')]
    alignment = list(SeqIO.parse("Alignments/aligned_templates.fasta", "fasta"))
    print(alignment)
    aligned_target = next(rec.seq for rec in alignment if rec.id == pdb_id + "_" + chain_id)
    aligned_template = next(rec.seq for rec in alignment if rec.id == "GLB1_CALSO")

    mapping = {}
    target_idx = 0 
    atom_idx = 0   
    
    for t_char, p_char in zip(aligned_target, aligned_template):        
        if t_char != '-' and p_char != '-':
            if atom_idx < len(ca_atoms):
                mapping[target_idx] = ca_atoms[atom_idx]
        
        if t_char != '-':
            target_idx += 1 
        if p_char != '-':
            atom_idx += 1   
            
    return mapping, structure

def superimpose_structures(target_structure, template_structure, mapping):
    target_atoms = [mapping[i] for i in sorted(mapping.keys())]
    template_atoms = [mapping[i] for i in sorted(mapping.keys())]

    print(target_atoms)
    print(template_atoms)