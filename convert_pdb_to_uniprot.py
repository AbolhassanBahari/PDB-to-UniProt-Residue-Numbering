from Bio import PDB
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Function to read the UniProt sequence from a FASTA file
def read_uniprot_sequence(fasta_file):
    with open(fasta_file, 'r') as file:
        seq_record = SeqIO.read(file, "fasta")
        return str(seq_record.seq)

# Function to extract the sequence from the PDB file
def extract_pdb_sequence(pdb_file):
    # Load the PDB structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB', pdb_file)
    
    pdb_sequence = ""
    
    # Extracting sequence for the first chain (assuming single chain or first chain)
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                pdb_sequence = str(poly.get_sequence())
                
    return pdb_sequence

# Function to align the PDB sequence with UniProt sequence using PairwiseAligner
def align_sequences(pdb_sequence, uniprot_sequence):
    aligner = PairwiseAligner()
    alignments = aligner.align(pdb_sequence, uniprot_sequence)
    
    # Return the best alignment (first alignment)
    return alignments[0]

# Function to map PDB residue numbers to UniProt sequence numbers based on alignment
def convert_pdb_to_uniprot(pdb_file, uniprot_fasta_file):
    # Read UniProt sequence from the provided FASTA file
    uniprot_sequence = read_uniprot_sequence(uniprot_fasta_file)

    # Extract PDB sequence
    pdb_sequence = extract_pdb_sequence(pdb_file)

    # Align the sequences
    alignment = align_sequences(pdb_sequence, uniprot_sequence)
    
    # Print alignment details (optional)
    print(alignment)

    # Create a dictionary to store mapping of PDB index to UniProt position
    pdb_to_uniprot_mapping = {}

    # Iterate over the alignment and create the correct mapping
    pdb_idx, uniprot_idx = 0, 0
    for pdb_residue, uniprot_residue in zip(alignment[0], alignment[1]):
        if pdb_residue != '-':  # If it's not a gap in the PDB sequence
            if uniprot_residue != '-':  # Ensure we are not mapping a gap in UniProt
                pdb_to_uniprot_mapping[pdb_idx] = uniprot_idx
            pdb_idx += 1
        if uniprot_residue != '-':
            uniprot_idx += 1

    # Load the PDB structure again to modify residue numbers
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB', pdb_file)

    # Iterate over chains and residues in the PDB structure
    residue_idx = 0
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            
            for poly in polypeptides:
                for residue in poly:
                    if residue_idx in pdb_to_uniprot_mapping:
                        # Get the corresponding UniProt residue number
                        uniprot_pos = pdb_to_uniprot_mapping[residue_idx] + 1  # +1 for 1-based indexing
                        
                        # Modify the residue number in the PDB to match the UniProt sequence number
                        residue.id = (' ', uniprot_pos, ' ')
                        residue_idx += 1

    # Output the updated PDB file with new residue numbering
    io = PDB.PDBIO()
    io.set_structure(structure)
    output_file = f"converted_{pdb_file}"
    io.save(output_file)  # Saving the file with new name
    print(f"Updated PDB file saved as {output_file}")

# Usage
pdb_file = '/mnt/data/5ikr_renumbered.pdb'  # Path to your PDB file
uniprot_fasta_file = '/mnt/data/P35354.fasta'  # Path to your UniProt FASTA file
convert_pdb_to_uniprot(pdb_file, uniprot_fasta_file)
