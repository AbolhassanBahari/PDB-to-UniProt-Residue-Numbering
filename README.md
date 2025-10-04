# PDB-to-UniProt-Residue-Numbering
A simple Python script that renumbers the residues of a PDB file to match the corresponding residue numbers in the UniProt sequence.
It uses sequence alignment to ensure that the residue numbers in the PDB file align correctly with the UniProt sequence, making it easier to compare PDB structures with their corresponding protein sequences.

# Features:

Extracts sequences from PDB files and UniProt FASTA files.
Aligns PDB sequence with UniProt sequence using pairwise sequence alignment.
Renumbers the residues in the PDB file based on the UniProt sequence.
Handles gaps in sequences properly (e.g., insertions or deletions).
Outputs a new PDB file with correctly aligned residue numbers.

# Requirements:

Python 3.x
Biopython library

You can install the required dependencies using pip:

```pip install biopython```

# How It Works:

Extracts Sequence from PDB: The script extracts the sequence of residues from a PDB file. It assumes a single chain or the first chain in the PDB file.

Extracts UniProt Sequence: The UniProt sequence is provided in FASTA format.

Sequence Alignment: The PDB and UniProt sequences are aligned using the PairwiseAligner from Biopython's Align module.

Residue Numbering: The script updates the residue numbers in the PDB file to match the numbering of the aligned UniProt sequence.

Outputs Updated PDB: The PDB file is saved with the new residue numbering.

# Usage:

Download the necessary files:

PDB File: The PDB file that you want to modify.

UniProt FASTA File: The FASTA file containing the UniProt sequence of the protein.

Run the Script:

Ensure you have Python and the required libraries installed (Biopython).

Save the script to a Python file (e.g., convert_pdb_to_uniprot.py).

Run the script as follows:
```python convert_pdb_to_uniprot.py```

Script Parameters:

```pdb_file: The path to the input PDB file.```


```uniprot_fasta_file: The path to the FASTA file containing the UniProt sequence.```

# Notes:

The script currently handles only the first chain in the PDB file. If you have a multi-chain PDB file, modifications may be necessary to handle multiple chains.

The sequence alignment may have gaps (-) in the aligned sequences, especially in cases of insertions or deletions in one sequence compared to the other. The script handles these gaps and aligns the residue numbers correctly.

# License:

This project is licensed under the MIT License
