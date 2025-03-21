"""
process_pdb.py
18 February 2025

Adds chain information to pdb files 
python3 process_pdb.py <input PDB file path>

"""
from Bio.PDB import PDBParser, PDBIO, Select
import pandas as pd
import sys
import os

import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', DeprecationWarning)

input_pdb=sys.argv[1]  # Path of input pdb file
pdb = os.path.splitext(os.path.basename(input_pdb))[0]
output_pdb = os.path.dirname(input_pdb)+"/"+pdb+"_out.pdb"

# print(input_pdb)
# print(pdb)
# print(output_pdb)

def add_chain_information(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        current_chain = 'A'
        previous_residue_number = None

        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_number = int(line[22:26].strip())
                if previous_residue_number is not None and residue_number == 1 and previous_residue_number - residue_number > 0:
                    current_chain = chr(ord(current_chain) + 1)
                previous_residue_number = residue_number
                line = line[:21] + current_chain + line[22:]

            if line.startswith("TER"):
                current_chain = "A"
                previous_residue_number = None
            outfile.write(line)

add_chain_information(input_pdb,output_pdb)




