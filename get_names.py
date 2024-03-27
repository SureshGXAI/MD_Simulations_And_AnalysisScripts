from Bio import PDB
import os, sys
inpfile = sys.argv[1]
pdb = inpfile + ".pdb"
name = pdb[:3]

p = PDB.PDBParser()
s = p.get_structure(name, pdb)

y = s.get_residues()
for x in y:
    print(x.get_id()[1], x.resname)
