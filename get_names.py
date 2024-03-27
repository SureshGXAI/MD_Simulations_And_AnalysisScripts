from Bio import PDB

pdb = "8fza_preQ0-complex.pdb"
name = pdb[:3]

p = PDB.PDBParser()
s = p.get_structure(name, pdb)

y = s.get_residues()
for x in y:
    print(x.get_id()[1], x.resname)
