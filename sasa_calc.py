######################################################
# SASA calculation for pdb structures
# packages required: freesasa, biopython
# python3 sasa_calc.py <pdb file>
######################################################

import os, sys
import freesasa
from Bio import PDB

inputfile = sys.argv[1]
outfile = open(inputfile+".txt", 'a')

# get the residue names and ids 
pdb = inputfile+".pdb"
name = pdb[:3]
p = PDB.PDBParser().get_structure(name, pdb)
TotalResidues = p.get_residues()

# Reading the structure file and calculating SASA
structure = freesasa.Structure(inputfile+".pdb")
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result, structure)

outfile.write("Total SASA: "+str(result.totalArea())+"\n")

for key in area_classes:
    outfile.write(str(key) +" : "+ str(area_classes[key])+"\n")


# Writing the pdb file with SASA as B-factor and Atomic radii as occupancy
result.write_pdb(inputfile+'.sasa.pdb')

#result = freesasa.calc(structure, freesasa.Parameters({'algorithm' : freesasa.LeeRichards, 'n-slices' : 100}))

# SASA for individual residues in the system
for res in TotalResidues:
    idx = res.get_id()[1]
    rdx = res.resname
    selections = freesasa.selectArea(('%s, resn %s' %(rdx, rdx), '%d, resi %d' %(idx, idx)), structure, result)
    for key1 in selections:
        if key1.isdigit():
            outid = key1 
            outsasa = selections[key1]
        else:
            outresn = key1
    outfile.write(str(outid)+" "+str(outresn)+" "+str(outsasa)+"\n")
        
outfile.close()
