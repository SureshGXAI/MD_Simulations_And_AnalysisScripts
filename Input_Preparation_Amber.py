#!/usr/bin/env python

import os, sys
from prody import *
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO
import pypdb


def get_pdb_components(pdb_id):
    pdb = parsePDB(pdb_id)
    protein = pdb.select('protein')
    ligand = pdb.select('not protein and not water')
    return protein, ligand


def process_ligand(ligand, res_name):
    output = StringIO()
    sub_mol = ligand.select(f"resname {res_name}")
    chem_desc = pypdb.describe_chemical(f"{res_name}")
    sub_smiles = None
    for item in chem_desc.get('pdbx_chem_comp_descriptor', []):
        if item.get('type') == 'SMILES':
            sub_smiles = item.get('descriptor')
            break
    template = AllChem.MolFromSmiles(sub_smiles)
    writePDBStream(output, sub_mol)
    pdb_string = output.getvalue()
    rd_mol = AllChem.MolFromPDBBlock(pdb_string)
    new_mol1 = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
    new_mol = Chem.AddHs(new_mol1, addCoords=True, addResidueInfo=True)
    return new_mol


def write_pdb(protein, pdb_name):
    output_pdb_name = f"{pdb_name}_protein.pdb"
    writePDB(f"{output_pdb_name}", protein)
    print(f"wrote {output_pdb_name}")


def write_sdf(new_mol, pdb_name, res_name):
    outfile_name = f"{pdb_name}_{res_name}_ligand.sdf"
    writer = Chem.SDWriter(f"{outfile_name}")
    writer.write(new_mol)
    print(f"wrote {outfile_name}")
    Chem.MolToPDBFile(new_mol, f"{pdb_name}_{res_name}_ligand.pdb")


def main(pdb_name):
    protein, ligand = get_pdb_components(pdb_name)
    write_pdb(protein, pdb_name)

    res_name_list = list(set(ligand.getResnames()))
    for res in res_name_list:
        new_mol = process_ligand(ligand, res)
        write_sdf(new_mol, pdb_name, res)
        cmd3 = f'antechamber -i {pdb_name}_{res}_ligand.sdf -fi sdf -nc 0 -o {pdb_name}_ligand.mol2 -fo mol2 -c bcc -s 1 -pf yes -at gaff -m 1 -rn {res}'
        cmd4 = f'parmchk2 -i {pdb_name}_ligand.mol2 -f mol2 -o {pdb_name}_ligand.frcmod'
        os.system(cmd3)
        os.system(cmd4)
        receptor = f"{pdb_name}_protein.pdb"
        ligand = f"{pdb_name}_{res}_ligand.pdb"
        padding = 10
        cmd5 = f'cat {receptor} {ligand} >comp.pdb'
        cmd51 = 'sed -i "s/END/TER/g" comp.pdb'
        os.system(cmd5)
        os.system(cmd51)
        water_box(receptor, ligand, padding, f'{res}')


def water_box(receptor, ligand, padding, res):
    initial = f'''source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff
loadoff atomic_ions.lib
loadamberparams frcmod.ionsjc_tip3p
loadamberparams {pdb_name}_ligand.frcmod
{res} = loadmol2 {pdb_name}_ligand.mol2
saveoff {res} {pdb_name}_ligand.lib
saveamberparm {res} {pdb_name}_ligand.prmtop {pdb_name}_ligand.rst7

loadoff {pdb_name}_ligand.lib

REC = loadpdb {receptor}
LIG = loadpdb {ligand}
#COMPL = combine{{REC LIG}}
COMPL = loadpdb comp.pdb
#saveAmberParm LIG ligand.prmtop ligand.inpcrd
#saveAmberParm REC receptor.prmtop protein.inpcrd
#saveAmberParm COMPL complex.prmtop complex.inpcrd
solvatebox COMPL TIP3PBOX {padding}
savepdb COMPL solv.pdb
#saveamberparm COMPL solv.prmtop solv.inpcrd
quit'''

    with open('initial.in', 'w') as f:
        f.write(initial)
    cmd = "tleap -f initial.in"
    os.system(cmd)

    with open('countions_for_saltconc1.tcl', 'w') as counteri:
        counteri.write("""set saltConcentration 0.154
mol delete all
mol load pdb solv.pdb 
set sel [atomselect top "water and noh"];
set nWater [$sel num];
$sel delete
if {$nWater == 0} {
    error "ERROR: Cannot add ions to unsolvated system."
    exit
}
set all [ atomselect top all ]
set charge [measure sumweights $all weight charge]
set intcharge [expr round($charge)]
set chargediff [expr $charge - $intcharge]
if { ($chargediff < -0.01) || ($chargediff > 0.01) } {
    error "ERROR: There is a problem with the system. The system does not seem to have integer charge."
    exit
}
puts "System has integer charge: $intcharge"
set cationStoich 1
set anionStoich 1
set cationCharge 1
set anionCharge -1
set num [expr {int(0.5 + 0.0187 * $saltConcentration * $nWater)}]
set nCation [expr {$cationStoich * $num}]
set nAnion [expr {$anionStoich * $num}]
if { $intcharge >= 0 } {
    set tmp [expr abs($intcharge)]
    set nCation [expr $nCation - round($tmp/2.0)]
    set nAnion  [expr $nAnion + round($tmp/2.0)] 
    if {$intcharge%2!=0} {
    set nCation [expr $nCation + 1]}
    puts "System charge is positive, so add $nCation cations and $nAnion anions"
} elseif { $intcharge < 0 } {
    set tmp [expr abs($intcharge)]
    set nCation [expr $nCation + round($tmp/2.0)]
    set nAnion  [expr $nAnion - round($tmp/2.0)]
    if {$intcharge%2!=0} { 
    set nAnion [expr $nAnion + 1]}
    puts "System charge is negative, so add $nCation cations and $nAnion anions"
}
if { [expr $intcharge + $nCation - $nAnion] != 0 } {
    error "ERROR: The calculation has gone wrong. Adding $nCation cations and $nAnion will not result in a neutral system!"
    exit
}
puts "\n";
puts "Your system already has the following charge: $intcharge"
puts "Your system needs the following ions to be added in order to be \
neutralized and have a salt concentration of $saltConcentration M:"
puts "\tCations of charge $cationCharge: $nCation"
puts "\tAnions of charge $anionCharge: $nAnion"
puts "The total charge of the system will be [expr $intcharge + $nCation - $nAnion]."
puts "\n";
exit """)

    cmd6 = 'vmd -dispdev text -eofexit <countions_for_saltconc1.tcl >ions.log'
    os.system(cmd6)

    with open("ions.log",'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'Cations of charge 1' in line:
                cations = str(line.split(':')[1].strip())
            elif 'Anions of charge -1' in line:
                anions = str(line.split(':')[1].strip())
            else:
                pass

    ions_rand = f'addIonsRand COMPL Na+ {cations} Cl- {anions} 5'

    second = f'''source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff
loadoff atomic_ions.lib
loadamberparams frcmod.ionsjc_tip3p
loadamberparams {pdb_name}_ligand.frcmod
{res} = loadmol2 {pdb_name}_ligand.mol2
saveoff {res} {pdb_name}_ligand.lib
saveamberparm {res} {pdb_name}_ligand.prmtop {pdb_name}_ligand.rst7

loadoff {pdb_name}_ligand.lib

REC = loadpdb {receptor}
LIG = loadpdb {ligand}
#COMPL = combine{{REC LIG}}
COMPL = loadpdb comp.pdb
#saveAmberParm LIG ligand.prmtop ligand.inpcrd
#saveAmberParm REC receptor.prmtop protein.inpcrd
saveAmberParm COMPL complex.prmtop complex.inpcrd
solvatebox COMPL TIP3PBOX {padding}
{ions_rand}
savepdb COMPL solv.pdb
saveamberparm COMPL solv.prmtop solv.inpcrd
quit'''

    with open('second.in', 'w') as f:
        f.write(second)
    cmd00 = "tleap -f second.in"
    os.system(cmd00)
    os.system('rm initial.in second.in countions_for_saltconc1.tcl *.log sqm.* ')
    return 


if __name__ == "__main__":
    if len(sys.argv) == 2:
        pdb_name = sys.argv[1]
        main(pdb_name)

    else:
        print("Usage: {sys.argv[1]} pdb_id", file=sys.stderr)

