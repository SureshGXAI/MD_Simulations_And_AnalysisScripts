import os, sys
import pandas as pd 
from matplotlib import pyplot as plt
import numpy as np


def rmsf_calc(topology, trajectory, ligname):
    rmsf = """set top [lindex $argv 0]
set trj [lindex $argv 1]
set ligname [lindex $argv 2]
mol new $top
mol addfile $trj type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

package require pbctools
pbc wrap -all -compound res -center bb -centersel "protein or nucleic"

set frame0 [atomselect top "protein and backbone and noh" frame 0]
set lig [atomselect top "resname $ligname and noh"]
set prot [atomselect top "protein and backbone and name CA"]

set nf [molinfo top get numframes]

for {set i 1 } {$i < $nf } { incr i } {
    set sel [atomselect top "protein and backbone and noh" frame $i]
    set all [atomselect top all frame $i]
    $all move [measure fit $sel $frame0]
}

puts [open rmsf_protein.dat w] "[measure rmsf $prot first 1 last -1 step 100]"
puts [open rmsf_ligand.dat w] "[measure rmsf $lig first 1 last -1 step 100]"
puts [open ligand_names.txt w] "[$lig get name]"

quit"""
    with open("rmsf_calc.tcl", 'w') as w:
        w.write(rmsf)
    os.system(f'vmd -dispdev text -e rmsf_calc.tcl -args {topology} {trajectory} {ligname}')
    os.system("rm rmsf_calc.tcl")

    return 



def write_file(inpfile, outfile, gap):
    out = open(outfile, 'w')
    with open(inpfile, 'r') as f1:
        for lines in f1.readlines():
            words = lines.strip().split()
            for k in range(len(words)):
                out.write(str(int(k)+int(gap)+1) +' '+ str(words[k])+ '\n')
    f1.close()
    out.close()

    return 

def plot_for_receptor(inpfile1):
    rmsd = pd.read_csv(f'{inpfile1}', sep='\s+', header=None, names=["Time", "RMSF"])

    plt.rcParams["figure.figsize"] = [12.00, 4.50]
    plt.rcParams["figure.autolayout"] = True
    plt.plot(rmsd['Time'], rmsd['RMSF'], color='blue', lw=1.0, label="Receptor", marker=".")

    plt.legend(loc='upper center')
    plt.xlabel('Residue Number', fontsize=15)
    plt.ylabel('RMSF ($\AA$)', fontsize=15)
    plt.xlim(rmsd['Time'].min(), rmsd['Time'].max())
    plt.ylim(0, rmsd['RMSF'].max())
    plt.xticks(np.arange(rmsd['Time'].min(), rmsd['Time'].max()+1, 50.0))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(str(inpfile1.split('.')[0])+"_rmsf.png")
    return


def plot_for_ligand(inpfile1):
    rmsd = pd.read_csv(f'{inpfile1}', sep='\s+', header=None, names=["Time", "RMSF"])

    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    plt.plot(rmsd['Time'], rmsd['RMSF'], color='blue', lw=1.0, label="Ligand", marker="o")

    plt.legend(loc='upper center')
    plt.xlabel('Atom Number', fontsize=15)
    plt.ylabel('RMSF ($\AA$)', fontsize=15)
    plt.xlim(0, rmsd['Time'].max()+1)
    plt.ylim(0, rmsd['RMSF'].max()+1)
    plt.xticks(np.arange(0, rmsd['Time'].max()+1, 2.0))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.savefig(str(inpfile1.split('.')[0])+"_rmsf.png")
    return


def ligand_2d_image_generation(pdbfile):
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdmolfiles as rdmd
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = rdmd.MolFromPDBFile(f'{pdbfile}')
    tmp=AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(1500, 1500)
    d.drawOptions().addAtomIndices = True
    rdMolDraw2D.PrepareAndDrawMolecule(d,mol)
    d.FinishDrawing()
    d.WriteDrawingText('ligand_2d.png')
    return 


def plotting(inpfile1, gap):
    if 'ligand' in inpfile1:
        write_file(inpfile1, 'out_'+inpfile1, gap)
        ligand_2d_image_generation(ligandfile)
        plot_for_ligand('out_'+inpfile1)
    elif 'protein' in inpfile1:
        write_file(inpfile1, 'out_'+inpfile1, gap)
        plot_for_receptor('out_'+inpfile1)
    else:
        pass

    return 



def main():
    rmsf_calc(topology, trajectory, ligname)
    plotting("rmsf_ligand.dat", -1)
    plotting("rmsf_protein.dat", pgap)

    return 


if __name__ == "__main__":
    if len(sys.argv) > 5:
        topology = sys.argv[1]
        trajectory = sys.argv[2]
        ligname = sys.argv[3]
        pgap = sys.argv[4]
        ligandfile = sys.argv[5]

        main()
    else:
        print("Usage: python3 plot_rmsf.py topology trajectory ligname pgap ligandfile")
