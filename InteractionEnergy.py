import os, sys
import pandas as pd
import MDAnalysis as mda
import warnings 

warnings.filterwarnings('ignore')

def InteractionEnergy(topology, trajectory, receptor_resid, ligand):
    selection1 = f"protein and resid {receptor_resid}"
    selection2 = f"resname {ligand}"
    outputname = f"IE_{receptor_resid}_{ligand}.dat"

    IE = f"""mol new {topology}
        mol addfile {trajectory} type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        set prot [atomselect top "{selection1}"]
        set ligand [atomselect top "{selection2}"]
        source /usr/local/lib/vmd/plugins/noarch/tcl/namdenergy1.3/namdenergy.tcl
        namdenergy -exe namd3 -elec -vdw -sel $ligand $prot -ofile "{outputname}" -tempname "{outputname}_temp" -ts 10 -timemult 2 -stride 1 -switch  7.5 -cutoff 9 -p {topology}
        quit
    """

    with open("InteractionEnergy.tcl", "w") as w:
        w.write(IE)

    os.system(f"vmd -dispdev text -e InteractionEnergy.tcl")
    os.system("rm InteractionEnergy.tcl")
    tot = pd.read_table(outputname, sep="\s+")["Total"]
    out = list(tot)
    #    os.system(f'rm -r {outputname}')

    return



inpfile = sys.argv[1]
ligname = sys.argv[2]

topology = inpfile.split(".")[0] + ".prmtop"
trajectory = inpfile.split(".")[0] + "_output.xtc"

u = mda.Universe(inpfile)
binding_pocket = u.select_atoms(f'protein and byres around 3.5 resname {ligname}').residues.resids

for resid in list(binding_pocket):
    print(resid)
    InteractionEnergy(f'{topology}', f'{trajectory}', f'{resid}', f'{ligname}')


print("All the calculations are done !!!")
print("Happy Ending")
