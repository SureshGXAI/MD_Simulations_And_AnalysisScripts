import os, sys
import pandas as pd
import MDAnalysis as mda
import warnings
warnings.filterwarnings('ignore')


def InteractionEnergy(topology, xtc, ligname):
    selection1 = f"protein or resname {ligname}"
    outputname = f"IE_{ligname}.dat"

    IE = f"""mol new {topology}
        mol addfile {xtc} type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        set complex [atomselect top "{selection1}"]
        source /usr/local/lib/vmd/plugins/noarch/tcl/namdenergy1.3/namdenergy.tcl
        namdenergy -exe namd3 -elec -vdw -sel $complex -ofile "{outputname}" -tempname "{outputname}_temp" -ts 10 -timemult 2 -stride 1 -switch  7.5 -cutoff 9 -p {topology}
        quit
    """

    with open("InteractionEnergy.tcl", "w") as w:
        w.write(IE)

    os.system(f"vmd -dispdev text -e InteractionEnergy.tcl")
    os.system("rm InteractionEnergy.tcl")
    tot = pd.read_table(outputname, sep="\s+")["Total"]
    frame = pd.read_table(outputname, sep="\s+")["Frame"]
#    x, y = (list(t) for t in zip(*sorted(zip(x, y))))
    tot, frame = (list(t) for t in zip(*sorted(zip(tot, frame))))
    framenum = frame[0]

    write_indivdual_frames(topology, xtc, framenum, ligname)
    return 


def write_indivdual_frames(topology, xtc, frame, ligname):

    u = mda.Universe(topology, xtc)
    for i, ts in enumerate(u.trajectory):
        if int(i) == int(frame):
            u_protein = u.select_atoms(f'protein or resname {ligname}')
            proteinfile = f'lowest_energy_complex_{i}.pdb'

            with mda.Writer(proteinfile, u_protein.n_atoms) as w:
                w.write(u_protein)

    return

def main(topology, trajectory, ligname):
    InteractionEnergy(topology, trajectory, ligname)
    return 

if __name__ == "__main__":

    if len(sys.argv) > 3:

        topology = sys.argv[1]
        trajectory = sys.argv[2]
        ligname = sys.argv[3]
        main(topology, trajectory, ligname)
    else:
        print("usuage: python3 lowest_energy_config_pick.py topology trajectory ligandname")
