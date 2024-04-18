##################################################################################
# Analysis of the MD trajectories obtained from MD simulations
# First version: Suresh Gorle, 16 April 2024
# It does trajectory wrapping, writing complex only trajectory, 
# Tools required: pandas, seaborn, scipy, matplotlib, parmed. MDAnalysis, VMD, NAMD
#
##################################################################################

import os, sys
import MDAnalysis as mda
import warnings
warnings.filterwarnings('ignore')
from matplotlib import pyplot as plt
from MDAnalysis.analysis import rms, align, pca, contacts
from MDAnalysis import transformations
import parmed as pmd
import pandas as pd
import seaborn as sns
from scipy.stats import norm


def wrap_dcd(topology, trajectory, wrap_trj):
        wrap_protocol = f'''mol new {topology}
mol addfile {trajectory} type xtc first 1 last -1 step 1 filebonds 1 autobonds 1 waitfor all
package require pbctools
pbc wrap -all -compound res -center bb -centersel "protein or nucleic"
animate write dcd {wrap_trj}
quit'''

        with open('wrap.tcl', 'w') as f:
            f.write(wrap_protocol)

        os.system(f'vmd -dispdev text -e wrap.tcl')
        os.system("rm wrap.tcl")
        return 

def convert_dcd_to_xtc(pdb, dcd):
    u = mda.Universe(pdb, dcd)
    protein = u.select_atoms("protein or resname RES")
    system = u.select_atoms("all")
    with mda.Writer("wrapped_trajectory.xtc", system.n_atoms) as w:
        for ts in u.trajectory:
            w.write(system)
    return 

def RemoveWaterIons(pdb, xtc):
    u = mda.Universe(pdb, xtc)
    protein = u.select_atoms("protein or resname RES")
    system = u.select_atoms("all")
    with mda.Writer("trimmed_trajectory.xtc", protein.n_atoms) as w:
        for ts in u.trajectory:
            w.write(protein)
    return

def split_complex(pdb):
    u = mda.Universe(pdb)
    pl = u.select_atoms("protein or resname RES")
    p = u.select_atoms("protein")
    l = u.select_atoms("resname RES")

    mda.Writer("complex_only.pdb", pl.n_atoms).write(pl)
    mda.Writer("protein_only.pdb", p.n_atoms).write(p)
    mda.Writer("ligand_only.pdb", l.n_atoms).write(l)
    return 


def write_indivdual_frames(topology, xtc):
    if not os.path.exists('for_frames'):
        os.mkdir('for_frames')

    u = mda.Universe(topology, xtc)
    for i, ts in enumerate(u.trajectory):
        u_protein = u.select_atoms('protein')
        proteinfile = f'for_frames/protein_{i}.pdb'

        u_ligand = u.select_atoms('resname RES')
        ligandfile = f'for_frames/ligand_{i}.pdb'

        with mda.Writer(proteinfile, u_protein.n_atoms) as w:
            w.write(u_protein)

        with mda.Writer(ligandfile, u_ligand.n_atoms) as w:
            w.write(u_ligand)
    return 



def rmsd(topology, xtc):
    trj = mda.Universe(topology, xtc)
    trj.trajectory[0] # set to first frame
    rmsd_analysis = rms.RMSD(trj, select='backbone', groupselections=['name CA', 'protein'])
    rmsd_analysis.run()

    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:], columns=['Backbone', 'C-alphas', 'Protein'], index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df.index.name = 'Time (ps)'
    rmsd_df.plot()
    plt.savefig("rmsd.png")
    return 

def radius_of_gyration(topology, xtc):
    trj = mda.Universe(topology, xtc)
    time = []
    rgyr = []
    for ts in trj.trajectory:
        time.append(trj.trajectory.time)
        rgyr.append(trj.select_atoms("protein").radius_of_gyration())

    rgyr_df = pd.DataFrame(rgyr, columns=['Radius of gyration (A)'], index=time)
    rgyr_df.index.name = 'Time (ps)'

    rgyr_df.plot(title='Radius of gyration')
    plt.savefig("rgyr.png")
    return 

def RMSF(topology, xtc):
    trj = mda.Universe(topology, xtc)
    average = align.AverageStructure(trj, trj, select='protein and name CA', ref_frame=0).run()
    ref = average.results.universe
    aligner = align.AlignTraj(trj, ref, select='protein and name CA', in_memory=True).run()
    c_alphas = trj.select_atoms('protein and name CA')
    R = rms.RMSF(c_alphas).run()
    plt.plot(c_alphas.resids, R.results.rmsf)
    plt.xlabel('Residue number')
    plt.ylabel('RMSF ($\AA$)')
    plt.savefig("rmsf.png")
    return 

def Native_contacts(topology, xtc):
    trj = mda.Universe(topology, xtc)
    receptor = trj.select_atoms('protein')
    ligand = trj.select_atoms('resname RES')

    ca1 = contacts.Contacts(trj,
                        select=("protein", "resname RES"),
                        refgroup=(receptor, ligand),
                        radius=4.5,
                        method='hard_cut').run()

    ca1_df = pd.DataFrame(ca1.results.timeseries,
                      columns=['Frame',
                               'Contacts from first frame'])
    ca1_df.plot(x='Frame')
    plt.ylabel('Fraction of contacts')
    plt.savefig("fraction_contacts.png")
    return 

def PCA(topology, xtc):
    trj = mda.Universe(topology, xtc)
    aligner = align.AlignTraj(trj, trj, select='backbone', in_memory=True).run()
    pc = pca.PCA(trj, select='backbone', align=True, mean=None, n_components=None).run()
    plt.plot(pc.cumulated_variance[:10])
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative variance')
    plt.savefig("pcavariance.png")

    transformed = pc.transform(trj.select_atoms('backbone'), n_components=3)
    df = pd.DataFrame(transformed, columns=['PC{}'.format(i+1) for i in range(3)])
    df['Time (ps)'] = df.index * trj.trajectory.dt

    g = sns.PairGrid(df, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df)))
    g.map(plt.scatter, marker='.')
    plt.savefig("three_pca_components.png")
    return 


def hbond_occupancy(topology, xtc, typ):

    hbond=f'''set trj {xtc}
set prmtop {topology}
mol new $trj type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $prmtop type {typ} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set top [atomselect top "all"]
set prot [atomselect top "protein"]
set lig [atomselect top "resname RES"]

source /usr/local/lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl
hbonds -sel1 $prot -sel2 $lig -type unique -writefile yes -dist 3.5 -ang 120 -polar yes -detailout hbond-details_occupancy.dat
exit
'''

    with open('occupancy.tcl', 'w') as w:
        w.write(hbond)

    os.system(f'vmd -dispdev text -e occupancy.tcl')
    os.system("rm occupancy.tcl")
    return 


def InteractionEnergy(topology, xtc, receptor_resid, ligand):
    selection1 = f'protein and resid {receptor_resid}'
    selection2 = f'resname {ligand}'
    outputname = f'IE_{receptor_resid}_{ligand}.dat'

    IE = f'''mol new {topology}
        mol addfile {xtc} type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        set prot [atomselect top "{selection1}"]
        set ligand [atomselect top "{selection2}"]
        source /usr/local/lib/vmd/plugins/noarch/tcl/namdenergy1.3/namdenergy.tcl
        namdenergy -exe namd3 -elec -vdw -sel $ligand $prot -ofile "{outputname}" -tempname "{outputname}_temp" -ts 10 -timemult 2 -stride 1 -switch  7.5 -cutoff 9 -p {topology} 
        quit
    '''

    with open('InteractionEnergy.tcl', 'w') as w:
        w.write(IE)

    os.system(f'vmd -dispdev text -e InteractionEnergy.tcl')
    os.system("rm InteractionEnergy.tcl")
    tot = pd.read_table(outputname, sep='\s+')['Total']
    out = list(tot)
    os.system(f'rm -r {outputname}')

    return 


def NormalDistribution(inpfile):

    fig = plt.subplots(figsize=(12, 6))

    hbond = pd.read_csv(f'{inpfile}', sep='[ ]', header=None, names=["Time", "Distance"])
    data = norm.rvs(hbond['Distance'])
    ax = sns.distplot(data, bins=50, kde=True, color='green', hist_kws={"linewidth": 15,'alpha':1})

    plt.xlabel('Hbond Distance ($\AA$)', fontsize=15)
    plt.ylabel('Frequency', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(0, 6)
    plt.ylim(0, 1)
    plt.savefig('hbond.jpg')
    return 

def time_series_plot(inpfile):
    rmsd = pd.read_csv(f'{inpfile}', sep='[ ]', header=None, names=["Time", "RMSD"])

    rmsd['NewTime'] = rmsd['Time'] * 0.002
    plt.plot(rmsd['NewTime'], rmsd['RMSD'], color='blue', lw=2.0)
    plt.xlabel('Time (ns)', fontsize=15)
    plt.ylabel('RMSD ($\AA$)', fontsize=15)
    plt.xlim(0, 50)
    plt.ylim(0, 5)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('rmsd.jpg')
    return

def extract_binding_pocket(inpfile):
    pocket = f'''mol load pdb {inpfile}.pdb
set complex [atomselect top "(protein and same residue as within 5 of resname RES)"]
set ligand [atomselect top "resname RES"]
$complex writepdb {inpfile}_pocket.pdb
$ligand writepdb {inpfile}_ligand.pdb'''

    with open("vmd_pocket_extraction.tcl", 'w') as w:
        w.write(pocket)
    os.system(f'vmd -dispdev text -e vmd_pocket_extraction.tcl')
    os.system("rm vmd_pocket_extraction.tcl")
    return 


def main(topology, trajectory, wrap_traj, inpfile):
    #Removing the pbc effects 
#    wrap_dcd(topology, trajectory, wrap_traj)

    #Converting from dcd to xtc format
#    convert_dcd_to_xtc(str(inpfile)+".pdb", wrap_traj)
#    os.system("rm wrapped.dcd")

    #Removing water and ions from the xtc file 
    ipdb = str(inpfile)+".pdb"
    dcd = "wrapped_trajectory.xtc"
    RemoveWaterIons(ipdb, dcd)

    #Seperating complex, protein and ligand from initial structure
    split_complex(ipdb)

    xtc = "trimmed_trajectory.xtc"
    top = "complex_only.pdb"
    #Calculating the dynamics of protein
    RMSF(top, xtc)
    rmsd(top, xtc)
    radius_of_gyration(top, xtc)
 #   PCA(top, xtc)

    #calculating interactions between protein and ligand
    Native_contacts(top, xtc)
    hbond_occupancy(top, xtc, "pdb")
    return 


if __name__ == "__main__":
    if len(sys.argv) > 2:
        inpfile = sys.argv[1] 
        topology = str(inpfile)+".prmtop"
        trajectory = str(sys.argv[2])+".xtc"
        wrap_traj = "wrapped.dcd"
        main(topology, trajectory, wrap_traj, inpfile)
    else:
        print("Usage: python3 recenter.py {topologyfile} {trajectoryfile}")
