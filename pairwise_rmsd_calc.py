import os, sys
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')


def pairwise_rmsd(topology, trajectory):
    u = mda.Universe(topology, trajectory)

    aligner = align.AlignTraj(u, u, select='name CA',
                          in_memory=True).run()

    matrix = diffusionmap.DistanceMatrix(u, select='name CA').run()
    plt.imshow(matrix.results.dist_matrix, cmap='viridis')
    plt.xlabel('Frame')
    plt.ylabel('Frame')
    plt.colorbar(label=r'RMSD ($\AA$)')
    plt.savefig('pairwise_rmsd.png')

    return

def pairwise_rmsd_two_trjs(topology1, trajectory1, topology2, trajectory2):
    u1 = mda.Universe(topology1, trajectory1)
    u2 = mda.Universe(topology2, trajectory2)
    prmsd = np.zeros((len(u2.trajectory),  
                  len(u1.trajectory)))

    for i, frame_u2 in enumerate(u2.trajectory):
        r = rms.RMSD(u1, u2, select='name CA',
                 ref_frame=i).run()
        prmsd[i] = r.results.rmsd[:, -1]

    plt.imshow(prmsd, cmap='viridis')
    plt.xlabel('Frame (First trj1)')
    plt.ylabel('Frame (Second trj2)')
    plt.colorbar(label=r'RMSD ($\AA$)')
    plt.savefig('pairwise_rmsd_twotrjs.png')

    return 


if len(sys.argv) == 3:

    topology = sys.argv[1]
    trajectory = sys.argv[2]
    pairwise_rmsd(topology, trajectory)

elif len(sys.argv) == 5: 
    topology1 = sys.argv[1]
    trajectory1 = sys.argv[2]
    topology2 = sys.argv[3]
    trajectory2 = sys.argv[4]

    pairwise_rmsd_two_trjs(topology1, trajectory1, topology2, trajectory2)

else:
    print("Usauage: python3 pairwise_rmsd_calc.py topology1 trajectory1 topology2 trajectory2")
