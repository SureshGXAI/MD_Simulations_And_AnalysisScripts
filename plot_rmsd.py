import os, sys
import pandas as pd
import matplotlib.pyplot as plt

def time_series_plot(inpfile1, inpfile2):
    rmsd = pd.read_csv(f'{inpfile1}', sep='\s+', header=None, names=["Time", "RMSD"])
    rmsd1 = pd.read_csv(f'{inpfile2}', sep='\s+', header=None, names=["Time", "RMSD"])

    rmsd['NewTime'] = rmsd['Time'] * 0.002
    rmsd1['NewTime'] = rmsd1['Time'] * 0.002
    plt.rcParams["figure.autolayout"] = True
    plt.plot(rmsd['NewTime'], rmsd['RMSD'], color='blue', lw=1.0, label="Protein")
    plt.plot(rmsd1['NewTime'], rmsd1['RMSD'], color='orange', lw=1.0, label="Ligand")

    plt.legend(loc='upper center')
    plt.xlabel('Time (ns)', fontsize=15)
    plt.ylabel('RMSD ($\AA$)', fontsize=15)
    plt.xlim(0, 50)
    plt.ylim(0, 4)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('rmsd.png')
    return



inpfile1 = sys.argv[1]
inpfile2 = sys.argv[2]

time_series_plot(inpfile1, inpfile2)
