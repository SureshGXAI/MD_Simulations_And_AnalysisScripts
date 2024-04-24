import os, sys
import parmed as pmd

inpfile = sys.argv[1]

parm = pmd.load_file(f'{inpfile}.prmtop', f'{inpfile}.inpcrd')
parm.save('gromacs.top', format='gromacs')
parm.save('gromacs.gro')

parm.save('charmm.psf')
parm.save('charmm.crd')
