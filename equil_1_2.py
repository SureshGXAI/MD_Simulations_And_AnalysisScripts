from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from copy import deepcopy
import os, sys

inpfile = sys.argv[1]
Timesteps = 250000
workDir = "."

# Function to add backbone position restraints
def add_backbone_posres(system, positions, atoms, restraint_force):
  force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * kilocalories_per_mole/angstroms**2
  force.addGlobalParameter("k", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(positions, atoms)):
    if atom.name in  ('CA', 'C', 'N'):
      force.addParticle(i, atom_crd.value_in_unit(nanometers))
  posres_sys = deepcopy(system)
  posres_sys.addForce(force)
  return posres_sys

# Set up 
prmtop = AmberPrmtopFile(inpfile + '.prmtop')
inpcrd = AmberInpcrdFile(inpfile + '.inpcrd')
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}

#System preperation 
system = prmtop.createSystem(nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometers, constraints=HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)
posres_sys = add_backbone_posres(system, inpcrd.positions, prmtop.topology.atoms(), 100)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
integrator.setConstraintTolerance(0.00001)
simulation = Simulation(prmtop.topology, posres_sys, integrator, platform)
simulation.reporters.append(
    StateDataReporter(
        stdout,
        10000,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)

#Minimize
print('Minimizing...')
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.step(3500)

#WarmUp with a NVT run.  Slowly warm up temperature - every 1000 steps raise the temperature by 5 K, ending at 300 K
simulation.context.setVelocitiesToTemperature(5*kelvin)
print('Warming up the system...')
T = 5
mdsteps = 50000
for i in range(60):
  simulation.step(int(mdsteps/60) )
  temperature = (T+(i*T))*kelvin 
  integrator.setTemperature(temperature)


#NPT equilibration, reducing backbone constraints
mdsteps = 250000
barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
simulation.context.reinitialize(True)
print('Running NPT equilibration...')
for i in range(100):
  simulation.step(int(mdsteps/100))
  simulation.context.setParameter('k', (float(99.02-(i*0.98))*kilocalories_per_mole/angstroms**2))
simulation.context.setParameter('k', 0)

# save the equilibration results to file : state is platform independent but less precise, checkpoint file
simulation.saveState(('eq.state'))
simulation.saveCheckpoint(('eq.chk'))
