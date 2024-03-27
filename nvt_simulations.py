from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from copy import deepcopy
import os, sys

inpfile = sys.argv[1] 
temperature = TEMP
Timesteps = 250000
workDir = f'Dir_{temperature}'

os.mkdir(f'Dir_{temperature}')

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

# Load checkpoint
# reset step and time counters
simulation.loadCheckpoint(('eq.chk'))
eq_state = simulation.context.getState(getVelocities=True, getPositions=True)
positions = eq_state.getPositions()
velocities = eq_state.getVelocities()

simulation.context.setVelocitiesToTemperature(300*kelvin)
integrator.setTemperature(TEMP)
simulation.context.setPositions(positions)

mdsteps = 250000

# append reporters
simulation.reporters.append(
    XTCReporter((f'{workDir}/10ns_traj.xtc'), 1000))
simulation.reporters.append(
    StateDataReporter(
         stdout,
         10000,
         step=True,
         potentialEnergy=True,
         temperature=True,
         progress=True,
         remainingTime=True,
        speed=True,
        totalSteps=Timesteps,
        separator='\t'
        )
    )

# run production simulation fpr 10ns
print('Running Production...')
simulation.step(Timesteps)
print('Done!')
simulation.saveState((f'{workDir}/sim_10ns.state'))
simulation.saveCheckpoint((f'{workDir}/sim_10ns.chk'))
