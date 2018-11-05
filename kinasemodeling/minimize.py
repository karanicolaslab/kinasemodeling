from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmoltools.forcefield_generators import gaffTemplateGenerator
from glob import glob
import os

# Written by G. A. Ross, S. K. Albanese, J. Gao, K. Palani

#------Functions used to clean and minimizing with openmm------#
def get_upper_location(pdb_location, back=2):
    """
    Get the path of the directory that's up from the path specified from pdb_location.
    """
    path = pdb_location.split('/')[0:-back]
    upper_path = ''
    for location in path:
        upper_path += location + '/'
    return upper_path


def openmm_minimize(pdb_filename, pdbname, out_folder='minimized', gpu=False, solvate=False):
    """
    Minimize a supplied system with openmm. Can handle small molecules as long as CONECT
    records are supplied.
    """
    # Initialize forcefield with small molecule capabilities
    forcefield = ForceField('data/gaff.xml','tip3p.xml','amber99sbildn.xml')
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)

    # Use modeller to remove unwanted residues
    pdb = PDBFile(pdb_filename)
    model = Modeller(pdb.topology, pdb.positions)

    # Remove unwanted molecules

    # Add waters in a cubic box
    if solvate == True:
        model.addSolvent(forcefield, padding=1.0*nanometers)

    # Create the system with a cheap electrostatic cutoff
    system = forcefield.createSystem(model.topology, nonbondedMethod=CutoffNonPeriodic)

    # Minimize system with a placeholder integrator
    integrator = VerletIntegrator(0.001*picoseconds)
    if gpu == True:
        platform = Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed'}
        simulation = Simulation(model.topology, system, integrator, platform, properties)
    else:
        simulation = Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)
    simulation.minimizeEnergy()

    # Print PDB
    positions = simulation.context.getState(getPositions=True).getPositions()
    out_directory = get_upper_location(pdb_filename)
    try:
        os.mkdir(out_directory + out_folder)
    except OSError:
        pass
    PDBFile.writeFile(simulation.topology, positions, open(out_directory + out_folder + '/' + pdbname + '-minimized.pdb', 'w'))

if __name__ == "__main__":

    import logging
    logging.basicConfig(filename='explicit_solvent_minimization.log',level=logging.DEBUG)
    logging.captureWarnings(True)

    #------Cleaning and minimizing with openmm------#

    # Running through refined structures and seeing if they work
    filelist = ['./data/pdb/minimize_test.pdb']
    filename= './data/pdb/minimize_test.pdb'
    pdb = 'minimize_test.pdb'
    for file in filelist:
        logging.info('Refining structure {0}'.format(pdb))
        try:
            openmm_minimize(filename, pdb, out_folder='minimized', gpu=False, solvate=False)
            logging.info('Structure {0} COMPLETED.'.format(pdb))
            logging.info('Original file located in {0}'.format(pdb, filename))
            logging.info('-----------------------------')
        except Exception as detail:
            logging.error(detail)
            logging.info('Structure {0} BROKEN'.format(pdb))
            logging.info('Original file located in {0}'.format(pdb, filename))
            logging.info('-----------------------------')
