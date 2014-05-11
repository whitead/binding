import sys
from simulation import *


def run(resname, shortname, anion_atoms, cation_atoms, cset=[1,5,10], pset=[1, 10, 50]):
    for ci in cset:
        for pi in pset:
            sim = Simulation('run_{}_{}_{}'.format(shortname, ci, pi),
                             'structures/acetate.gro', 
                             'ACT', 
                             'structures/{}.gro'.format(resname), 
                             shortname,
                             'structures/topology.top',
                             pressure=pi, 
                             cation_number=ci,
                             cation_atoms=cation_atoms,
                             anion_atoms=anion_atoms)
            sim.emin()
            sim.equilibrate()
            sim.production()
            os.chdir('..')

if(len(sys.argv) != 2):
    print 'must specify simulation'

if(sys.argv[1] == 'GUA'):
    run('guanidine', 'GUA', [1,2], [8])
elif(sys.argv[1] == 'MAM'):    
    run('methylamine', 'MAM', [1,2], [4])
elif(sys.argv[1] == 'BET'):
    run('betaine', 'BET', [1,2], [1])
else:
    print 'invalid argument'
