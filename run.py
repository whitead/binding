import sys
from simulation import *


def run(resname, shortname, anion_atoms, cation_atoms, index, cset=[1,5,10], pset=[1, 10, 50]):
    i = 0
    box_dim=[3,3,3]
    for ci,bi in zip(cset, box_dim):
        for pi in pset:
            if(i != index):
                i += 1
                continue
            i += 1
            sim = Simulation('run_{}_{}_{}'.format(shortname, ci, pi),
                             'structures/acetate.gro', 
                             'ACT', 
                             'structures/{}.gro'.format(resname), 
                             shortname,
                             'structures/topology.top',
                             pressure=pi, 
                             cation_number=ci,
                             cation_atoms=cation_atoms,
                             anion_atoms=anion_atoms,
                             box_dim=bi)
        
            try:
                sim.emin()
                sim.equilibrate()
                sim.production()
                sim.analyze()
            except RuntimeError:
                print "run_{}_{}_{} failed".format(shortname, ci, pi)                
            finally:
                os.chdir('..')
            

if(len(sys.argv) != 3):
    print 'must specify simulation and index'

index = int(sys.argv[2])
if(sys.argv[1] == 'GUA'):
    run('guanidine', 'GUA', [1,2], [2,6,7], index=index)
elif(sys.argv[1] == 'MAM'):    
    run('methylamine', 'MAM', [1,2], [4], index=index)
elif(sys.argv[1] == 'BET'):
    run('betaine', 'BET', [1,2], [1], index=index)
else:
    print 'invalid argument'
