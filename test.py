import unittest, os, shutil
from simulation import *

class TestSimulationSimple(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestSimulationSimpleDir', 
                              'structures/acetate.gro', 
                              'ACT', 
                              'structures/guanidine.gro', 
                              'GUA',
                              'structures/topology.top',
                              1, 
                              1, 
                              [1],
                              [1],
                              [1,2,3],
                              initialize=False)

    def tearDown(self):
        os.chdir('..')
        shutil.rmtree('TestSimulationSimpleDir')
    
    def test_pack(self):
        self.sim._pack('acetate.gro', 'ACT', 'guanidine.gro', 'GUA', 1, 3)
        self.assertTrue(os.path.exists('dry_mix.gro'))

    def test_solvate(self):
        self.sim._pack('acetate.gro', 'ACT', 'guanidine.gro', 'GUA', 1, 3)
        self.sim._solvate(3)
        self.assertTrue(os.path.exists('dry_mix.gro'))
        self.assertTrue(os.path.exists('mix.gro'))

    def test_different_numbers(self):
        self.sim._pack('acetate.gro', 'ACT', 'guanidine.gro', 'GUA', 5, 3)
        self.assertTrue(os.path.exists('dry_mix.gro'))


        

class TestSimulation(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestSimulationDir', 
                              'structures/acetate.gro', 
                              'ACT', 
                              'structures/guanidine.gro', 
                              'GUA',
                              'structures/topology.top',
                              1, 
                              1,
                              [1],
                              [1],
                              [1,2,3])

    def tearDown(self):
        os.chdir('..')        
#        shutil.rmtree('TestSimulationDir')
    
    def test_emin(self):
        self.sim.emin()
        self.assertTrue(os.path.exists('emin.gro'))
        
    def test_ionc(self):
        self.sim.ionc = 0.1
        self.sim.emin()
        self.assertTrue(os.path.exists('emin.gro'))

    def test_equil(self):
        self.sim.emin()
        self.sim.equilibrate()
        self.assertTrue(os.path.exists('equil.gro'))

    def test_equil(self):
        self.sim.emin()
        self.sim.equilibrate()
        self.sim.production()
        self.assertTrue(os.path.exists('prod.gro'))

        

if __name__ == '__main__':
    unittest.main()
