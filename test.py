import unittest, os, shutil
from simulation import *

class TestSimulationSimple(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestSimulationSimpleDir', 
                              'test_files/acetate.gro', 
                              'ACT', 
                              'test_files/guanidine.gro', 
                              'GUA',
                              'test_files/topology.top',
                              1, 
                              1, 
                              False)

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
        

class TestSimulation(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestSimulationDir', 
                              'test_files/acetate.gro', 
                              'ACT', 
                              'test_files/guanidine.gro', 
                              'GUA',
                              'test_files/topology.top',
                              1, 
                              1, 
                              True)

    def tearDown(self):
        os.chdir('..')        
#        shutil.rmtree('TestSimulationDir')
    
    def test_emin(self):
        self.sim.emin()
        self.assertTrue(os.path.exists('mix.gro'))
        

if __name__ == '__main__':
    unittest.main()
