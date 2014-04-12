import unittest, os, shutil
from simulation import *

class TestSimulationSimple(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestSimulationDir', 'test_files/acetate.pdb', 'test_files/guanidine.pdb', 1, 1, False)

    def tearDown(self):
        os.chdir('..')
        shutil.rmtree('TestSimulationDir')
    
    def test_pack(self):
        self.sim._pack('acetate.pdb', 'guanidine.pdb', 1, 3)
        self.assertTrue(os.path.exists('dry_mix.pdb'))

    def test_solvate(self):
        self.sim._pack('acetate.pdb', 'guanidine.pdb', 1, 3)
        self.sim._solvate(3)
        self.assertTrue(os.path.exists('dry_mix.gro'))
        self.assertTrue(os.path.exists('mix.gro'))
        

if __name__ == '__main__':
    unittest.main()
