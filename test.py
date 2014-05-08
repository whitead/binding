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
                              pressure=1, 
                              cation_number=1,
                              cation_atoms=[8],
                              anion_atoms=[1,2])

    def tearDown(self):
        os.chdir('..')        
        #shutil.rmtree('TestSimulationDir')
    
    def test_all(self):
        self.sim.emin()
        self.assertTrue(os.path.exists('emin.gro'))
        self.sim.equilibrate()
        self.assertTrue(os.path.exists('equil.gro'))
        self.sim.production()
        self.assertTrue(os.path.exists('prod.gro'))
        self.analyze()
        self.assertTrue(os.path.exists('analysis/pmf.dat'))
        self.assertTrue(os.path.exists('analysis/pmf.png'))


class TestMethylamine(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestMethylamineDir', 
                              'structures/acetate.gro', 
                              'ACT', 
                              'structures/methylamine.gro', 
                              'MAM',
                              'structures/topology.top',
                              pressure=1, 
                              cation_number=1,
                              cation_atoms=[4],
                              anion_atoms=[1,2])

    def tearDown(self):
        os.chdir('..')        
        #shutil.rmtree('TestMethylamineDir')

    def test_all(self):
        self.sim.emin()
        self.assertTrue(os.path.exists('emin.gro'))
        self.sim.equilibrate()
        self.assertTrue(os.path.exists('equil.gro'))
        self.sim.production()
        self.assertTrue(os.path.exists('prod.gro'))
        self.analyze()
        self.assertTrue(os.path.exists('analysis/pmf.dat'))
        self.assertTrue(os.path.exists('analysis/pmf.png'))


class TestBetaine(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation('TestBetaineDir', 
                              'structures/acetate.gro', 
                              'ACT', 
                              'structures/betaine.gro', 
                              'BET',
                              'structures/topology.top',
                              pressure=1, 
                              cation_number=1,
                              cation_atoms=[1],
                              anion_atoms=[1,2])

    def tearDown(self):
        os.chdir('..')        
        #shutil.rmtree('TestBetaineDir')



    def test_prod(self):
        self.sim.emin()
        self.assertTrue(os.path.exists('emin.gro'))
        self.sim.equilibrate()
        self.assertTrue(os.path.exists('equil.gro'))
        self.sim.production()
        self.assertTrue(os.path.exists('prod.gro'))
        self.analyze()
        self.assertTrue(os.path.exists('analysis/pmf.dat'))
        self.assertTrue(os.path.exists('analysis/pmf.png'))    


if __name__ == '__main__':
    unittest.main()
