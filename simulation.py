import logging, os, shutil, datetime, subprocess


ION_CONC=100#mM
BOX_DIM=3#nm
PACKMOL='$HOME/packmol/packmol'
FF='oplsaa'
WATER='tip4p'
MDRUN='mdrun'

        
class Simulation:
    """
    A simulation that has a directory and is able to be executed
    
    Arguments:
    name -- the name to be used for prefixes
    anion_filename -- a pdb file with anion structure
    cation_filename -- a pdb file with cation structure
    pressure -- the pressure in atmospheres
    cation_number -- the number of cations to use in the simulation
    """
    def __init__(self, name, anion_filename, cation_filename, pressure, cation_number, initialize=True):
        """
        Creates the directory, puts the cation and anion together and solvates
        """
        self.name = name        
        self.current_structure, self.current_top, self.current_run = None, None, None
        self._setup_directory(anion_filename, cation_filename)
        if(initialize):
            self._pack(anion_filename, cation_filname, cation_number, BOX_DIM)
            self._solvate(BOX_DIM)

    def _putInDir(dirname):
        """
        A decorator that encloses a function inside a directory
        """
        def wrap(f):
            def mod_f(self, *args):
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                    #bring files
                for f in [self.current_structure, self.current_top, self.current_run]:
                    if(f is not None):
                        shutil.copyfile(f, os.join(dirname, os.path.basename(f)))            
                #go there
                os.chdir(dirname)
                try:            
                    f(self, *args)
                finally:
                    #make sure we leave
                    os.chdir(self.dir)
                    #bring back files
                    for f in [self.current_structure, self.current_top, self.current_run]:
                        if(f is not None):
                            shutil.copyfile(os.join(dirname, f),f)
            return mod_f
        return wrap
                

    @_putInDir('emin')
    def emin(self):
        """
        Energy minimize the system.        
        """
        #build input file
        input_file = 'emin.mdp'
        with open(input_file, 'w') as f:
            f.write('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = cg
            nsteps                   = 10000; max steps
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 0
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            
            ; NEIGHBORSEARCHING PARAMETERS = 
            rlist                    = 1.3
            cutoff-scheme            = Verlet
            
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulomb_type             = PME-Switch
            rcoulomb_switch          = 1.0
            rcoulomb                 = 1.0
            ; Method for doing Van der Waals = 
            vdw_type                 = Shift
            ; cut-off lengths        = 
            rvdw                     = 1
            dispcorr                 = ener
            ewald_rtol               = 1e-5

            '''.format())

        #prepare run
        self._exec_log('grompp', {'c':self.current_structure, 'f':input_file, 'p':self.current_top})
        #run
        emin_structure = 'emin_noion.gro'
        self._exec_log(MDRUN, {'c':emin_structure, 'f':input_file, 'p':self.current_top})
        self.current_structure = emin_structure
        self.current_run = 'topol.tpr'
        
    def _setup_directory(self, *to_copy):
        #build directory, start log
        if not os.path.exists(self.name):
            os.mkdir(self.name)
        self.dir = os.path.abspath(self.name)
        for f in to_copy:
            shutil.copyfile(f, os.path.join(self.dir, os.path.basename(f)))
        os.chdir(self.dir)
        if os.path.exists('simulation.log'):
            logging.basicConfig(filename='simulation.log',level=logging.DEBUG)
            logging.info('Restarting simulation on {}'.format(datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')))
        else:
            logging.basicConfig(filename='simulation.log',level=logging.DEBUG)
            logging.info('Beginning simulation on {}'.format(datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')))
            

    def _pack(self, afile, cfile, cnumber, box):
        #build input file
        with open('pack.inp', 'w') as f:
            f.write('''
            tolerance 2.0
            filetype pdb
            output dry_mix.pdb
            
            structure {afile}
              number 1 
              inside box 0. 0. 0. {box} {box} {box}
            end structure
            
            structure {cfile}
              number {cnumber}
              inside box 0. 0. 0. {box} {box} {box}
            end structure
            '''.format(afile=afile, cfile=cfile, cnumber=cnumber, box=box*10))
            
        #pack
        self._exec_log('{} < pack.inp'.format(PACKMOL))
        self.current_structure = 'dry_mix.pdb'

    def _exec_log(self, string, arg_dic=None):
        """
        Smartly executes the given string and logs the results
        """
        if(arg_dic is not None):
            for k,v in arg_dic.iteritems():
                string += ' -{} {}'.format(k,v)
        try:
            process = subprocess.Popen(string, shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
                           
            out, err = process.communicate()
            retcode = process.returncode
            logging.debug(out)
            if(len(err) > 0):
                logging.warning(err)
            if retcode < 0:
                logging.error('{} failed with {}. stderr follow'.format(string, retcode))
                logging.error(err)
            else:
                logging.info('{} with arguments'.format(string))
        except OSError as e:
            print >>sys.stderr, 'Execution of {} failed:'.format(string), e

    def _solvate(self, box):
        args = {'ff':FF}
        args['water'] = WATER
        args['f'] = self.current_structure
        args['o'] = 'dry_mix.gro'
        self._exec_log('pdb2gmx', args)
        args = {'cp':'dry_mix.gro'}
        args['o'] = 'mix.gro'
        args['box'] = '{box} {box} {box}'.format(box=box)
        self._exec_log('genbox', args)
        self.current_structure = 'mix.gro'
        self.current_top = 'topol.top'
