import logging, os, shutil, datetime, subprocess, re, textwrap


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
    anion_filename -- a gro file with anion structure
    cation_filename -- a gro file with cation structure
    topology_filename -- a top file with itp include statements
    pressure -- the pressure in atmospheres
    cation_number -- the number of cations to use in the simulation
    """
    def __init__(self, 
                 name, 
                 anion_filename, 
                 anion_resname, 
                 cation_filename, 
                 cation_resname, 
                 topology_filename, 
                 pressure, 
                 cation_number, 
                 initialize=True):
        """
        Creates the directory, puts the cation and anion together and solvates
        """
        self.name = name        
        self.current_top = os.path.basename(topology_filename)
        self.current_structure, self.current_run = None, None
        
        to_copy = [anion_filename, cation_filename, topology_filename]
        #parse the topology file to find itp files to take
        self.itp_files = []
        with open(topology_filename) as f:
            dir = os.path.dirname(topology_filename)
            for line in f:
                m = re.match(r'#include "(\S+)\.itp"', line)
                if(m):
                    file = '{}.itp'.format(m.group(1))
                    if(os.path.exists(os.path.join(dir, file))):
                        self.itp_files.append(file)
                        to_copy.append(os.path.join(dir, file))

        
        self._setup_directory(*to_copy)
        if(initialize):
            self._pack(os.path.basename(anion_filename), anion_resname, os.path.basename(cation_filename), cation_resname, cation_number, BOX_DIM)
            self._solvate(BOX_DIM)

    def _putInDir(dirname):
        """
        A decorator that encloses a function inside a directory
        """
        def wrap(fxn):
            def mod_f(self, *args):
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                    #bring files
                for f in self.itp_files + [self.current_structure, self.current_top, self.current_run]:
                    if(f is not None):
                        shutil.copyfile(f, os.path.join(dirname, os.path.basename(f)))            
                #go there
                os.chdir(dirname)
                try:
                    fxn(self, *args)
                finally:
                    #make sure we leave
                    os.chdir(self.dir)
                    #bring back files
                    for f in self.itp_files + [self.current_structure, self.current_top, self.current_run]:
                        if(f is not None):
                            shutil.copyfile(os.path.join(dirname, f),f)
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
            f.write(textwrap.dedent('''
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
            
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulomb_type             = PME
            ; cut-off lengths
            rvdw                     = 1
            dispcorr                 = ener
            ewald_rtol               = 1e-5

            '''.format()))

        #prepare run
        self._exec_log('grompp', {'c':self.current_structure, 'f':input_file, 'p':self.current_top})
        #run
        emin_structure = 'emin_noion.gro'
        self.current_run = 'topol.tpr'
        self._exec_log(MDRUN, {'c':emin_structure, 's':self.current_run})
        self.current_structure = emin_structure
        
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
            
    def _pack(self, afile, ares, cfile, cres, cnumber, box):
        #update topology file
        with open(self.current_top, 'a') as f:
            f.write(textwrap.dedent('''
            [ system ]
            ;Name
            {name}

            [ molecules ]
            {ares:10s} {anumber:3d}
            {cres:10s} {cnumber:3d}
            '''.format(name=self.name, ares=ares, cres=cres, anumber=1, cnumber=cnumber)))
            
        #pack
        out_name = 'dry_mix.gro'
        args = {}
        args['p'] = self.current_top
        args['cp'] = afile
        args['ci'] = cfile
        args['o'] = out_name
        args['nmol'] = cnumber
        args['box'] = '{0} {0} {0}'.format(box)
        self._exec_log('genbox', args)
        self.current_structure = out_name

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
        args = {'cp':self.current_structure}
        out_name = 'mix.gro'
        args['o'] = out_name
        args['p'] = self.current_top
        args['cs'] = 'tip4p.gro'
        args['box'] = '{box} {box} {box}'.format(box=box)
        self._exec_log('genbox', args)
        self.current_structure = out_name
