import logging, os, shutil, datetime, subprocess, re, textwrap


ION_CONC=0#Molarity
BOX_DIM=3#nm
MDRUN_PROD='mpiexec.hydra -n 8 mdrun_mpi'
MDRUN='mpiexec.hydra mdrun_mpi'
SUM_HILLS='sum_hills'
EQUIL_TIME=2.5
PROD_TIME=40
DEBUG=True
HILL_HEIGHT=0.1
SIGMA=0.05

        
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
                 cation_atoms,
                 anion_atoms,
                 equil_time=EQUIL_TIME,
                 prod_time=PROD_TIME,
                 salt_conc=ION_CONC,
                 hill_height=HILL_HEIGHT,
                 sigma=SIGMA,
                 initialize=True):
        """
        Creates the directory, puts the cation and anion together and solvates
        """
        self.name = name        
        self.current_top = os.path.basename(topology_filename)        
        self.current_structure, self.current_run = None, None
        self.anion_resname = anion_resname
        self.cation_resname = cation_resname
        self.pressure = pressure
        self.equil_time = equil_time
        self.prod_time=prod_time
        self.ionc = salt_conc
        self.add_ions = False
        self.cation_number = cation_number
        if(salt_conc > 0 or cation_number != 1):
            self.add_ions = True
        self.cation_atoms = cation_atoms
        self.anion_atoms = anion_atoms
        self.hill_height = hill_height
        self.sigma = sigma
        self.plumed_files = ['COLVAR', 'HILLS', 'BIAS']
        
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
                for f in self.itp_files + self.plumed_files + [self.current_structure, self.current_top]:
                    if(f is not None and os.path.exists(f)):
                        shutil.copyfile(f, os.path.join(dirname, os.path.basename(f)))            
                #go there
                os.chdir(dirname)
                try:
                    fxn(self, *args)
                finally:
                    #make sure we leave
                    os.chdir(self.dir)
                    #bring back files
                    for f in self.itp_files + self.plumed_files + [self.current_structure, self.current_top]:
                        if(f is not None and os.path.exists(os.path.join(dirname, f))):
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
            integrator               = steep
            emtol                    = 0.01
            nsteps                   = 10000; max steps
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 5
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME
            rvdw            = 1
            rlist           = 1
            rcoulomb        = 1
            optimize-fft    = yes
            fourierspacing  = 0.12
            pme-order       = 4
            ewald-rtol      = 1e-5
            dispcorr                 = ener
            ;constraints              = h-bonds

            '''.format()))

        #prepare run
        self._exec_log('grompp', {'c':self.current_structure, 'f':input_file, 'p':self.current_top})
        #run
        emin_structure = 'emin_noion.gro' if self.add_ions else 'emin.gro'
        self.current_run = 'topol.tpr'
        if(DEBUG):
            shutil.copyfile(self.current_structure, emin_structure)
        else:
            self._exec_log(MDRUN, {'c':emin_structure, 's':self.current_run})
        self.current_structure = emin_structure

        #now add ions if needed
        if(self.add_ions):
            emin_structure = 'emin_noion.gro'
            arg_dic={'s':self.current_run, 
                     'o':emin_structure, 
                     'neutral':'', 
                     'conc':self.ionc}


            #get charge
            toadd = self.cation_number - 1
            if(toadd > 0):
                arg_dic['nn'] = toadd
            
            #get the group for SOL. Requires a read and kill to find group ordering
            match = self._read_and_kill(string='genion', arg_dic=arg_dic,
                                        read_fxn=lambda x: re.match('\s*Group\s*(\d+)\s*\(\s*SOL.*', x))
            sol_group = match.group(1)
            self._exec_log('genion', arg_dic, sol_group)
            self.current_structure = emin_structure



    @_putInDir('equil')
    def equilibrate(self):
        """
        Equilibrate the system in NPT
        """
        #build input file
        input_file = 'equil.mdp'
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = sd
            nsteps                   = {time}
            dt                       = 0.002
            
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 1
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            

            
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME
            rvdw            = 1
            rlist           = 1
            rcoulomb        = 1
            optimize-fft    = yes
            fourierspacing  = 0.12
            pme-order       = 4
            ewald-rtol      = 1e-5
            dispcorr                 = ener

            ;Temperature
            tcoupl                   = v-rescale
            tau_t                    = 2
            ref_t                    = 300
            tc-grps                  = System
            
            ;Pressure            
            pcoupl                   = berendsen
            tau_p                    = 2
            compressibility          = 4.5e-5
            ref_p                    = {pressure}

            ;Constraints
            constraints              = h-bonds

            '''.format(pressure=self.pressure * 1.01325, 
                       time=self.equil_time * 10**6 / 2.)))

        #prepare run
        self._exec_log('grompp', {'c':self.current_structure, 
                                  'f':input_file, 
                                  'p':self.current_top})
        #run
        equil_structure = 'equil.gro'
        self.current_run = 'topol.tpr'
        if(DEBUG):
            shutil.copyfile(self.current_structure, equil_structure)
        else:
            self._exec_log(MDRUN_PROD, {'c':equil_structure, 's':self.current_run})
        self.current_structure = equil_structure


    @_putInDir('prod')
    def production(self):
        """
        Production run in NVT with metadynamics
        """

        #First we need to convert the CV atom indices to match our system.
        #read through gro file and connect residue names with what we packed the system with
        def convert_indices(resname, indices):
            result = []
            with open(self.current_structure) as f:
                atom_index = 0            
                residue_index = 0
                for line in f:                
                    regex = '\s*(\d*){}\s*\w*\s*(\d*).*'.format(resname)
                    m = re.match(regex, line)
                    if(m):
                        if(m.group(1) != residue_index):
                            residue_index = m.group(1)
                            atom_index = 1
                        else:
                            atom_index += 1
                    if(atom_index in indices):
                        result.append(int(m.group(2)))
            return result
        cation_indices = convert_indices(self.cation_resname, self.cation_atoms)
        anion_indices = convert_indices(self.anion_resname, self.anion_atoms)
                            
        #build input file
        input_file = 'prod.mdp'
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = sd
            nsteps                   = {time}
            dt                       = 0.002
            
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 1000
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            

            
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME
            rvdw            = 1
            rlist           = 1
            rcoulomb        = 1
            optimize-fft    = yes
            fourierspacing  = 0.12
            pme-order       = 4
            ewald-rtol      = 1e-5
            dispcorr                 = ener

            ;Temperature
            tcoupl                   = v-rescale
            tau_t                    = 2
            ref_t                    = 300
            tc-grps                  = System
            
            ;Constraints
            constraints              = h-bonds

            '''.format(pressure=self.pressure * 1.01325, 
                       time=self.prod_time * 10**6 / 2.)))

        #plumed file
        plumed_file = 'plumed.dat'
        with open(plumed_file, 'w') as f:
            f.write(textwrap.dedent('''
            HILLS HEIGHT {hill_height} W_STRIDE 1000
            WELLTEMPERED SIMTEMP 300 BIASFACTOR 5
            PRINT W_STRIDE 500
            DISTANCE LIST <anion> <cation> SIGMA {sigma}
            GRID CV 1 MIN 0 MAX 3.5 NBIN 300 
            WRITE_GRID FILENAME BIAS W_STRIDE 10000
            anion->
            {anion}
            anion<-
            cation->
            {cation}
            cation<-
            ENDMETA
            '''.format(hill_height=self.hill_height, 
                       sigma=self.sigma, 
                       cation=reduce(lambda x,y:'{} {}'.format(x,y),
                                     cation_indices),
                       anion=reduce(lambda x,y:'{} {}'.format(x,y),
                                    anion_indices))))
                    
            #prepare run
            self._exec_log('grompp',
                       {'c':self.current_structure, 
                        'f':input_file, 
                        'p':self.current_top})
        #run
        prod_structure = 'prod.gro'
        self.current_run = 'topol.tpr'
        if(DEBUG):
            shutil.copyfile(self.current_structure, prod_structure)
        else:
            self._exec_log(MDRUN_PROD, 
                           {'c':prod_structure, 
                            's':self.current_run,
                            'plumed':plumed_file})
        self.current_structure = prod_structure

    @_putInDir('analysis')
    def analyze(self):
        #get PMF
        self._exec_log(SUM_HULLS,{
            'file': 'HILLS',
            'ndw':1,
            'ncv': 1,
            'kt': 2.494,
            'out':'pmf.dat'})

        #plot that sucker
        import matplotlib as plt
        import numpy as np
        plt.figure()
        data = np.genfromtxt('pmf.dat')
        plt.plot(data[:,1], data[:,2])
        plt.xlabel('r [nm]')
        plt.ylabel(r'$\Delta$ A [kJ/mol]')
        plt.savefig('pmf.png')
        
    def _setup_directory(self, *to_copy):
        #build directory, start log
        if not os.path.exists(self.name):
            os.mkdir(self.name)
        self.dir = os.path.abspath(self.name)
        for f in to_copy:
            shutil.copyfile(f, os.path.join(self.dir, os.path.basename(f)))
        os.chdir(self.dir)
        #remove other logging handlers
        log = logging.getLogger()
        for hdlr in log.handlers:
            log.removeHandler(hdlr)        
        fileh = logging.FileHandler('simulation.log', 'a')
        fileh.setLevel(logging.DEBUG)
        log.addHandler(fileh)
        logging.info('Beginning simulation {}'.format(datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')))
            
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

    def _read_and_kill(self, string, read_fxn, arg_dic=None):
        """
        Execute the given string, reads its input and kills it
        
        Arguments:
        string --- the command to execute
        arg_doc --- additional arguments to add to the command string
        read_fxn -- a function which is called on each line of the output. It should return None when
                        it needs to process additional lines. The first none None value is returned
                        and the process is terminated.
        """
        if(arg_dic is not None):
            for k,v in arg_dic.iteritems():
                string += ' -{} {}'.format(k,v)
        try:
            
            from threading import Thread
            from Queue import Queue, Empty
            #wait for process to stop. hack
            def enqueue_output(out, queue):
                for line in iter(out.readline, b''):
                    queue.put(line)
                out.close()

            process = subprocess.Popen(string, shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)

            q = Queue()
            #SO FUCKING ANNOYING, gromacs uses stderr for communicating. So dumb
            t = Thread(target=enqueue_output, args=(process.stderr, q))
            t.daemon=True
            t.start()

            logging.debug('Reading from cmd {}'.format(string))
            while process.poll() is None:
                try: 
                    line = q.get(timeout=1)
                    logging.debug('line: {}'.format(line))
                    v = read_fxn(line)
                    if(v is not None):
                        process.terminate()
                        return v
                except Empty: continue


        except OSError as e:
            print >>sys.stderr, 'Execution of {} failed:'.format(string), e

        logging.error('Command {} failed'.format(string)) 
        raise RuntimeError()
            


    def _exec_log(self, string, arg_dic=None, input=None):
        """
        Smartly executes the given string and logs the results
        """
        if(arg_dic is not None):
            for k,v in arg_dic.iteritems():
                string += ' -{} {}'.format(k,v)
        try:
            process = subprocess.Popen(string, shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       stdin=subprocess.PIPE)
                           
            if(input is None):
                out, err = process.communicate()
            else:
                out, err = process.communicate(input)
            retcode = process.returncode
            logging.debug(out)
            if(len(err) > 0):
                logging.warning(err)
                logging.debug(out)
            if retcode < 0:
                logging.error('{} failed with {}. stderr follow'.format(string, retcode))
                logging.error(err)
                raise OSError('{} failed with {}. stderr follow'.format(string, retcode))
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
