from pytriqs.operators import *
from pytriqs.archive import *
#from pytriqs.gf.local import *
#from pytriqs.arrays import BlockMatrix, BlockMatrixComplex
import pytriqs.utility.mpi as mpi

from copy import deepcopy

try:
  #from triqs_ctint import SolverCore as Solver
  from triqs_ctint import Solver
except:
  if mpi.is_master_node():
    print "CTINT not installed"

from slave_run import slave_run

import copy

################################ IMPURITY #########################################

class solvers:
  class ctint:
    @staticmethod
    def initialize_solver(
      solver_data_package = None,  
      beta = None,
      nsites = None,
      n_iw = None,
      ntau = 100000, 
    ):
      if solver_data_package is None: solver_data_package = {}    

      gf_struct = [ (b, range(nsites)) for b in ['up','dn'] ]

      assert ntau>2*n_iw, "solvers.ctint.initialize_solvers: ERROR! ntau too small!!" 

      solver_data_package['constructor_parameters']={}
      solver_data_package['constructor_parameters']['beta'] = beta
      solver_data_package['constructor_parameters']['n_iw'] = n_iw
      solver_data_package['constructor_parameters']['n_tau'] = ntau
      solver_data_package['constructor_parameters']['gf_struct'] = gf_struct
      solver_data_package['tag'] = 'construct'

      if mpi.is_master_node(): print "solver_data_package:", solver_data_package  

      if mpi.size>1: solver_data_package = mpi.bcast(solver_data_package)

      return Solver( **solver_data_package['constructor_parameters'] )

    @staticmethod
    def run(
      solver, 
      U,      
      alpha=0.5,
      delta=0.1,
      n_cycles=20000,
      max_time = 5*60,
      solver_data_package = None,
      only_sign = False
    ):

      block_names = [name for name,g in solver.G0_iw]
      assert len(block_names)==2, "we need two blocks!!"
      N_states = len(solver.G0_iw[block_names[0]].data[0,0,:])      
      
      gf_struct = {}
      for bn in block_names:
        gf_struct[bn] = range(N_states)

      h_int = U * n(block_names[0],0)*n(block_names[1],0)
      for i in range(1,N_states):
        h_int += U * n(block_names[0],i)*n(block_names[1],i)

      N_s = 1      
      ALPHA = [ [ [ alpha + delta*(-1)**(s+sig) for s in range(N_s)] for i in range(N_states)] for sig in range(2) ]

      if solver_data_package is None:  solver_data_package = {}    

      solver_data_package['solve_parameters'] = {}
      solver_data_package['solve_parameters']['U'] = U
      #solver_data_package['solve_parameters']['alpha'] = ALPHA
      solver_data_package['solve_parameters']['n_cycles'] = n_cycles
      solver_data_package['solve_parameters']['max_time'] = max_time
      solver_data_package['solve_parameters']['length_cycle'] = 50
      solver_data_package['solve_parameters']['n_warmup_cycles'] = 2000
      solver_data_package['solve_parameters']['measure_M_tau'] = True
      solver_data_package['solve_parameters']['post_process'] = True
      solver_data_package['solve_parameters']['measure_histogram'] = True

      print solver_data_package['solve_parameters']
       
      solver_data_package['G0_iw'] = solver.G0_iw

      solver_data_package['tag'] = 'run'

      if mpi.size>1: 
         if mpi.is_master_node(): print "broadcasting solver_data_package!!"
         solver_data_package = mpi.bcast(solver_data_package)

      if mpi.is_master_node(): print "about to run "
      dct = deepcopy(solver_data_package['solve_parameters'])
      del dct['U']
      try:
        solver.solve(h_int = h_int, **dct )
      except Exception as e:
        A = HDFArchive('black_box','w')
        A['solver']=solver
        del A
        raise e
      if mpi.is_master_node(): print "average sign: ",solver.average_sign


    @staticmethod  
    def slave_run(solver_data_package, printout=True, additional_tasks = {}):
      internal_data = {}
      def construct(solver_data_package):
        if printout: print "[Node ",mpi.rank,"] constructing solvers!!!"
        internal_data['solver'] = Solver( **(solver_data_package['constructor_parameters']) )
        internal_data['gf_struct'] = solver_data_package['constructor_parameters']['gf_struct']

      def run(solver_data_package):
        solver = internal_data['solver']
        gf_struct = internal_data['gf_struct']

        solver.G0_iw << solver_data_package['G0_iw']
        U = solver_data_package['solve_parameters']['U']
        block_names = [name for name,g in solver.G0_iw]
        N_states = len(solver.G0_iw[block_names[0]].data[0,0,:])
        h_int = U * n(block_names[0],0)*n(block_names[1],0)
        for i in range(1,N_states):
          h_int += U * n(block_names[0],i)*n(block_names[1],i)   
        try:
          dct = deepcopy(solver_data_package['solve_parameters'])
          del dct['U']
          if printout: print "[Node ",mpi.rank,"] about to run..."
          solver.solve(h_int = h_int, **dct )

          if printout: print "[Node ",mpi.rank,"] finished running successfully!"

        except Exception as e:
          print "[Node ",mpi.rank,"] ERROR: crash during running solver" 

      tasks = {
        'construct': construct,
        'run': run
      }
      tasks.update(additional_tasks)
      slave_run(solver_data_package, printout=False, tasks = tasks)
         

    @staticmethod
    def dump(solver, archive_name, suffix=''):    
      dct = {
        'average_sign': solver.average_sign,
        'histogram': solver.histogram,
        'density': solver.density,
        'G_iw': solver.G_iw,
        'Sigma_iw': solver.Sigma_iw,
        'G0_iw': solver.G0_iw,
        'G0_shift_iw': solver.G0_shift_iw,
        'M_tau': solver.M_tau,
        'M_iw': solver.M_iw,

      }     
      A = HDFArchive(archive_name)
      A['solver%s'%suffix] = dct
      del A
