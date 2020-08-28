import pytriqs.utility.mpi as mpi

import data
from data import *
from data import data
from tail_fitters import *
from generic_loop import generic_action, generic_loop, mixer, converger, monitor
from getters import *
from impurity_solvers import solvers
from cautionaries import *

def magnetic_dmft_data( 
  beta,
  niw,
  ntau,
  blocks = ['up','dn']
):
  dt = data() 
  dt.beta = beta
  dt.T = 1./beta
  assert niw % 2 == 0, "must be even number"
  dt.niw = niw
  dt.n_iw = niw/2
  dt.ntau = ntau
  dt.blocks = blocks 
  dt.ns = {b: 0 for b in blocks}
  dt.mus = {b: 0 for b in blocks}  
  AddGfData(dt, ['G_imp_iw', 'G_loc_iw', 'Sigma_imp_iw', 'Gweiss_iw'], blocks, 1, dt.n_iw, beta, domain = 'iw', suffix='', statistic='Fermion')
  dt.iws = numpy.array([iw.value for iw in dt.G_imp_iw.mesh])   
  #AddGfData(dt, ['Sigma_imp_tau',' Gweiss_tau'], blocks, 1, ntau, beta, domain = 'tau', suffix='', statistic='Fermion')
  return dt
  
def get_Sigma(
  solver,
  Sigma_imp_iw,
  Gweiss_iw, 
  U,
  max_time=5*60,
  delta=0.1,
  solver_data_package=None
):
  solver.G0_iw << Gweiss_iw 
  solvers.ctint.run(
    solver=solver, 
    U=U,      
    alpha=0.5,
    delta=delta,
    n_cycles=100000000,
    max_time = max_time,
    solver_data_package = solver_data_package,
    only_sign = False
  )
  Sigma_imp_iw << solver.Sigma_iw

def magnetic_dmft_set_up_calc( dt, max_time=5*60, delta=0.1, solver_data_package=None ):
  dt.get_Gweiss = lambda: orbital_space_dyson_get_G0(dt.Gweiss_iw, dt.G_loc_iw, dt.Sigma_imp_iw)

  dt.get_Sigma = lambda: get_Sigma(dt.solver, dt.Sigma_imp_iw, dt.Gweiss_iw, dt.U, max_time, delta, solver_data_package)

  dt.get_G_loc = lambda: Hilbert(dt.G_loc_iw, dt.ws, dt.dos, dt.mus, dt.Sigma_imp_iw, blocks=dt.blocks)

  dt.get_mu = lambda: dt.mu

  def set_mu(mu):
    dt.mu = mu
    for b in dt.blocks:
      dt.mus[b] = dt.hs[b]+mu

  dt.set_mu = set_mu

  def get_n():  
    dt.get_G_loc()
    for b in dt.blocks:
      fit_fermionic_g_tail(dt.G_loc_iw[b], starting_iw=20.0,  max_order=5, overwrite_tail=True)
      dt.ns[b] = dt.G_loc_iw[b].density()[0,0].real
      print "inside get_n: ns:", dt.ns
    return sum(dt.ns.values()) 

  dt.get_n = get_n

def magnetic_dmft_set_params_and_initialize(
  dt,   
  U,
  mu,
  n,
  fixed_n, 
  h, #Zeeman splitting
  p, q, ncells=1, # gauge field settings: B_z ~ p/q, q size of the unit cell, L size of the lattice = q*ncells
  ph_symmetry=False,
  t=-0.25,
  initial_guess='metal',
  filename=None,
  solver_data_package=None
):
  dt.U = U
  dt.t = t
  dt.h = h
  dt.p = p
  dt.q = q
  dt.L = q*ncells
  dt.hs = {'up': h, 'dn': -h}
  dt.set_mu(mu)
  dt.fixed_n = fixed_n
  dt.n=n 
  dt.ph_symmetry = ph_symmetry  
  dt.iteration = 0

  if filename is None: 
    filename = "dmft.U%g.T%g.h%g"\
                %(U,dt.T,h)
    if p>0: filename += ".p%d.q%d.L%d"%(p,q,L)
    if fixed_n: filename += ".n%g"%n
    else: filename += ".mu%g"%mu
    filename += "_from_%s"%initial_guess

  dt.archive_name = filename
  dt.dump = lambda dct: DumpData(dt, filename, Qs=[], exceptions=['solver'], dictionary=dct)
  dt.dump_final = dt.dump

  dt.solver = solvers.ctint.initialize_solver(
    solver_data_package = solver_data_package,  
    beta = dt.beta,
    nsites = 1,
    n_iw = dt.n_iw,
    ntau = max(dt.ntau, 100001)
  )
  if p==0:
    dt.ws, dt.dos = get_dos(0, t, 2000, wmax = 4.1*t, nw = 2001)
  else:
    dt.ws, dt.dos = get_gauge_field_dos(0, t, p*L/q, L, wmax = 4.1*t, nw = 2001, optimize_kys=True)

  print "Filling Sigma_imp_iw.."
  dt.Sigma_imp_iw << U*0.5-int(initial_guess=='atomic')*inverse(iOmega_n)
  
  print "Getting Gloc and n on the lattice: ",dt.get_n()
  print "Getting Gweiss.."
  dt.get_Gweiss() 
    
  print "Done initializing, about to dump..."
  dt.dump('initial')

def magnetic_dmft_actions(dt, accr):
  def impurity(dt):
    dt.get_Sigma()

  def lattice(dt):
    dt.iteration += 1

    if dt.fixed_n:
      search_for_mu( dt.get_mu, dt.set_mu, dt.get_n, dt.n, dt.ph_symmetry ) 
    else:     
      print "fixed mu calculation, doing G"
      dt.n = dt.get_n() 
      print "n(G_loc) =",dt.n

  def pre_impurity(dt):
    dt.get_Gweiss() 

  actions = [
    generic_action( 
      name = "pre_impurity",
      main = pre_impurity,
      mixers = [],#[lambda data, it: 0],
      cautionaries = [ ],#[lambda data, it: 0], allowed_errors = [],               
      printout = lambda data, it: 0,
      short_timings = True
    ),
    generic_action( 
      name = "impurity",
      main = impurity,
      mixers = [], #[lambda data, it: 0],
      cautionaries = [ lambda data, it: impose_real_valued_in_imtime_blockGf(data.Sigma_imp_iw) ],#[lambda data, it: 0], allowed_errors = [],            
      printout = lambda data, it: ( 
        solvers.ctint.dump(data.solver, dt.archive_name, suffix=it)
      ),
      short_timings = True
    ),
    generic_action( 
      name = "lattice",
      main = lattice,
      mixers = [],#[lambda data, it: 0],
      cautionaries = [],#[lambda data, it: 0], allowed_errors = [],               
      printout = lambda data, it: data.dump(it),# if (int(it[-3:])%5==0) else 0),
      short_timings = True
    )
  ]

  monitors = [
    monitor(
      monitored_quantity = lambda: dt.solver.average_sign, 
      h5key = "average_sign_vs_it", 
      archive_name = dt.archive_name
    ),
    monitor(
      monitored_quantity = lambda: dt.mu, 
      h5key = "mu_vs_it", 
      archive_name = dt.archive_name
    ),
    monitor(
      monitored_quantity = lambda: sum(dt.ns.values()), 
      h5key = "n_vs_it", 
      archive_name = dt.archive_name
    )
  ]\
  +[
    monitor(
      monitored_quantity = lambda q=q, b=b: vars(dt)[q][b], 
      h5key = q[:1]+b+"_vs_it", 
      archive_name = dt.archive_name
    )
   for q in ['ns','mus'] for b in dt.blocks
  ]

  convergers = [
    converger( monitored_quantity = lambda: dt.Sigma_imp_iw, accuracy=accr, func=None, archive_name=dt.archive_name, h5key='diffs_Sigma_imp_iw'),
    converger( monitored_quantity = lambda: dt.G_loc_iw, accuracy=accr, func=None, archive_name=dt.archive_name, h5key='diffs_G_loc_iw'),
  ]

  return actions, monitors, convergers

def magnetic_dmft_launcher(
  U,
  T,
  mu, #used as initial mu for mu search if fixed_n. otherwise, that's the mu.
  n, #desired n if fixed_n. otherwise whatever
  fixed_n, 
  h, #Zeeman splitting
  p, q, ncells=1, # gauge field settings: B_z ~ p/q, q size of the unit cell, L size of the lattice = q*ncells
  ph_symmetry=False,
  t=-0.25,

  initial_guess='metal',

  n_cycles=10000000,      
  delta = 0.5,        
  max_time_rules= [ [1, 5*60], [2, 20*60], [4, 80*60], [8, 200*60], [16,400*60] ], 
  time_rules_automatic=False, 
  exponent = 0.7, 
  overall_prefactor=1.0,
  no_timing = False,

  max_its = 20,
  min_its = 5, 
  accr = 1e-3,  

  iw_cutoff=40.0,

  filename=None
):
  solver_data_package={}

  if mpi.is_master_node():
    print "------------------- Welcome to MAGNETIC DMFT! -------------------------"
    beta = 1.0/T
    niw = 2*int(((iw_cutoff*beta)/numpy.pi-1.0)/2.0)
    ntau = 3*niw
    print "Automatic niw:",niw   

    if no_timing:
      max_time=-1
      print "no timing!!!"
    else:
      if time_rules_automatic:
        pref = ((beta/8.0)*U)**exponent #**1.2
        print "pref: ",pref 
        max_time = int(overall_prefactor*pref*5*60)
        print "max times automatic: ",max_time
      else:
        max_time = max_time_rules[0][1]
        print "max_time from rules: ",max_time

    dt = magnetic_dmft_data( 
      beta, niw, ntau
    )

    magnetic_dmft_set_up_calc( 
      dt, max_time,delta, solver_data_package 
    ) 

    magnetic_dmft_set_params_and_initialize(
      dt, U, mu, n, fixed_n, 
      h, #Zeeman splitting
      p, q, ncells, # gauge field settings: B_z ~ p/q, q size of the unit cell, L size of the lattice = q*ncells
      ph_symmetry, t, initial_guess, filename,
      solver_data_package
    )

    actions, monitors, convergers = magnetic_dmft_actions(dt, accr)

    magnetic_dmft = generic_loop(
      name = "Magnetic DMFT loop", 
      actions = actions,
      convergers = convergers,  
      monitors = monitors
    )

    magnetic_dmft.run(
      dt, 
      max_its = max_its,
      min_its = min_its, 
      max_it_err_is_allowed = 7,
      print_final = True,
      print_current = 1,
      start_from_action_index = 1
    )
    if mpi.size>1:
      solver_data_package['tag'] = 'exit'
      solver_data_package = mpi.bcast(solver_data_package)
    return dt
  else:
    solvers.ctint.slave_run(solver_data_package=solver_data_package, printout=False)




