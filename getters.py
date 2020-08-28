from pytriqs.archive import *
from pytriqs.gf import *
import numpy

from amoeba import amoeba

#------------------ basic -------------------------------------------#

def get_dos(eps,t,L, wmax = 1.1, nw = 2201):
    ks = numpy.linspace(0,2.*numpy.pi,L,endpoint=False)
    epsk = eps+2*t*(numpy.cos(ks)[None,:]+numpy.cos(ks)[:,None])
    ws = numpy.linspace(-wmax,wmax,nw,endpoint=True)
    dw = ws[1]-ws[0]
    ws-=dw/2.
    dos = 0*the_ws
    for e in epsk.nditer():
        dos[int(math.floor((e - ws[0])/dw))]+=1
    ws+=dw/2.
    return ws, dos/(dw*L**2)

def get_gauge_field_dos(eps,t,n,L, wmax = 1.1, nw = 2201, optimize_kys=True):
    #implements Eq.60 by constructing temporary H_k matrix to save memory.
    # USE THIS for large L.
    # check q=L case: does the spectrum depend on k_y? If not, optimize the code in this respect!
    gcd = computeGCD(n, L)
    p = n/gcd
    q = L/gcd
    print "L:",L,"n:",n,"gcd:",gcd,"p:",p,"q:",q
    ks = numpy.linspace(0,2.*numpy.pi,L,endpoint=False)
    kxs = ks[:gcd]
    kys = ks if ((q!=L)or(not optimize_kys)) else numpy.array([0])
    
    ws = numpy.linspace(-wmax,wmax,nw,endpoint=True)
    dw = ws[1]-ws[0]
    ws-dw/2. 

    dos = numpy.zeros((nw),dtype=numpy.int_)
    
    mat1 = numpy.roll(numpy.eye(q), -p, axis=1)
    mat2 = numpy.roll(numpy.eye(q), p, axis=1)
    for kxi,kx in enumerate(kxs):
        kxls = kx+numpy.arange(q)*(2.*numpy.pi/q)
        Hdiag = numpy.diag(eps+2.0*t*numpy.cos(kxls))
        for kyi,ky in enumerate(kys):       
            Hllp = Hdiag+t*(numpy.exp(1j*ky)*mat1+numpy.exp(-1j*ky)*mat2)    
            Es = numpy.linalg.eigvalsh(Hllp)
            for E in Es:
                dos[int(math.floor((E - ws[0])/dw))]+=1

    return ws+dw/2.,dos/(dw*(L*len(kys)))

#----------------------------------------------------------------------------#

#------------------ various Dyson -------------------------------------------#

def Hilbert(G_iw, ws, dos, mu, Sigma_iw, blocks=['up','dn']):
    dw = ws[1]-ws[0]
    for b in blocks:  
        iws = [iw.value for iw in g.mesh]
        G_iw[b].data[:,0,0] = dw*numpy.sum(
          dos[:,None]\
          / ( iws[None,:] + mu[b] - ws[:,None] - Sigma_iw[b].data[None,:,0,0] ),
          axis=(0)
        )           

def orbital_space_dyson_get_G(G_ij_iw, G0_ij_iw, Sigma_ij_iw):
    G_ij_iw << inverse(inverse(G0_ij_iw)-Sigma_ij_iw)

def orbital_space_dyson_get_Sigma(Sigma_ij_iw, G0_ij_iw, G_ij_iw):
    Sigma_ij_iw << inverse(G0_ij_iw)-inverse(G_ij_iw)

def orbital_space_dyson_get_G0(G0_ij_iw, G_ij_iw, Sigma_ij_iw, mu_shift = 0 ):
    G0_ij_iw << inverse(mu_shift + inverse(G_ij_iw)+Sigma_ij_iw)

#---------------------mu search--------------------------------------------------------#

def search_for_mu(get_mu, set_mu, get_n, n, ph_symmetry, accepted_mu_range=[-20.0,20.0]):  
  print "getters: search_for_mu:  n: ",n,", ph_symmetry",ph_symmetry, "accepted mu_range: ",accepted_mu_range

  if (n is None) or ((n==0.5) and ph_symmetry):
    print "no mu search to be performed! it is your duty to set the chemical potential to U/2. mu =",get_mu()
    print 'n on the lattice : ', get_n()
  else:
    def func(var, data=None):
      mu = var[0]        
      set_mu(mu)
      actual_n = get_n()
      val = 1.0-abs(n - actual_n)  
      print "amoeba func call: mu: %.2f desired n: %.2f actual n: %.2f val = "%(mu,n,actual_n),val
      if val != val: return -1e+6
      else: return val

    print "about to do mu search:"

    guesses = [get_mu(), 0.0, -0.1, -0.3, -0.4, -0.5, -0.7, 0.3, 0.5, 0.7]
    found = False  
    for l in range(len(guesses)):
      varbest, funcvalue, iterations = amoeba(
        var=[guesses[l]],
        scale=[0.01],
        func=func, 
        data = None,
        itmax=30,
        ftolerance=1e-2,
        xtolerance=1e-2,
        known_max = 1.0,
        known_max_accr = 5e-5
      )
      if ( varbest[0]>accepted_mu_range[0] 
           and 
           varbest[0]<accepted_mu_range[1]
         ) 
         and 
         ( abs(funcvalue-1.0)<1e-2 ): #change the bounds for large doping
        found = True 
        func(varbest)
        break 
      if l+1 == len(guesses):
        print "mu search FAILED: doing a scan..."

        mu_grid = numpy.linspace(accepted_mu_range[0],accepted_mu_range[1],50)
        func_values = [func(var=[mu]) for mu in mu_grid]
        print "func_values: "
        for i in range(len(mu_grid)):
          print "mu: ",mu_grid[i], " 1-abs(n-n): ", func_values[i]
        mui_max = numpy.argmax(func_values)
        print "using mu: ", mu_grid[mui_max]
        set_mu(mu_grid[mui_max])  
        get_n()
           
    if found:
      print "guesses tried: ", l  
      print "mu best: ", varbest
      print "1-abs(diff n - data.n): ", funcvalue
      print "iterations used: ", iterations

