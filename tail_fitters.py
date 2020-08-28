import numpy
from pytriqs.gf import *

def fit_fermionic_sigma_tail(
    Q, 
    starting_iw=14.0, 
    no_loc=False, 
    hartree=None, 
    weight=None,
    max_order=5, 
    overwrite_tail=True
):
    Nc = len(Q.data[0,:,0])    
    if no_loc:
        known_moments = numpy.zeros((2,Nc,Nc),dtype=numpy.complex)    
    elif (weight is None) and (hartree is None):
        known_moments = numpy.zeros((0,Nc,Nc),dtype=numpy.complex) 
    elif weight is None:
        known_moments = numpy.zeros((1,Nc,Nc),dtype=numpy.complex) 
        known_moments[0,:,:] = hartree*numpy.eye(Nc)
    else:
        known_moments = numpy.zeros((2,Nc,Nc),dtype=numpy.complex) 
        known_moments[0,:,:] = hartree*numpy.eye(Nc)
        known_moments[1,:,:] = weight*numpy.eye(Nc)
      
    nmax = Q.mesh.last_index()
    nmin = int(((starting_iw*Q.beta)/numpy.pi-1.0)/2.0) 
    before = Q.data[-1,0,0]
    tail = Q.fit_hermitian_tail_on_window(nmin, nmax, known_moments, n_tail_max=nmax-nmin, expansion_order=max_order)
    #print "tail:",tail
    if overwrite_tail: Q.replace_by_tail_in_fit_window(tail[0])
    if (Q.data[-1,0,0]==before) and overwrite_tail:
        print "[ Node", mpi.rank, "]: fit_fermionic_sigma_tail: WARNING: negative part of the tail may not have been overwritten!!!"
    return tail

def fit_fermionic_g_tail(
    Q, 
    starting_iw=14.0, 
    max_order=5, 
    overwrite_tail=True
):
    fit_fermionic_sigma_tail(
        Q, 
        starting_iw=starting_iw, 
        no_loc=False, 
        hartree=0., 
        weight=1.,
        max_order=max_order, 
        overwrite_tail=overwrite_tail
    )
