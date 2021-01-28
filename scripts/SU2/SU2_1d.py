# general modules -----------------------------
import sys,os
import numpy as np
import math
import argparse
import csv

# plottig modules and setting --------------------
import matplotlib.pyplot as plt
from matplotlib import rc
rc('axes', linewidth=3)
rc('text', usetex=True)
rc('font', size=15, family='serif', weight='bold')
rc('lines',linewidth=3,color='blue')
rc('xtick.major', pad=10, size=5, width=3)
rc('ytick.major', pad=10, size=5, width=3)

# user basis modules -----------------------------
from quspin.operators import hamiltonian, commutator                                # operators to implement
from quspin.basis.user import user_basis     # Hilbert space user basis
from quspin.basis import tensor_basis,spinful_fermion_basis_1d,spin_basis_1d
from quspin.basis.user import pre_check_state_sig_32,op_sig_32,map_sig_32           # user basis data types
from numba import carray,cfunc # numba helper functions
from numba import uint32,int32 # numba data types
from time import time # timing package

# read in parameters
parser=argparse.ArgumentParser()
parser.add_argument('--N', help='number of sites')
parser.add_argument('--t', help='fermion hopping')
parser.add_argument('--m', help='fermion mass')
parser.add_argument('--g', help='coupling')
parser.add_argument('--OMP', help='(optional) no. openmp threads (default 1)',default=1)
parser.add_argument('--MKL', help='(optional) no. of MKL threads (default 1)',default=1)
parser.add_argument('--BC', help='boundary conditions OBC (O) or PBC (P), default open',default='O')

args=parser.parse_args()

# OMP and MKL parameters
os.environ['KMP_DUPLICATE_LIB_OK']='True' # uncomment this if omp error occurs on OSX for python 3
os.environ['OMP_NUM_THREADS']=str(int(args.OMP)) # number of OpenMP threads
os.environ['MKL_NUM_THREADS']=str(int(args.MKL)) # number of MKL threads

N=int(args.N)                     # number of sites
tf=float(args.t)                   # hopping parameter
mf=float(args.m)                   # mass
g=float(args.g)                   # coupling
#site = int(args.site)           # DEBUG only
#site1 = int(args.site1)           # DEBUG only

if args.BC=='O':
    print("# boundary conditions: OBC")
elif args.BC=='P':
    sys.exit("# ERROR: PBC not yet supported")
else:
    sys.exit("# ERROR boundary conditions are not specified")


# output parameters
folder = 'data_SU2/'
folder_save = folder + 'eigenvalues_N'+str(N)+'_g'+str(g)+'_tf'+str(tf)+'_mf'+str(mf)

if not os.path.exists(folder_save):
    os.makedirs(folder_save)

# info file
infoF = open(folder_save+'/info.txt', 'w')
for arg in vars(args):
    infoF.write(arg)
    infoF.write("\t")
    infoF.write(str(getattr(args, arg)))
    infoF.write("\n")
infoF.close()


# --------------------------------------------------------------------------- #
#   BUILDING BASIS - Fermion basis with two components
# ----------------------------------------------------------------------------#
ti = time()

# ------- FERMION BASIS ------------------------------------------------------#
# Basis: use internal fermion basis of quspin
fermion_basis = spinful_fermion_basis_1d(L=N,Nf=(N//2,N//2))
#fermion_basis = spinful_fermion_basis_1d(L=N,Nf=(1,0))

print("# Number of states (half filling), OBC ", fermion_basis.Ns)
print("# TIME: building Fermion basis took {0:.4f} sec".format(time()-ti))

# -------------------------------------------------------------#
#   BUILD HAMILTONIAN
# -------------------------------------------------------------#
ti = time()

# FERMION BASIS ---------------------------------------------- #
static = []
dynamic = []

# hopping term
if np.abs(tf)>0.0:
    ht = []
    ht_conj = []
    for n in range(N-1):
        ht.append([-tf,n,n+1])
        ht_conj.append([tf,n,n+1])
    
    static.append(["+-|",ht])
    static.append(["-+|",ht_conj])
    static.append(["|+-",ht])
    static.append(["|-+",ht_conj])
    
    print(static)
  
# mass term
if np.abs(mf)>0:
    hm = []
    for n in range(N):
        hm.append([mf*(-1.0)**n,n,n])
        
    static.append(["+-|",hm])
    static.append(["|+-",hm])

H =  hamiltonian(static, dynamic, basis=fermion_basis, dtype=np.float64)

print("# TIME: building fermion Hamiltonian took {0:.4f} sec".format(time()-ti))

# -------------------------------------------------------------#
#   Diagonalization
# -------------------------------------------------------------#
ti = time()

E,V = H.eigh()

print("# TIME: Diagonalization in {0:.4f} sec".format(time()-ti))

with open(folder_save+'/eigenvalues.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(np.arange(0,len(E)),E))

