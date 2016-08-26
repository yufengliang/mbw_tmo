# Initialization

# Yufeng Liang, August 2016 LBNL

# Import Modules
import numpy as np
# from numpy import linalg as la
"""
Important: scipy.linalg seems to be much faster than numpy.linalg
scipy.linalg.det is 6 times faster numpy.linalg.det for 500 x 500 complex matrices
Should I change everything from numpy to scipy instead ?
"""
import scipy as sp
from scipy import linalg as la
import sys
from datetime import datetime
from datetime import timedelta
import pdb
import os

# Parallel Computing

ismpi = True

# This is not robust for every mpi enviroment
if 'SLURM_JOB_ID' in os.environ:

	try:
		from mpi4py import MPI
		print "Running with mpi4py"
	except ImportError:
		ismpi = False
		print "Cannot import mpi4py. Running the serial mode."

else:
	ismpi = False
	print "This is not an mpi environment. Running the serial mode."

if ismpi:
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()
else:
	size = 1
	rank = 0

# Constants

CPLX = np.complex
REAL = np.real
INT = np.int
INT8 = np.dtype('i1')
norm_thr = 0.99
olp_thr = 1e-2
bro_thr = 0.95
xi_thr = 0.05
zero_thr = 1e-30
small_thr = 1e-8
spec_margin = 10 # range for plotting the spectra: margin beyond the range of the eigenstate energies
spec_dener = 0.01 # energy grid step

# Exit the code and report errors
def exit_code(*arg):
	
	# print the error message
	print arg[0]

	if ismpi:
		comm.Abort(0)
	else:
		sys.exit(0)
	
# Record the current time
def timeit():
	nowtime = str(datetime.now().time())
	print nowtime
