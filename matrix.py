# matrix.py

# Useful matrix operations 

# Yufeng Liang, August 2016 LBNL

from init import *

# qr decomposition for finding the orthogonal complement

# swap row
def swap_row(A, r1, r2):
	row_tmp = A[r1].copy(); A[r1] = A[r2].copy(); A[r2] = row_tmp.copy()

# row-based 
def find_oc_vector(A, *arg):

	N = A.shape[0]

	if len(arg) >= 1:
		# switch this row with the last one
		r = arg[0]
		swap_row(A, r, N - 1)

	# perform QR decomposition to A to find its orthogonal complement to the first N - 1 row vectors
	Q, R = la.qr(A.transpose())
	
	# The magnitude of the oc_vector is set to be the same as the product 
	# of the norm of the first N - 1 orthogonalized vectors
	oc_vector = abs(sp.prod(sp.diagonal(R[0 : N, 0 : N]))) * Q[:, N - 1]

	if len(arg) >= 1:
		# switch it back
		swap_row(A, r, N - 1)

	return oc_vector

