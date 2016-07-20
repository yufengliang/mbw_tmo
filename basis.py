#!/usr/bin/python

import numpy as np
import sys

fmost = open("most_n.out", "r")

ind = [];
phi2 = [];
ener = [];

for l in fmost:
	ind.append(int(l.split()[0]))
	phi2.append(float(l.split()[1]))
	ener.append(float(l.split()[2]))
ns_tot = len(list(fmost))
fmost.close()

# number of states used to construct the matrix 
ns = 864
nspin = 2
spin_div = 1122 # record the state that divide spin up and down state: 0 .. spin_div - 1 (down), spin_div, 2 * spin_div (up) 
norb = 432 * 2 # number of localized orbitals
states = [[0j for j in range(norb)] for i in range(ns)]

print "Starting processing ..."

fatm = open("atomic_proj_digitized.out", "r")

iorb = 0
for l in fatm:
	iorb += 1
	lsp = l.split()
	for i in range(ns):
		states[i][iorb] = float(lsp[2 * (ind[i] - 1)]) + float(lsp[2 * (ind[i] - 1) + 1]) * 1j

print "Finish reading states"
print 'Number of local oribitals (non-spin-degenerate): {0:5d}'.format(iorb)
norb = iorb

# It only takes 1s or so to read 864 most important states
fatm.close()

# convert it into a python array
states_arr = np.array(states)

# check norm
#print "Checking norm of each state ... "
#norm_threshold = 0.99
#print "Output state with a norm < {0:6f}".format(norm_threshold)
#for i in range(norb * 2):
#	snorm = np.linalg.norm(states_arr[i])
#	#states_arr[i] = states_arr[i] / snorm
#	if snorm < norm_threshold:
#		print i + 1, snorm, ener[i]
# 

m = 864
# check orthogonality, looking at the upper trigonal
print "state i  ener[i]  state j energy[j] overlap"
for i in range(m):
	for j in range(i + 1, m):
		# orbital overlap
		olp = np.inner(states_arr[i].conjugate(), states_arr[j])
		# spin overlap
		solp = 1
		if nspin == 2:
			if (ind[i] - 1) / spin_div != (ind[j] - 1) / spin_div:
				solp = 0
		olp = abs(olp) ** 2 * solp
		if olp > 1e-4:
			print i, ener[i], j, ener[j], olp
		
