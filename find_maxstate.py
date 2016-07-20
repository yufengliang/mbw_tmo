#!/usr/bin/python

# Find n states that have the largest "overlap" with the subspace spanning 
# over the chosen n local orbitals
# The local orbitals are defined in subspace.dat
# the overlap is defined as the summation of |coef|^2

# input the orbitals
fsp = open('subspace.dat', 'r')

orbital = []
orb_ind = []

for line in fsp:
	orbital.append( line.split()[0] )
	orb_ind.append( int(line.split()[1]) )

norb = len(orbital)
orbset = set(orb_ind)
print orbset
fsp.close()

# find the maximally overlapped states
fout = open('pdos_digitized.out', 'r')

n = 2 * norb
# initialization
maxe   = [0 for i in range(n)];  maxolp = list(maxe);
maxind = [int(0) for i in range(n)];

ind = int(0)
for line in fout:
	ind += 1
	lsp =  line.split()
	e = float( lsp[0] )
	ns = int( lsp[2] )
	orb_ind_s = [int(i) for i in lsp[3:3 + ns]]
	orb_olp_s = [float(i) for i in lsp[3 + ns:3 + 2 * ns]]
	olp = 0.0
	for i in range(ns): 
		if orb_ind_s[i] in orbset: olp += orb_olp_s[i]
	for i in range(n):
		if olp > maxolp[i]:
			maxolp[i] = olp
			maxe[i] = e
			maxind[i] = ind
			break
			
fout.close()

print maxe
print maxind
print maxolp



