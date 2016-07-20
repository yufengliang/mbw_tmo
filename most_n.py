#!/usr/bin/python

import sys

# Find 2 * n most important states and output their indices
# 2 = spin 2

infile = "pdos_digitized.out"
n = 432

if len(sys.argv) == 3:
	infile = sys.argv[0]
	n = int(sys.argv[1])
elif len(sys.argv) != 1:
	print "usage: most_n.py"
	print "usage: most_n.py infile n[# of states]"
	sys.exit()
 
fout = open(infile, "r")

phi2 = []
ener = []
for line in fout:
	phi2.append( float( line.split()[1] ) )
	ener.append( float( line.split()[0] ) )

# This needs to be optimized
ind = sorted(range(len(phi2)), key = lambda k : phi2[k], reverse = True)

for i in range(len(phi2)):
	print ind[i] + 1, phi2[ind[i]], ener[ind[i]]

fout.close()
