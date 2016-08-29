# mbxas.py -- main

# Obtaining the many-body matrix elements for XAS calculations

# Yufeng Liang, August 2016 LBNL

from init import *
from spec_gen import generate_spectrum
from matrix import find_oc_vector

# Input subroutine

def input_var(var_name, *arg):

	try:
		val = os.environ[var_name]
		# only do this if a boolean is detected
		if val.lower() in ['true', '.true.']:
			val = True
		elif val.lower() in ['false', '.false.']:
			val = False
		return val

	except KeyError:
		if len(arg) == 0:
			raise
		return arg[0]

# ------------------------------------------------------------------------------------

# Input user-defined variables
# There should be a more automatic way mapping external variables onto python ones
# I don't like this way of input

sigma = float(input_var('SIGMA', 0.2))
enerlo = float(input_var('ELOW'))
enerhi = float(input_var('EHIGH'))
eshift = float(input_var('ESHIFT'))
nener = float(input_var('NENER'))
maxfn = int(input_var('maxfn', 1))
only_do_maxfn = input_var('only_do_maxfn', False)
nbnd_i = int(input_var('nbnd_i'))
nbnd_f = int(input_var('nbnd_f', nbnd_i))
nspin = int(input_var('nspin', 1))
n_local_orbital = int(input_var('n_local_orbital'))
nelect = int(input_var('nelect'))
pol = int(input_var('pol'))
e_lo_thr=float(input_var('e_lo_thr'))
e_hi_thr=float(input_var('e_hi_thr'))
det_thr=float(input_var('det_thr', 1e-3))
throw_away_thr=float(input_var('throw_away_thr', 1e-6))
use_advanced_qr=input('use_advanced_qr', True)

# Determine if advanced qr decompositions (qr_insert, qr_delete ...) will be used
if use_advanced_qr and advanced_qr:
	do_advanced_qr = True
	if rank == 0:
		print "We will do advanced qr decompositions. "
else:
	do_advanced_qr = False
	if rank == 0:
		print "We will NOT do advanced qr decompositions. "

# Given that we've known about everything below
info_gs = np.load('info_gs.npy')
info_ox = np.load('info_ox.npy')
spins_gs = np.load('spins_gs.npy')
spins_ox = np.load('spins_ox.npy')
states_gs = np.load('states_gs.npy')
states_ox = np.load('states_ox.npy')
xi = np.load('xi.npy')

ener_gs = [info_gs[i][2] for i in range(len(info_gs))]
ind_gs = sorted(range(len(ener_gs)), key = lambda k : ener_gs[k])

e_vbm_gs = ener_gs[ind_gs[nelect - 1]]
e_cbm_gs = ener_gs[ind_gs[nelect]]

if rank == 0:
	print "e(VBM) gs = ", e_vbm_gs
	print "e(CBM) gs = ", e_cbm_gs

ener_ox = [info_ox[i][2] for i in range(len(info_ox))]
ind_ox = sorted(range(len(ener_ox)), key = lambda k : ener_ox[k])

e_vbm_ox = ener_ox[ind_ox[nelect - 1]]
e_cbm_ox = ener_ox[ind_ox[nelect]]

if rank == 0:
	print "e(VBM) ox = ", e_vbm_ox
	print "e(CBM) ox = ", e_cbm_ox

# ------------------------------------------------------------------------------------

"""
Turns out you can do this summation by calculating just one determinant

\sum_{c} A^f_c w_c

Only the last column is replaced as c changes. Just sum up all the c columns weighted 
by w_c and place the summation in the last column. Then we don't need to loop over
c anymore ;-)

"""

# Construct the last column which is a sum of all the c columns weighted by wc
xi_c = np.zeros([nbnd_f, nspin], CPLX)

# sum of the initial-state oscillator strengths
# This can be affected by ( ~25% ) by the GS procedure performed on the states !!! 
# os_sum_gs should be equal to nspin. GS can increase that by ~25% !!!
os_sum_gs = 0.0

I0 = [{} for i in range(0, nspin)]

for j in range(nelect, nbnd_i):

	c = ind_gs[j]
	if c > nbnd_i: continue

	# define the spin channel of the photoelectron in the initial state
	spinc = 0
	if nspin == 2: spinc = spins_gs[c]

	wc = states_gs[c, pol - 1]

	# This doesn't help at all
	# if abs(wc) < det_thr: continue

	os_sum_gs += abs(wc) ** 2

	for i in range(nbnd_f):

		f = ind_ox[i]
		if f > nbnd_f: continue

		xi_c[i, spinc] += xi[f, c] * wc

	print spinc, c
	# Non-interacting spectrum
	I0[spinc].update({c: sp.array([float(ener_gs[c] - e_cbm_ox), abs(wc) ** 2])})

if rank == 0:

	spec0, os_sum_gs = generate_spectrum(If = I0, enerlo = enerlo, enerhi = enerhi, \
	spec_dener = (enerhi - enerlo) / nener, sigma = sigma, eshift = eshift)

	fspec0 = "spec0.dat"

	print "Non-interacting spectra done. "

	print "Output to ", fspec0
	
	np.savetxt(fspec0, spec0, delimiter = ' ')

# Construct the (N + 1) * (N + 1) matrix
xi_mat = np.zeros([nelect + 1, nelect + 1], CPLX)

spcount = 0
sp_thr = 0.01
# Construct the valence overlap matrix
for i in range(nelect + 1): # i goes over the final occupied orbitals

	for j in range(nelect + 1): # j goes over the initial occupied orbitals

		xi_mat[i, j] = xi[ind_ox[i], ind_gs[j]]
		
		if abs(xi_mat[i, j]) < sp_thr: # check the sparsity of the matrix
			spcount += 1	

# row: just the first nelect columns, excluding the photoelectron column
row_norm = np.array([la.norm(xi[i, ind_gs[0 : nelect]]) for i in range(nbnd_f)])

# ------------------------------------------------------------------------------------

"""
 parallelization over f(1) only

 structure of the loop

 for ic1 in range(nelect, nbnd_f): # these two are the outer loops. Only parallelize this one

	for iv2 in range(nelect - 1, -1, -1): # the recursive routine begins here
		for ic2 in range(ic1 + 1, nbnd_f):
 			for iv3 in range(iv2 - 1, -1, -1):
				for ic3 in range(ic2 + 1, nbnd_f):
					......

"""

def fgen(ndepth, maxv, minc, efinal, f_config, row_prod, oc_vector):
	
	global iter_rank

	# You may consider put last_v on the argument list
	last_v = maxv + 1
	e_last_v = ener_ox[ind_ox[last_v]]

	# save the current status
	xi_mat_v_tmp = xi_mat[last_v].copy()

	#
	for ic in range(minc, nbnd_f):

		icind = ind_ox[ic]

		if icind >= nbnd_f: continue

		# energy filter: lower bound not apply here
		enew = ener_ox[icind] - e_last_v + efinal

		if enew > e_hi_thr: break

		# norm filter: this could be quite loose
		row_prod_ic = row_prod * row_c_norm[icind]

		if row_prod_ic < det_thr: continue

		if rank == iter_rank:
			# do it !

			if ndepth == 1:
				print rank, ic
			
			# Construct the new configuration
			f_config_ic = f_config + str(ic) + ' '

			# Construct the row vector at ic
			ic_vec = sp.append(xi[icind, ind_gs[0 : nelect]], xi_c[ic, ispin])

			# replace the occupied (maxv + 1) orbital chosen in the last loop
			# with a chosen empty (ic) orbital.
			xi_mat[last_v] = ic_vec.copy()

			# evaluate when only reach the energy threshold
			if enew >= e_lo_thr:

				# Only interested in f^(maxfn) ?
				if not only_do_maxfn or only_do_maxfn and ndepth == maxfn:

					# Here comes our determinant !
	
					# This is the old way
					# Afc_wc = la.det(xi_mat) * wc
					# A new way using the orthogonal complement (I don't care the phase):
					#Intensity = abs(la.det(xi_mat)) ** 2
					Intensity = abs(sp.dot(oc_vector.conjugate(), ic_vec)) ** 2

					if sp.sqrt(Intensity) > throw_away_thr:
			
						# Data structure of If[ispin]
						# {'ic1 iv2 ic2 ...': [energy, Intensity]}
						# Eachi configuration is unique in a spin manifold If[ispin]
						
						If[ispin][f_config_ic] = sp.array([enew, Intensity])
			
						# Checking this or not depends on the parallelization
						#if sp.sqrt(If) > det_thr:
			
							#print len(If), efinal, f_config, estimate, np.abs(Afc_wc)


			# If hasn't reached specified depth, then continue to the next recursion

			if ndepth < maxfn:
	
				for iv in range(maxv, -1, -1):

					ivind = ind_ox[iv]

					if ivind >= nbnd_f: continue

					# Don't go too deep into the valence bands
					if e_cbm_ox - ener_ox[ivind] + enew > e_hi_thr: break

					try:
						row_prod_iv = row_prod_ic / row_c_norm[ivind]

					except ZeroDivisionError:
						# I will need to improve this. But basically this is impossible
						exit_code("fgen: row_prod filter is broken. ")

					# Only interested in f^(maxfn) ?
					if not only_do_maxfn or only_do_maxfn and ndepth == maxfn - 1:
	
						### Here comes another expensive calculation: find the oc_vector for each (ic, iv)
						oc_vector_new = find_oc_vector(xi_mat, iv)
	
					# Construct the new configuration
					f_config_iv = f_config_ic + str(iv) + ' '
					#print f_config_iv

					fgen(ndepth = ndepth + 1, maxv = iv - 1, minc = ic + 1, efinal = enew, \
					f_config = f_config_iv, row_prod = row_prod_iv, oc_vector = oc_vector_new)
		
		# endif rank == iter_rank

		# Only parallelize this
		if ndepth == 1:
			iter_rank += 1
			iter_rank %= size

	xi_mat[last_v] = xi_mat_v_tmp.copy()
	
# end fgen

# Record Af from all possible final-state indices f
# If is an array of dictionaries (len = nspin)
If = [{} for i in range(0, nspin)]

iter_rank = 0

# loop over initial-state spin
for ispin in range(0, nspin):

#	pdb.set_trace()

	# set the photoelectron column (only loop over spins)
	xi_mat[:, nelect] = xi_c[0 : nelect + 1, ispin].copy()

	"""
	Let's first estimate which rows (fn) are important for 
	a given empty state (c = ind_gs[j]) in the initial state.
	We do this by calculating the norm of each line [i, (0 : nelect, c)]

	"""
	
	# perform QR decomposition to xi_mat to find its orthogonal complement to the first N - 1 row vectors
	xi_mat_Q, xi_mat_R = la.qr(xi_mat.transpose())
	
	try:
		oc_vector = abs(sp.prod(sp.diagonal(xi_mat_R))) / abs(xi_mat_R[nelect, nelect]) * xi_mat_Q[:, nelect]
	except ZeroDivisionError:
		print "Orthogonality catastrophe occurs, my friend. "
		oc_vector = xi_mat_Q[:, nelect] * 0.0

	# an array that stores the norm of each row
	# Now add the first nelect column with the c column
	row_c_norm = REAL(np.sqrt(row_norm ** 2 + np.abs(xi_c[:, ispin]) ** 2))

	# This is the product of the norm of the lowest nelect states
	row_prod = np.prod(row_c_norm[ind_ox[0 : nelect]])

	fgen(ndepth = 1, maxv = nelect - 1, minc = nelect, efinal = 0.0, \
	f_config = '', row_prod = row_prod, oc_vector = oc_vector)

#print os_sum_gs

print "process ", rank, "done"

# Generate the spectra

# Do it on each core, which is faster than doing it just on the root

spec, os_sum = generate_spectrum(If = If, enerlo = enerlo, enerhi = enerhi, \
spec_dener = (enerhi - enerlo) / nener, sigma = sigma, eshift = eshift)

# Gather results if mpi
if ismpi:

	comm.barrier()

	# Reduce the spectrum

	# comm.reduce might not be robust for python on every platform
	spec_all = comm.reduce(spec[:, 1 : nspin + 1], op = MPI.SUM)
	os_sum_arr = sp.array([os_sum])
	os_sum_all = comm.reduce(os_sum_arr, op = MPI.SUM)

	if rank == 0:
		spec[:, 1 : nspin + 1] = spec_all.copy()
		os_sum = os_sum_all[0]

	# Gather the stick
	If_gather = comm.gather(If, root = 0) # Don't ever do this in the rank == 0 block

	# Only gather to pid = 0
#	if rank == 0:
#
#		print "Begin to gather If..."
#	
#		If_all = {}
#
#		# This is extremely inefficient !!!
#		for iter_If in If_gather:
#
#			for icommon in set(iter_If) & set(If_all):
#				If_all[icommon][1 : nspin + 1] += iter_If[icommon][1 : nspin + 1]
#
#			for iupdate in set(iter_If) - set(If_all):
#				If_all[iupdate] = iter_If[iupdate].copy()
#		
#		print "Done gather If."

else:
	If_all = If

# Output
if rank == 0:

	#sp.save('If_stick', If_all)
		
	if only_do_maxfn:
		fspec = "spec.only_maxfn_" + str(maxfn) + ".dat"
	else:
		fspec = "spec.maxfn_" + str(maxfn) + ".dat"

	#print "Number of significant A^f's that need calculation: ", len(If_all)

	print "Sum-rule test: ", os_sum / os_sum_gs

	print "Interacting spectra with ", maxfn, " order of shakeup effects done. "

	print "Output to ", fspec
	
	np.savetxt(fspec, spec, delimiter = ' ')
