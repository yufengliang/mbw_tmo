# mbxas.py -- main

# Obtaining the many-body matrix elements for XAS calculations

# An even faster version.

# We only need to calculate the Nelect x Nelect determinant, ONCE !!!

# Yufeng Liang, Sep 2016 LBNL

from init import *
from spec_gen import generate_spectrum, generate_spectrum_Af
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
#det_thr=float(input_var('det_thr', 1e-3))
I_thr=float(input_var('I_thr', 5e-5))
#throw_away_thr=float(input_var('throw_away_thr', det_thr))
I_dump_thr=float(input_var('I_dump_thr', I_thr))
use_advanced_qr=input_var('use_advanced_qr', True)
output_stick=input_var('output_stick', False)
only_sp_and_bse=input_var('only_sp_and_bse', False)
bse_strength=float(input_var('bse_strength', 1.0))
nbnd_bse=int(input_var('nbnd_bse', nbnd_i))

# Determine if advanced qr decompositions (qr_insert, qr_delete ...) will be used
if use_advanced_qr and advanced_qr:
	if only_do_maxfn:
		do_advanced_qr = False
		if rank == 0:
			print("We will NOT do advanced qr decompositions because only_do_maxfn = True. ")
	else:
		do_advanced_qr = True
		if rank == 0:
			print("We will do advanced qr decompositions. ")
else:
	do_advanced_qr = False
	if rank == 0:
		print("We will NOT do advanced qr decompositions. ")


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
	print("e(VBM) gs = ", e_vbm_gs)
	print( "e(CBM) gs = ", e_cbm_gs)

ener_ox = [info_ox[i][2] for i in range(len(info_ox))]
ind_ox = sorted(range(len(ener_ox)), key = lambda k : ener_ox[k])

e_vbm_ox = ener_ox[ind_ox[nelect - 1]]
e_cbm_ox = ener_ox[ind_ox[nelect]]

if rank == 0:
	print("e(VBM) ox = ", e_vbm_ox)
	print("e(CBM) ox = ", e_cbm_ox)

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

	# Non-interacting spectrum: Is this alignment OK ?
	I0[spinc].update({str(j): sp.array([float(ener_gs[c] - e_cbm_ox), abs(wc) ** 2])})

if rank == 0:

	## Single-particle spectrum
	spec0, os_sum_gs = generate_spectrum(If = I0, enerlo = enerlo, enerhi = enerhi, \
	spec_dener = (enerhi - enerlo) / nener, sigma = sigma, eshift = eshift)

	fspec0 = "spec0.dat"

	print("Non-interacting spectra done. ")

	print("Output to ", fspec0)
	
	np.savetxt(fspec0, spec0, delimiter = ' ')

	## BSE spectrum

	# Construct the BSE Hamitonian
	print("Constructing the BSE Hamiltonian ... ")

	ener_mat = sp.matrix(sp.diag(ener_ox[0 : nbnd_f]))
	# xi is not perfectly unitary !!! xi.conjugate().transpose() * xi is not a perfect identity !!!
	# The good news is that H_BSE is almost Hermitian

	# Should we do qr decomposition of xi ? Looks like we shouldn't !
	xiq, xir = la.qr(xi)
	H_BSE_full = sp.matrix(xi.conjugate().transpose()) * ener_mat * sp.matrix(xi)

	# Cancel out the small non-Hermitian part
	H_BSE_full = (H_BSE_full + H_BSE_full.conjugate()) / 2.0

	H_BSE = sp.zeros([nbnd_bse - nelect, nbnd_bse - nelect])
	H_BSE = H_BSE + 1j * H_BSE

	# In the initial-state basis set
	for j1 in range(nelect, nbnd_bse):

		if ind_gs[j1] > nbnd_i: continue

		for j2 in range(nelect, nbnd_bse):

			if ind_gs[j2] > nbnd_i: continue

			# H_bse = H_i + bse_strength * (H_f - H_i)
			H_BSE[j1 - nelect, j2 - nelect] = bse_strength * complex(H_BSE_full[ind_gs[j1], ind_gs[j2]]) \
						        + (1 - bse_strength) * int(j1 == j2) * ener_gs[ind_gs[j1]]

	# Solve the BSE: use only the empty orbitals
	print("Solving the BSE ... ")

	E_BSE, Af_BSE = la.eig(H_BSE)

	#
	I_BSE = [{} for i in range(0, nspin)]

	for j1 in range(nelect, nbnd_bse):
		
		if ind_gs[j1] > nbnd_i: continue

		# define the spin channel of the photoelectron in the initial state
		# Now j1 is in energy order
		spinc1 = 0
		if nspin == 2: spinc1 = spins_gs[j1]

		Aj1 = 0
		
		for j2 in range(nelect, nbnd_bse):

			if ind_gs[j2] > nbnd_i: continue

			#spinc2 = 0
			#if nspin == 2: spinc2 = spins_gs[j2]
			#if spinc1 != spinc2: continue

			Aj1 += states_gs[ind_gs[j2], pol - 1] * Af_BSE[j2 - nelect, j1 - nelect].conjugate()

		I_BSE[spinc1].update({str(j1):sp.array([E_BSE[j1 - nelect] - e_cbm_ox, abs(Aj1) ** 2])})

	spec_bse, os_sum_bse = generate_spectrum(If = I_BSE, enerlo = enerlo, enerhi = enerhi, \
	spec_dener = (enerhi - enerlo) / nener, sigma = sigma, eshift = eshift)

	fspec_bse = "spec_bse.dat"

	print("BSE spectrum done. ")

	print("Output to ", fspec_bse)
	
	np.savetxt(fspec_bse, spec_bse, delimiter = ' ')


if only_sp_and_bse: exit_code("Only want single-particle and BSE spectra. Done. ")

# Construct the xi matrix: now we need all the coefficients for the fast algorithm
xi_mat = np.zeros([nbnd_f, nelect + 1], CPLX)

# Construct the valence overlap matrix
for i in range(nbnd_f): # i goes over the final occupied orbitals

	ind_i = ind_ox[i]

	if ind_i >= nbnd_f: continue
	
	xi_mat[i, 0 : nelect] = xi[ind_i, ind_gs[0 : nelect]]


iter_rank = 0

nIf = 0 # number of If calculated
#nsIf = 0 # number of significant If

# loop over initial-state spin
for ispin in range(0, nspin):

#	pdb.set_trace()

	# set the photoelectron column (only loop over spins)
	xi_mat[:, nelect] = xi_c[:, ispin].copy()

	"""
	Let's first estimate which rows (fn) are important for 
	a given empty state (c = ind_gs[j]) in the initial state.
	We do this by calculating the norm of each line [i, (0 : nelect, c)]

	"""
	
	# Calculate the mother determinant
	det_mom = la.det(xi_mat[0 : nelect + 1, :])
	xi_mat_tmp = xi_mat.copy()

	# if the mother determinant is too small, then refine the last row vector
	if abs(det_mom) < small_thr:

		xi_mat_q, xi_mat_r = la.qr(xi_mat.transpose())
		xi_mat_tmp[nelect] = xi_mat_q[:, nelect].transpose().copy()
		det_mom = la.det(xi_mat_tmp[0 : nelect + 1, :])

	# Construct the zeta-matrix
	xi_mat_inv = la.inv(xi_mat_tmp[0 : nelect + 1, :])

	# Now the zeta matrix include the (N + 1)th orbital !
	zeta_mat = sp.matrix(xi_mat[nelect : nbnd_f, :]) * sp.matrix(xi_mat_inv)

	# Check the sparsity of the zeta-matrix
	sparse_thr = 1e-4

	"""
	Record the coordinates and values of all significant matrix elements
	in such an order
	0 0 4 0 1 
	0 6 5 3 2
	0 7 0 0 0
	Then we need to flip our zeta_mat from left to right, and then tranpose, 
	This looks sloppy. Maybe there is a more elegant way
	"""
	zeta_coord = sp.array(sp.where(abs(sp.fliplr(zeta_mat).transpose()) > sparse_thr))
	# Then I need to convert the coordinates like
	# i = j', j = N - i'
	zeta_coord[0], zeta_coord[1] = zeta_coord[1].copy(), zeta_coord[0].copy()
	zeta_coord[1] = sp.ones(len(zeta_coord[1])) * nelect - zeta_coord[1]
	zeta_coord = zeta_coord.astype(int)

	"""
	Carry out breath-first search for nontrivial spectral contributions 
	from all possible configurations. 

	For each configuration that has been obtained and stored in the queue, the search procedure 
	matches one significant matrix element in the zeta matrix and calculates the contributions of 
	this contribution to a child configuration. The matrix element is matched in an ascending row order
	so as to avoid double-counting, which means the new matrix element can only be chosen
	from the row vectors below the rows that form the parent configuration. 

	"""

	## Now produce the f(1) configurations as the parent configurations
	
	# Format of Af
	# {'c1 v2 c2 v3 c3...': [energy, amplitude/intensity]}

	Af = {}
	Af[''] = sp.array([0.0, det_mom])

	# determinant theshold 
	det_thr = sp.sqrt(abs(I_thr))

	for ndepth in range(1, maxfn + 1):

		# Record Af from all possible final-state indices f
		# If is an array of dictionaries (len = nspin)
		Af_new = {}
		
		"""
		The index of a configuration: f_config = 'c1 v2 c2 v3 c3 ...'
		iv1 is omitted because iv1 \equiv nelect (the lowest-energy conduction state)
		"""

		for f_config in Af:

			# Basic information of the parent configuration

			print(f_config)
			if ndepth > 1:
				
				# valence indices
				f_config_v = [int(ind) for ind in f_config.split()][1 :: 2]
				# the minimum of the valence indices
				f_config_minv = min(f_config_v + [nelect])
				# conduction indices
				f_config_c = [int(ind) for ind in f_config.split()][ :: 2]

			else:
			
				f_config_v = []
				f_config_minv = nelect + 1 # the (N + 2)th state
				f_config_c = []

			f_config_c_set = set(f_config_c)

			# energy of the parent configuration
			f_ener = Af[f_config][0]

			## peform the search in descending column order
			# I find doing ascending row order will screw up things !!!

			# Make sure it is in descending column order 	
			low_izeta = bisect.bisect_left(-zeta_coord[1], -(f_config_minv - 1))

			for izeta in range(low_izeta, len(zeta_coord[0])):

				# zeta_coord[0] = 0 corresponds c = nelect, the (N+1)th state 
				new_c = zeta_coord[0][izeta] + nelect

				new_v = zeta_coord[1][izeta]

				# ndepth = 1 is special
				if ndepth == 1 and new_v < nelect: break

				# Make sure c doesn't appear twice in a configuration
				if new_c in f_config_c_set: continue

				# Energy filter
	
				# don't go too deep into the valence band
				if ener_ox[ind_ox[nelect]] - ener_ox[ind_ox[new_v]] + f_ener > e_hi_thr: break
				
				enew = ener_ox[ind_ox[new_c]] - ener_ox[ind_ox[new_v]] + f_ener
				if enew > e_hi_thr: continue

				"""
				Find out the sign for the child configuration

				{     c1, new_v}   ...  {     c1, v_{n-1}}   ...  {     c1, v1}
				...............................................................
				{  new_c, new_v}   ...  {  new_c, v_{n-1}}   ...  {  new_c, v1}
				...............................................................
				{c_{n-1}, new_v}   ...  {c_{n-1},   new_v}   ...  {c_{n-1}, v1}

				sign = (-1) ** i, i is the order new_c locates in c1, c2, ..., c_{n-1}
				"""

				insert_pos = 0

				v_sgn = 1
				if ndepth > 1:

					# Reverse the descending valence index sequence and determine the where
					# to insert the new_v index. 
					insert_pos =  bisect.bisect_left(sp.array(f_config_c), new_c)
					v_sgn *= (-1) ** insert_pos
				
				# Construct the new configuration

				f_config_new = ''

				for idepth in range(0, ndepth):

					# add c_{idepth} and insert new_c
					if insert_pos == idepth:
						f_config_new += str(new_c) + ' '
					else:
						f_config_new += str(f_config_c[idepth - int(idepth > insert_pos)]) + ' '

					# add v_{idepth}
					if idepth < ndepth - 2: 
						f_config_new += str(f_config_v[idepth]) + ' '
					elif idepth == ndepth - 2:
						f_config_new += str(new_v) + ' '

				#if f_config_new == '626 622 657 ':
				#	i_am_here = 1
				### Calculate the contribution to the child configuration
				# = sign x det(n-1) x new_matrix_element
				new_contribution = v_sgn * Af[f_config][1] * zeta_mat[zeta_coord[0][izeta], zeta_coord[1][izeta]]

				# Too small ? throw it away !				
				if abs(new_contribution) < det_thr: continue

				# Check if this configuration already exists
				if f_config_new in Af_new:

					Af_new[f_config_new][1] += new_contribution

				else:

					Af_new[f_config_new] = sp.array([enew, new_contribution])
				
			# end for izeta	

		# end for f_config

		# adapt det_thr
		det_thr /= 1.0

		# gather all Af_new, record the results, and copy Af_new into A_f
		if ismpi:

			comm.barrier()
		# end if ismpi
		else:

			Af = Af_new.copy()
			del Af_new

		# Turn Af into a spectrum
		spec_tmp, os_sum_tmp = generate_spectrum_Af(Af = Af, enerlo = enerlo, enerhi = enerhi, \
		spec_dener = (enerhi - enerlo) / nener, sigma = sigma, eshift = eshift, nspin = nspin, ispin = ispin)

		if ispin == 0 and ndepth == 1:
			spec = spec_tmp.copy()
			os_sum = os_sum_tmp
		else:
			spec[:, ispin + 1] += spec_tmp[:, ispin + 1]
			os_sum += os_sum_tmp

	# end for ndepth

	#iter_rank += 1
	#iter_rank %= size

#print os_sum_gs

print("process ", rank, "done")

# Output
if rank == 0:

	if only_do_maxfn:
		fspec = "spec.only_maxfn_" + str(maxfn) + ".dat"
	else:
		fspec = "spec.maxfn_" + str(maxfn) + ".dat"

	#print "Number of significant A^f's that need calculation: ", len(If_all)

	print("Sum-rule test: ", os_sum / os_sum_gs)

	print("Interacting spectra with ", maxfn, " order of shakeup effects done. ")

	print("Output to ", fspec)
	
	np.savetxt(fspec, spec, delimiter = ' ')
	
	if output_stick:
		sp.save('If_stick', If_all)
		
