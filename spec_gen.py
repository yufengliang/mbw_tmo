# spec_gen.py

# Yufeng Liang, July 2016 LBNL

# Generate the spectra

from init import *

# Lorentzian
def lorentzian(ediff, I, sigma):
	return I / np.pi * sigma / (ediff ** 2 + sigma ** 2)

def generate_spectrum(If, enerlo, enerhi, spec_dener, sigma, eshift):

	nspin = len(If)
	ener = np.arange(enerlo, enerhi, spec_dener)
	spec = np.zeros([len(ener), nspin + 1])
	spec[:, 0] = ener

	os_sum = 0

	for ispin in range(0, nspin):

		for iconf in If[ispin]:
		
			ener_f = If[ispin][iconf][0]
		
			I_f = If[ispin][iconf][1]
			
			os_sum += I_f
		
			for s in range(len(ener)):
		
				spec[s, ispin + 1] += lorentzian(ener[s] - (ener_f + eshift), I_f, sigma)
	
	return spec, os_sum
	


