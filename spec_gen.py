# spec_gen.py

# Yufeng Liang, July 2016 LBNL

# Generate the spectra

from init import *

# Lorentzian
def lorentzian(ediff, I, sigma):
	return I / np.pi * sigma / (ediff ** 2 + sigma ** 2)

def generate_spectrum(If, nspin, enerlo, enerhi, spec_dener, sigma, eshift):

	ener = np.arange(enerlo, enerhi, spec_dener)
	spec = np.zeros([len(ener), nspin + 1])
	spec[:, 0] = ener

	os_sum = 0

	for iconf in If:
	
		ener_ic = float(If[iconf][0])
	
		if nspin == 2:
		
			os_sum += If[iconf][1] + If[iconf][2]
	
		else:
	
			os_sum += If[iconf][1]
		
		for s in range(len(ener)):
	
			if nspin == 2:
				spec[s, 1] += lorentzian(ener[s] - (ener_ic + eshift), If[iconf][1], sigma)
				spec[s, 2] += lorentzian(ener[s] - (ener_ic + eshift), If[iconf][2], sigma)
			else:
				spec[s, 1] += lorentzian(ener[s] - (ener_ic + eshift), If[iconf][1], sigma)

	return spec, os_sum
	


