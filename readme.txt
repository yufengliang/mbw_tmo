The code used to produce the many-body spectra for transition metal oxides as in
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.096402

The single-particle wavefunctions and energies are imported from Quantum Espresso.
The code relies on the projwfc.x that expands solid-state wavefunctions in terms of atomic wavefunctions,
and make use of these wavefunctions to calculate the overlap matrix, \xi_{ij}.

mbxas_fast.py is the main code to obtain the \xi matrix, and perform the breadth-first search
for the many-body states that are relevant for producing x-ray spectra.

This is an ad-hoc solution to produce many-body spectra for transition metal oxides,
where the orbitals involved in the near-edge transitons can be perfectly described
by local atomic basis sets. Therefore, the many-body spectra are only limited to a 
few eVs above onset.

Please refer to 
https://github.com/yufengliang/mbxaspy
for a complete and generic simulation package for many-body x-ray spectra
