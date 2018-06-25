### Phonon SED
PhononSED.f90 Fortran-90 code calculates Projected Phonon Spectral Energy Density (SED) from molecular dynamics atomic velocity data and phonon eigenvectors.

fitter.py Python code calculates phonon lifetimes by fitting Lorentzians to the SED data

Compile with
`make`

Run as:
`./PhononSED.x < PhononSED.inp`
`python fitter.py`

To compile a parallel version with MPI used
`make parallel`

To run the parallel version, used
`mpirun -stdin all -np *number_of_processors* ./PhononSED.x < RDXPhononSED.inp`

references:
* J. M. Larkin, Ph.D. thesis, Carnegie Mellon University, 2013
* Larkin, et al., *Phys. Rev. B* **81**, 081411(R) (2010)
* A. J. H. McGaughey and M. Kaviany, *Phys. Rev. B* **69**, 094303 (2004).
* J. E. Turney, E. S. Landry, A. J. H. McGaughey, and C. H. Amon, *Phys. Rev. B* **79**, 064301 (2009).
