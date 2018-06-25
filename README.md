### Phonon SED
The *PhononSED* Fortran-90 code calculates phonon projected Spectral Energy Densities (SEDs) from molecular dynamics atomic velocity data and phonon eigenvectors calculated with GULP.

Compile with
`make`

Run as:
`./PhononSED.x < PhononSED.inp`
`python fitter.py`

To compile a parallel version with MPI use
`make parallel`

To run the parallel version, use
`mpirun -stdin all -np *number_of_processors* ./PhononSED.x < RDXPhononSED.inp`

*fitter.py* is an example python code calculates phonon lifetimes by fitting Lorentzians to the SED data. A newer version has been developed in Matlab.


### References
* J. M. Larkin, Ph.D. thesis, Carnegie Mellon University, 2013
* Larkin, et al., *Phys. Rev. B* **81**, 081411(R) (2010)
* A. J. H. McGaughey and M. Kaviany, *Phys. Rev. B* **69**, 094303 (2004).
* J. E. Turney, E. S. Landry, A. J. H. McGaughey, and C. H. Amon, *Phys. Rev. B* **79**, 064301 (2009).
