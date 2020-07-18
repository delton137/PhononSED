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

*fitter.py* is an example Python code which calculates phonon lifetimes by fitting Lorentzians to the SED data. A newer version has been developed in Matlab.

### Attribution 

If you use this code, please cite this work:

G. Kumar, F. G. VanGessel, D. C. Elton, and P. W. Chung. “[Phonon Lifetimes and Thermal Conductivity of the Molecular Crystal α-RDX](https://www.cambridge.org/core/journals/mrs-advances/article/phonon-lifetimes-and-thermal-conductivity-of-the-molecular-crystal-rdx/14B1FC4424D8C4A659589DC535DBB5A7)”, *MRS Advances*, **4**, 2191 (2019) 

(an arXiv preprint of the work is available [here](https://arxiv.org/abs/1904.12038).)

You may also wish to read and cite: 

F. G. VanGessel, G. Kumar, D. C. Elton, and P. W. Chung, “[A Phonon Boltzmann Study of Microscale Thermal Transport in α-RDX Cook-Off](https://arxiv.org/abs/1808.08295)”, *Proceedings of the 16th International Detonation Symposium*, Cambridge MD, USA, July 2018. 

### Further references
* J. M. Larkin, Ph.D. thesis, Carnegie Mellon University, 2013
* Larkin, et al., *Phys. Rev. B* **81**, 081411(R) (2010)
* A. J. H. McGaughey and M. Kaviany, *Phys. Rev. B* **69**, 094303 (2004).
* J. E. Turney, E. S. Landry, A. J. H. McGaughey, and C. H. Amon, *Phys. Rev. B* **79**, 064301 (2009).
