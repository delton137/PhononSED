!-------------------------------------------------------------------------
!-- PhononSED : Phonon Spectral Energy Density calculator
!-- Ref: Jason M. Larkin "Vibrational Mode Properties of Disordered Solids
!-- from High-Performance Atomistic Simulations and Calculations", Ph.D. thesis,
!-- Carnegie Melon University, 2013
!-------------------------------------------------------------------------
! Copyright (c) 2017 Daniel C. Elton
!
! This software is licensed under The MIT License (MIT)
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-------------------------------------------------------------------------------------

Program PhononSED
 use main_vars
 use lun_management
 use dans_timer
 use eig_projec
 use InputOutput

 call start_timer("total")

 !read input file
 call read_input_file

 !read in k-vectors we will be working with
 call read_eigvector_file

 !read in velocities data
 call read_velocities_file(lunvel)

 !test project
 call eigen_projection(eig_vecs(1,:,:))


 !main loop

 !calculate SED

 !print SED

 !end main loop

 call print_timing_report(6)

end program PhononSED
