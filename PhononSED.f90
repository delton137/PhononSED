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
 use mpi
 use main_vars
 use lun_management
 use dans_timer
 use eig_project
 use InputOutput
 integer :: k

 !---------------------- Start MPI and find number of nodes and pid --------------
 !call MPI_Init(ierr)
 if (.not. (ierr .eq. 0)) write(*,*) "WARNING: MPI did not intialize correctly."
 call MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
 write(*, '(a,i4,a)') "running on ", Nnodes, " nodes"

 if(pid .eq. 0)  call start_timer("total")

 call read_input_file

 !read in k-vectors we will be working with
 write(*, *) "reading eigenvector file... "
 call read_eigvector_file

 !read in velocities data
 write(*, *) "reading velocities file... "
 call read_LAAMPS_files

 !calculate/allocate variables
 call calculate_frequencies_and_smoothing

 !calculate equilibrium unit cell coordinates
 call r_unit_cell_coords

 !Calculate SEDs for different eigenvectors
 do ik = 1, Nk
     do ie = 1, Neig
         write(*, '(a,i4,a,i4,a,i4,a,i4)') "Doing eigenvector", ieig, " of ", Neig, " and k vector", ik, " of ", Nk
         call eigen_projection_and_SED(eig_vecs(ik, ie, :, :), all_SED_smoothed(ik, ie, :))
     enddo
 enddo

 !print SED
 if(pid .eq. 0)  call print_SED()

 !print timing report
 if(pid .eq. 0)  call print_timing_report(6)

end program PhononSED
