!-------------------------------------------------------------------------
!-- PhononSED : Phonon Spectral Energy Density calculator
!-- Ref: Jason M. Larkin "Vibrational Mode Properties of Disordered Solids
!-- from High-Performance Atomistic Simulations and Calculations", Ph.D. thesis,
!-- Carnegie Melon University, 2013
!-------------------------------------------------------------------------
! Copyright (c) 2017-2018 Daniel C. Elton
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
#ifdef parallel
 use mpi
#endif
 use main_vars
 use lun_management
 use dans_timer
 use eig_project
 use InputOutput
 implicit none
 integer :: k

 !---------------------- Start MPI and find number of nodes and pid --------------
#ifdef parallel
 call MPI_Init(ierr)
 if (.not. (ierr .eq. 0)) write(*,*) "WARNING: MPI did not intialize correctly."
 call MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
 if (pid .eq. 0) write(*, '(a,i4,a)') "running on ", Nnodes, " nodes"
#endif

 if(pid .eq. 0)  call start_timer("total")

 call read_input_files

 if(pid .eq. 0)  call stop_timer("reading files")

 !calculate/allocate variables
 call calculate_frequencies_and_smoothing

 !Calculate SEDs for different eigenvectors


#ifdef parallel
  Nbatches = Neig/Nnodes

  if (.not. (pid .eq. 0)) then
    allocate(a_eig_vec(3, Natoms))
    allocate(a_SED_smoothed(NPointsOut))
  endif

  do ik = 1, Nk
    if (pid .eq. 0) write(*,*) "doing k vector", ik, "of", Nk
    do bat = 0, Nbatches - 1
       if (pid .eq. 0) then
         !-------- masternode send out jobs -----------------------------------
          do i = 1, Nnodes - 1
              ie = bat*Nnodes + i
              write(*,*)  "proccessor ", i, " doing eigenvector", ie, " of ", Neig
              Call MPI_Send(eig_vecs(ik, ie, :, :), 3*Natoms, MPI_DOUBLE_COMPLEX, i, 0, MPI_COMM_WORLD, ierr)
          enddo

          !-------- masternode job ---------------------------------------------
          ie = bat*Nnodes + Nnodes

          write(*, *)   "proccessor ", pid, " doing eigenvector", ie, " of ", Neig

          call eigen_projection_and_SED(eig_vecs(ik, ie, :, :), all_SED_smoothed(ik, ie, :), ik)

          !-------- masternode recieve work done  ------------------------------
           do i = 1, Nnodes - 1
              ie = bat*Nnodes + i

              call MPI_Recv(all_SED_smoothed(ik, ie, :), NPointsOut, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)

          enddo
      else
         !--------- slavenode part ---------------------------------------------
         call MPI_Recv(a_eig_vec, 3*Natoms,  MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, status2, ierr)

         call eigen_projection_and_SED(a_eig_vec, a_SED_smoothed, ik)

         call MPI_Send(a_SED_smoothed, NPointsOut, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

      endif !pid .eq. 0

      !Call MPI_Barrier(MPI_COMM_WORLD, ierr)

    enddo !bat =0, Nbatches - 1
 enddo !ik = 1, Nk

#else
!-----------------------------------------------------------------------------------------------
!--------------------------------------- serial version ---------------------------------------
!-----------------------------------------------------------------------------------------------
!Calculate SEDs for different eigenvectors
 do ik = 1, Nk
   do ie = 1, Neig
       write(*, '(a,i4,a,i4,a,i4,a,i4)')  "Doing k vector", ik, " of ", Nk, " and eigenvector", ie, " of ", Neig

       if (BTEMD) then
           call BTE_MD(eig_vecs(ik, ie, :, :), all_SED_smoothed(ik, ie, :), ik, ie, all_corr_fns(ik, ie, :))
       else
           call eigen_projection_and_SED(eig_vecs(ik, ie, :, :), all_SED_smoothed(ik, ie, :), ik)
       endif
   enddo
 enddo


#endif

 !print SED
 if(pid .eq. 0)  then
    call print_SED()

    if (BTEMD) call print_corr_fns()

    call print_timing_report(6)

  endif

  call io_close_all() !close ALL open file units

#ifdef parallel
 call MPI_Barrier(MPI_COMM_WORLD, ierr)
 call MPI_Finalize(ierr)
#endif

end program PhononSED
