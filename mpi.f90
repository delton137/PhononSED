!----------------------------------------------------------------------------------------
! This is a module of stub functions and variables which
! allow MPI calls to be executed even if the code is compiled with a serial (traditional) compiler
!
! MPI_Send and MPI_Recv may be called with various objects (arrays, scalars, etc) in the 'input' argument
! The code is written assuming that the 'input' argument is an array, but this may not be convenient
! The use of the 'unlimited' polymorphism feature of F2003 does not appear to solve this problem
!
! Obviously this way of allowing serial code compilation may be viewed as inelegant, but it does allow
! For consistent use of MPI calls (ie, for timing, getting the number of nodes (=1), etc) during serial code execution
!
! Alternatively, this module may be removed and preprocessor declarations (ifdef mpi .. endif) may be placed around
! around all MPI calls instead.
!
!-------------------------------------------------------------------------------------
! Copyright (c) 2014-2015 Daniel C. Elton
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
module mpi
Integer, parameter :: MPI_STATUS_SIZE=1
Integer, parameter :: MPI_COMM_WORLD=1
Integer, parameter :: MPI_DOUBLE_PRECISION=8

 contains


subroutine MPI_Init(ierr)

end subroutine MPI_Init



subroutine MPI_Finalize(ierr)

end subroutine MPI_Finalize




subroutine MPI_Send(input, counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr)
 integer ::  counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr
 double precision, intent(inout), dimension(counti) :: input
 !class(*) :: input
end subroutine MPI_Send


subroutine MPI_Recv(input, counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, status2, ierr)
 integer  ::  status2(MPI_STATUS_SIZE)
 integer   ::  counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr
 double precision, intent(inout), dimension(counti) :: input
 !class(*)  :: input
end subroutine MPI_Recv



subroutine MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
	Nnodes = 1
end subroutine MPI_Comm_size




subroutine MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
	integer, intent(out) :: pid
	pid = 0
end subroutine MPI_Comm_rank




subroutine MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine



function MPI_Wtime()
	double precision ::  MPI_Wtime
        call cpu_time(MPI_Wtime)
end function


end module
