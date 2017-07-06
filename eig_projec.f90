!-------------------------------------------------------------------------------------
! eigvec_projection
!-------------------------------------------------------------------------------------
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
module eig_projec
    use main_vars
    implicit none

 contains

!-----------------------------------------------------------------------
!----------------- project velocities onto eigenvector k --------------
!-----------------------------------------------------------------------
subroutine eigen_projection(k)
 implicit none
 real(8), dimension(Natoms, 3), intent(in) :: k !!Eigenvector to project onto
 real(8), dimension(:,:), allocatable :: r
 real(8) :: part1

 allocate(qdot(Ntimesteps, Natoms, 3))
 allocate(r(Nunitcells, 3))

 write(*, *) velocities(1,1, :)

 !projection
 do t = 1, Ntimesteps
     do j = 1, Nunitcells
         do i = 1, AtomsPerUnitCell
             ia = i + (j-1)*AtomsPerUnitCell
             r(j, :) = (/ 0, 0, 0 /)
             part1 = MassPrefac(ia)*dot_product(velocities(t, ia, :), k(ia, :))
             qdot(t, ia, :) = part1*exp(  dcmplx(0, 1)*dot_product(k(ia, :), r(j, :))  )
         enddo
     enddo
 enddo

end subroutine eigen_projection



!-----------------------------------------------------------------------
!----------------- calculate Spectral Energy Density  -----------------
!-----------------------------------------------------------------------
subroutine calculate_SED()
 implicit none





end subroutine calculate_SED


!------------------------------------------------------------------------
!-------  Compute DFT (discrete Fourier transform) using four1.f
!------------------------------------------------------------------------
subroutine calc_DFT(input, output, freqs, timestep, tread)
 Implicit none
 integer, intent(in) :: tread
 double precision, dimension(tread), intent(in) :: input
 double precision, intent(in) :: timestep
 double precision, dimension(tread), intent(out) :: freqs, output
 complex, dimension(:), allocatable :: transformed
 integer :: trun, i

 trun = 2**(    floor( dlog(  dble(tread)  )/dlog(2d0)  )  + 1 )

 if (.not. allocated(transformed)) then
 	allocate(transformed(0:2*trun-1))
 elseif (.not. (2*trun .eq. size(transformed))) then
	deallocate(transformed)
	allocate(transformed(0:2*trun-1))
 endif

 transformed = 0
 transformed(0:tread-1) = cmplx(input)

 call four1(transformed,2*trun,1)

 output = dble(transformed(0:tread-1))/(2d0*trun)

 do i = 0, tread-1
	freqs(i+1) = dble(i)/(timestep*2d0*tread)
 enddo

end subroutine calc_DFT


!------------------------------------------------
!----- block averaging for smoothing -----------
!------------------------------------------------
function block_average(input, N)
 integer, intent(in) :: N  !block size
 double precision, dimension(:), intent(in) :: input
 double precision, dimension(N)  :: block_average
 integer :: i, BlockSize

 BlockSize = floor(real(size(input))/N)

 if (BlockSize .eq. 0) then
	block_average = input
	return
 endif

 do i = 1, N
	block_average(i) = sum( input((i-1)*BlockSize+1:i*BlockSize) ) /BlockSize
 enddo

end function block_average

end module eig_projec
