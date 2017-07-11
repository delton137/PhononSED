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
module eig_project
    use main_vars
    implicit none

 contains

!-----------------------------------------------------------------------
!----------------- project velocities onto eigenvector k --------------
!-----------------------------------------------------------------------
subroutine eigen_projection_and_SED(k)
 implicit none
 real(8), dimension(Natoms, 3), intent(in) :: k !!Eigenvector to project onto
 real(8) :: part1
 integer :: length, BlockSize, PointsAvailable


 !-------------- projection ------------------------------
 do t = 1, Ntimesteps
     do j = 1, Nunitcells
         do i = 1, AtomsPerUnitCell
             ia = i + (j-1)*AtomsPerUnitCell
             r(j, :) = (/ 0, 0, 0 /)
             part1 = MassPrefac(i)*dot_product(velocities(t, ia, :), k(ia, :))
             qdot(t, ia, :) = part1*exp(  dcmplx(0, 1)*dot_product(k(ia, :), r(j, :))  )
         enddo
     enddo
 enddo

 !write(*,*) qdot(1:5, 1, 1)


!----------------- calculate Spectral Energy Density -----
 SED = 0.0

 length = Ntimesteps!size(qdot(:,1,1))

 do ia = 1, Natoms
     do ix = 1, 3
         call calc_DFT_squared(qdot(:,ia,ix), oneSED, length)
         SED = SED + oneSED
     enddo
 enddo

 SED = 2d0*(1d0/(2d0*3.14159d0))*SED !!*(1d0/(length*timestep))

!------------ Construct frequencies if they haven't already -------

if (.not.(allocated(spectrum_freqs))) then
    allocate(spectrum_freqs(Ntimesteps))

    do i = 0, Ntimesteps-1
      spectrum_freqs(i+1) = dble(i)/(timestep*2d0*length)
    enddo

    spectrum_freqs = spectrum_freqs/(ps2s*Cspeed) !convert to 1/cm

    allocate(freqs_smoothed(NPointsOut))

    BlockSize = floor(real(length/NPointsOut))

    do i = 1, NPointsOut
       freqs_smoothed(i) = sum(spectrum_freqs((i-1)*BlockSize+1:i*BlockSize) ) /BlockSize
    enddo
endif

!------------ figure out variables for smoothing --------------
 if (.not.(allocated(SED_smoothed))) then
    MinFreqOut =  spectrum_freqs(2) !smallest possible frequency

    PointsAvailable = floor(MaxFreqOut/MinFreqOut) !max number possible

    if (PointsAvailable .lt. NPointsOut) NPointsOut = PointsAvailable

    allocate(SED_smoothed(NPointsOut))

    BlockSize = floor(real(Length)/real(NPointsOut))
 endif

 !------- block averaging / smoothing --------------------------
 do i = 1, NPointsOut
   SED_smoothed(i) = sum(SED((i-1)*BlockSize+1:i*BlockSize) ) /BlockSize
 enddo


end subroutine eigen_projection_and_SED



!------------------------------------------------------------------------
!-------  Compute DFT (discrete Fourier transform) using four1.f
!------------------------------------------------------------------------
subroutine calc_DFT_squared(input, output, tread)
 Implicit none
 integer, intent(in) :: tread
 double complex, dimension(tread), intent(in) :: input
 double precision, dimension(tread), intent(out) ::  output
 complex, dimension(:), allocatable :: transformed
 integer :: trun

 !zero padding
 trun = 2**(    floor( dlog(  dble(tread)  )/dlog(2d0)  )  + 1 )

 if (.not. allocated(transformed)) then
 	allocate(transformed(0:2*trun-1))
 elseif (.not. (2*trun .eq. size(transformed))) then
	deallocate(transformed)
	allocate(transformed(0:2*trun-1))
 endif

 transformed = 0
 transformed(0:tread-1) = input

 call four1(transformed,2*trun,1)

 output=dble(transformed(0:tread-1)*CONJG(transformed(0:tread-1)))/(2d0*trun)  !!**2  !! normalization doesn't seem necessary here

end subroutine calc_DFT_squared

end module eig_project
