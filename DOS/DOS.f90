!-----------------------------------------------------------------------------------
! code hacked together to calculate velocity-velocity autocorrelation function
! and its Forier transform, called the "density of states"
! Daniel C. Elton, 2015 - 2018
!-----------------------------------------------------------------------------------
program DOS
 use math
 implicit none
 !--------------------------- user specified inputs -----------------------------------
 double precision, parameter :: timestep = 10*0.00002418  !timstep in ps
 integer  :: NumPointsOut = 6000
 integer, parameter :: Natoms = 24
 integer, parameter :: Ntimesteps = 44844
 logical, parameter :: has_comment_line = .false. !whether the xyz file has a comment line (2nd line before coordinates)
 double precision, parameter :: Cspeed=3.00d10 ! cm/s
 double precision, parameter :: ps2s=1d-12 ! 1fs in s
 !--------------------------- end inputs ----------------------------------------------
 integer, parameter :: lunvel = 15
 integer, parameter :: lun = 16
 double precision, dimension(:,:,:), allocatable :: velocities
 integer :: trun, tcor, ix, i, ia,  t, n,  PointsAvailable
 double precision :: omega, IR, magn, avgMag, MinFreqOut, MaxFreqOut
 double precision, dimension(:), allocatable ::  output, DFT, allfreqs
 double precision, dimension(:), allocatable :: spectra_smoothed, allfreqs_smoothed
 double complex, dimension(:), allocatable :: aux1, vcross, ACF

 integer :: NumAvgOver

 open(lunvel, file='HfS2.vel', status='old', action='read')

 allocate(velocities(Ntimesteps, Natoms, 3))

 do t = 1, Ntimesteps
      if (mod(t,1000).eq.0) write(*,*) "reading: ", t, "of", Ntimesteps

      read(lunvel,*)
      if (has_comment_line) read(lun,*) !comment line

      do ia = 1, Natoms
         read(lunvel,*) (velocities(t, ia, ix), ix=1,3)
      enddo
 enddo

 allocate(output(Ntimesteps))! Stores the cross correlation in reciprocal space
 allocate(ACF(Ntimesteps))   ! stores the autocorrelation Function (real space)
! allocate(Vcross(0:Ntimesteps-1))   ! stores the autocorrelation Function (real space)
! allocate(aux1(0:Ntimesteps-1))
 allocate(allfreqs(Ntimesteps))
 allocate(DFT(Ntimesteps))   ! stores the autocorrelation Function (real space)

 ACF=0
!find correlation function
 do ia = 1, Natoms
   do ix = 1, 3
 		call calc_corr_function(velocities(:, ia, ix), output)
		ACF = ACF + output
	enddo
enddo

! Do i = 1, Natoms!
  !  do ix=1,3
  !       ! Call to the direct FFT
    !     aux1=cmplx(velocities(:, i , ix))
  !     call four1(aux1,Ntimesteps,-1)
  !	     !Build the cross correlation in fourier space
  !       Vcross=aux1*conjg(aux1)
	 !      ! Inverse transform CrossV
  !       call four1(Vcross, Ntimesteps, 1)
  !       Vcross=Vcross/real(size(Vcross))
  !       !Store the autocorrelation function
  !       ACF=ACF+Vcross
!     enddo
 !enddo

 ACF = ACF/dble(Natoms)

 !! Save the correlation function to file--------------------
 open(lun, file= "out_vel_vel_corr_function.dat")
 do i = 1, Ntimesteps
  	write(lun,*) i*timestep, real(ACF(i))/real(ACF(1))
 enddo
 close(lun)
 !-----------------------------------------------------------

 call calc_DFT(ACF, DFT, allfreqs, timestep, size(ACF))

 !allfreqs = allfreqs/(ps2s*Cspeed) !convert to cm^-1

 allfreqs = (allfreqs/(ps2s))/(10d0**12) !convert to THz


 MinFreqOut =  allfreqs(2)
 MaxFreqOut = MinFreqOut*Ntimesteps
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)

 if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable

 allocate(allfreqs_smoothed(NumPointsOut))
 allocate(spectra_smoothed(NumPointsOut))

 allfreqs_smoothed = block_average(allfreqs(1:PointsAvailable),  NumPointsOut)
 spectra_smoothed  = block_average(DFT(1:PointsAvailable), NumPointsOut)

 !! Save the DOS to file--------------------
 open(lun, file= "out_DOS.dat", action="write")
 do i = 1, NumPointsOut !PointsAvailable
!  write(lun,'(2g12.4)')  allfreqs(i), DFT(i)
	write(lun,'(2g12.4)')  allfreqs_smoothed(i), spectra_smoothed(i)
 enddo
 close(lun)
 !!----------------------------------------


end program DOS
