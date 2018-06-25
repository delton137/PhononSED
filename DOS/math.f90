module math

contains

!-----------------------------------------------------------------------
subroutine calc_corr_function(input, output)
   Implicit none
   double precision, dimension(:), intent(in) :: input
   double precision, dimension(:), intent(out) :: output
   complex, dimension(:), allocatable :: transformed
   integer :: trun, i, tread

   tread = size(input)
   trun = 2**(    floor( dlog(  dble(tread)  )/dlog(2d0)  )  + 1	  )

   if (.not.(size(input) .eq. size(output))) then
   write(*,*) "ERROR: Corr input and output arrays not same size!!"
   stop
   endif
   if (.not. allocated(transformed)) then
     allocate(transformed(0:2*trun-1))
   elseif (.not. (2*trun .eq. size(transformed))) then
   deallocate(transformed)
   allocate(transformed(0:2*trun-1))
   endif

   transformed = 0
   transformed(0:tread-1) = cmplx(input)

   call four1(transformed , 2*trun,1)
   transformed = transformed*conjg(transformed)
   call four1(transformed , 2*trun,-1)


   output = dble(transformed(0:tread-1))!/(2*trun)

   !normalization 1
   do i = 1, tread-1
     output(i) = output(i)/(tread-i)
   enddo

end subroutine calc_corr_function



!------------------------------------------------------------------------
!----------------  Compute DFT four1.f ---------------------------------
!------------------------------------------------------------------------
subroutine calc_DFT(input, output, freqs, timestep, tread)
 Implicit none
 integer, intent(in) :: tread
 double complex, dimension(tread), intent(in) :: input
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
 transformed(0:tread-1) = input

 call four1(transformed,2*trun,1)

 output = dble(transformed(0:tread-1))/(2d0*trun)

 do i = 0, tread-1
	freqs(i+1) = dble(i)/(timestep*2d0*trun)
 enddo

end subroutine calc_DFT



!------------------------------------------------
!----- block averaging ------------------------
!------------------------------------------------
function block_average(input, N)
 integer :: N
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


end module math
