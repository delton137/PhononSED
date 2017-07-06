!-------------------------------------------------------------------------------------
! Dan's timing module
!-------------------------------------------------------------------------------------
! Copyright (c) 2015 Daniel C. Elton
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
module dans_timer
Implicit None
type timer
	private
	character(len=2000) :: name
	logical :: STATE=.false.
	logical :: USED=.false.
	integer :: nCalls=0
	double precision :: seconds=0d0
	double precision :: turnontime=0d0
end type timer

integer,parameter :: MAX_NUM_TIMERS=50
type(timer), save :: timers(MAX_NUM_TIMERS)

 contains

!-------------------------------------------------------------------
!--------------------- start/create a timer -----------------------
!-------------------------------------------------------------------
Subroutine start_timer(name)
Implicit None
 character(len=*), intent(in) :: name
 integer :: i

  i = find_timer(name)

  !if timer is not found, create it
  if (i .eq. -1) then
	do i = 1, MAX_NUM_TIMERS
		if (timers(i)%USED .eqv. .false.) then
			timers(i)%name = trim(name)
			timers(i)%USED = .true.
			exit
		endif
	enddo
  endif

  if (timers(i)%STATE .eqv. .true.) then
	write(*,*) "timer_mod: WARNING: There is a mistake in the timing scheme. ", name, " timer is alread on!"
  else
	timers(i)%STATE = .true.
	timers(i)%turnontime = get_CPU_time()
  endif

EndSubroutine start_timer

!-------------------------------------------------------------------
!------------------ stop a timer ----------------------------------
!-------------------------------------------------------------------
Subroutine stop_timer(name)
Implicit None
 character(len=*), intent(in) :: name
 integer :: i

  i = find_timer(name)

 if (i .eq. -1) then
	write(*,*) "timer_mod: ERROR: cannot find timer:  ", name
 endif

 if (timers(i)%STATE .eqv. .false.) then
	write(*,*) "timer_mod: WARNING There is a mistake in the timing scheme. ", name, " timer is alread off!"
 else

	timers(i)%STATE = .false.
	timers(i)%seconds = timers(i)%seconds + get_CPU_time() - timers(i)%turnontime
 endif

EndSubroutine stop_timer


!-------------------------------------------------------------------
!------------------ Print timing report ---------------------------
!-------------------------------------------------------------------
Subroutine print_timing_report(iun)
 Implicit None
 integer, intent(in) :: iun
 double precision :: seconds
 integer :: i

 write(iun,*) "#-----------  timing report ----------------------------"
 do i = 1, MAX_NUM_TIMERS
	if (timers(i)%USED .eqv. .true.) then
		!stop timer if its on
		if (timers(i)%STATE .eqv. .true.) then
			timers(i)%STATE = .false.
			timers(i)%seconds = timers(i)%seconds + get_CPU_time() - timers(i)%turnontime
		endif

	seconds = timers(i)%seconds
	write(iun,'(a29, a1, i5, a7, i3, a9, f16.4, a8)') trim(timers(i)%name),":", &
			int(real(seconds)/3600), " hours ",  &
			int(mod(seconds,3600d0)/60), " minutes ", &
			real(mod(seconds,60d0)), " seconds"
	endif
 enddo
! write(iun,*) "#-----------  end timing report ----------------------"

EndSubroutine print_timing_report


!-------------------------------------------------------------------
!------------------ get time from a timer -------------------------
!-------------------------------------------------------------------
Subroutine get_time(name,seconds)
Implicit None
 character(len=*), intent(in)  :: name
 double precision, intent(inout) :: seconds
 integer :: i

  i = find_timer(name)

  seconds = timers(i)%seconds

EndSubroutine get_time

!-------------------------------------------------------------------
!------------------ find timer index ------------------------------
!-------------------------------------------------------------------
Integer function find_timer(name)
	Implicit None
	character(len=*), intent(in) :: name
	integer j

	Do j = 1, MAX_NUM_TIMERS
		if (trim(timers(j)%name) .eq. trim(name)) then
			find_timer = j
			return
		endif
	EndDo

	find_timer = -1

EndFunction


!-------------------------------------------------------------------
!----------------portable CPU_time() -------------------------------
!-------------------------------------------------------------------
Double precision function get_CPU_time()
	Implicit none
!	double precision :: seconds
	integer, dimension(8) :: values

#ifdef mpi
	get_CPU_time = MPI_Wtime()
#else
	!call CPU_time(seconds)       !very compiler/system dependent
    call date_and_time(values=values)  !more portable intrinsic
 	get_CPU_time = values(3)*24*3600.d0 + values(5)*3600.d0 + values(6)*60.d0 + values(7) + values(8)*0.001d0
#endif

endfunction get_CPU_time

EndModule dans_timer
