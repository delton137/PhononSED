!-----------------------------------------------------------------------------------
! Fortan-90 module for logical unit management
! By using logical unit management, one ensures logical units for reading/writing
! do not conflict with those used by this code or interfaced code
!
! Inspiration for this comes from code by Alberto Garcia and Richard Maine
!-----------------------------------------------------------------------------------
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
Module lun_management
 implicit none
 integer, save 	    :: stdout=6, stderr=0
 integer   	    	:: iostat
 integer, parameter :: min_lun=10, max_lun=99
 integer, parameter :: nunits = max_lun - min_lun+1
 logical, save      :: lun_is_free(min_lun:max_lun)
 logical 	        :: used

 !initially we assume all units between 10 and 99 are free
 data lun_is_free /nunits*.true./

 contains


!---------------------------------------------------
!Looks for a free unit and assigns it to lun
!---------------------------------------------------
subroutine io_assign(lun)
 integer, intent(inout) :: lun

 do lun= min_lun, max_lun
	if (lun_is_free(lun)) then
		inquire(unit=lun, opened=used, iostat=iostat)

		if (iostat .ne. 0) used = .true.

        if (used .eqv. .true.) then
			lun_is_free(lun) = .false.
		else
			return
		endif

	endif
 enddo

 end subroutine io_assign


!---------------------------------------------------
!Assign & open a unit, and backup old files
!---------------------------------------------------
subroutine io_open(lun,filename,APPEND,REPLACE)
 integer, intent(inout) :: lun
 character(len=*), intent(in) :: filename
 logical, optional :: APPEND
 logical, optional :: REPLACE
 logical :: EXISTS

 call io_assign(lun)

 inquire(file=filename, exist=EXISTS)

 if (EXISTS) then
	if (present(APPEND) .and. APPEND .eqv. .true.) then
		open(lun, file=filename, status="old", position="append", action="write")
	else if (present(REPLACE) .and. REPLACE .eqv. .true.) then
		open(lun, file=filename, status="replace", action="write")
	else
		call system('mv '//filename//' '//filename//'_backup.dat')
		open(lun, file=filename, action="write")
	endif
 else
	open(lun, file=filename, action="write", status="unknown")
 endif

EndSubroutine io_open


!---------------------------------------------------
!Close a unit
!---------------------------------------------------
subroutine io_close(lun)
 integer, intent(in) :: lun

 close(lun)
 if ((lun .ge. min_lun) .and. (lun .le. max_lun)) then
	lun_is_free(lun) = .true.
  else
	write(*,*) "WARNING program closed a lun that was never created in lun management."
  endif

end subroutine io_close


!---------------------------------------------------
!Close all open units
!---------------------------------------------------
subroutine io_close_all()
 integer :: i

 do i = min_lun,max_lun
	if (lun_is_free(i) .eqv. .true.) then
		close(i)
	endif
 enddo

end subroutine io_close_all


!---------------------------------------------------
!Returns the stdout unit (default = 6)
!---------------------------------------------------
integer function io_get_stdout()
	io_get_stdout = stdout
end function


!---------------------------------------------------
!Reserve a unit (useful for compatibility with legacy code)
!eg. call io_reserve(15) before open(15) appears
!---------------------------------------------------
subroutine io_reserve(lun)
      integer, intent(in) :: lun
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) then
	 write(*,*) "Error in logical unit scheme. Cannot reserve unit ", lun, " since it is in use."
      else
      	 if (lun .ge. min_lun .and. lun .le. max_lun) lun_is_free(lun) = .false.
      endif
end subroutine io_reserve

end module 
