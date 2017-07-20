!-------------------------------------------------------------------------
! File I/O routines
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
module InputOutput
 use main_vars
 use lun_management
 implicit None

 contains

!------------------------------------------------------------
!----------- Read input file parameters --------------------
!------------------------------------------------------------
subroutine read_input_file
 implicit none
 integer :: luninp, Nlines, EoF


 !call io_assign(luninp)
 !open(luninp, file="PhononSED.inp", status='old', action='read')
 luninp = 5 !read from terminal standard in
 read(luninp, *) model
 read(luninp, *) fileheader
 read(luninp, *) Nunitcells
 read(luninp, *) Neig
 read(luninp, *) Nk
 read(luninp, *) fvel
 read(luninp, *) feig
 read(luninp, *) fcoord
 read(luninp, *) Ntimesteps
 read(luninp, *) timestep
 read(luninp, *) READALL
 read(luninp, *) NPointsOut
 read(luninp, *) MaxFreqOut
 !call io_close(luninp)

 if (model == 'RDX') then
    write(*,*) "Model is RDX"
    AtomsPerMolecule = 21
    MoleculesPerUnitCell = 8
    AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !168
    Natoms = Nunitcells*AtomsPerUnitCell

    allocate(MassPrefac(AtomsPerUnitCell))
    allocate(freqs(Nk, Neig))
    allocate(eig_vecs(Nk, Neig, Natoms, 3))

    !build array of masses
    do i = 1, MoleculesPerUnitCell
        idx = (i-1)*AtomsPerMolecule + 1
        MassPrefac(idx+0:idx+2)   = MC ! Carbon
        MassPrefac(idx+3:idx+8)   = MN ! Nitrogen
        MassPrefac(idx+9:idx+14)  = MO ! Oxygen
        MassPrefac(idx+15:idx+20) = MH ! Hydrogen
    enddo
    MassPrefac = sqrt(MassPrefac/real(Nunitcells))

    lattice_vector = (/ 13.18200, 11.5740, 10.709 /)
 endif

 !------------- read length of velocities file ----------------------
 if (READALL) then
    write(*,*) "Reading size of file ..."

    Nlines = 0

    call io_assign(lun)
    open(lun, file=fvel, status='old', action='read')

    do
        read(lun, *, iostat=EoF)
        if (.not.(EoF .eq. 0)) then
            exit
        else
            Nlines = Nlines + 1
        endif
    enddo
    Ntimesteps = floor(real(Nlines)/real(Natoms+9)) - 1
    write(*, *) "There are ", Ntimesteps, " timesteps in the file"

    call io_close(lun)

endif

 allocate(qdot(Ntimesteps, Natoms, 3))
 allocate(r(Natoms, 3))
 allocate(SED(Ntimesteps))
 allocate(oneSED(Ntimesteps))

 !-------------- set up r-vectors (to be filled in) ------

end subroutine read_input_file


!------------------------------------------------------------
!---------------- Read k-point file ------------------------
!------------------------------------------------------------
subroutine read_eigvector_file()
 integer :: Natoms_file, Neig_file

 call io_assign(luneig)
 open(luneig, file=feig, status='old', action='read')

 read(luneig, *) Natoms_file
 if (.not.(AtomsPerUnitCell .eq. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in eigenvector file ( ", Natoms_file," ) does", &
                 " not match expected number of atoms per unit cell (", &
                 AtomsPerUnitCell, ")"
    stop
 endif

 do ia = 1, AtomsPerUnitCell
     read(luneig, *)
 enddo

 read(luneig, *) !int
 read(luneig, *) Neig_file  !N kpoints
 if (Neig .gt. Neig_file) then
     write(*,*) "ERROR: N_k specified is larger than the number of  &
                 eigenvectors in the input file."
    stop
 endif
 do ik = 1, Nk
     read(luneig, *) !K point at ... in BZ

     do ie = 1, Neig
        read(luneig, *) !Mode    x
        read(luneig, *) freqs(ik, ie)
        do ia = 1, AtomsPerUnitCell
            read(luneig, *) eig_vecs(ik, ie, ia, :)
        enddo
     enddo
     do ie = Neig+1, Neig_file
        read(luneig, *) !Mode    x
        read(luneig, *) !freq
        do ia = 1, AtomsPerUnitCell
            read(luneig, *)
        enddo
     enddo


 enddo

end subroutine read_eigvector_file

!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
subroutine read_LAAMPS_files
 implicit none
 integer :: lun

     !------------- read velocities file ------
     call io_assign(lun)
     open(lun, file=fvel, status='old', action='read')

     allocate(velocities(Ntimesteps, Natoms, 3))
     do t = 1, Ntimesteps
         velocities(t, :, :) = one_frame(lun)
     enddo

     call io_close(lun)

     !------------- read coordinates file ------
     call io_assign(lun)
     open(lun, file=fcoord, status='old', action='read')

     allocate(coords(Natoms, 3))
     coords = one_frame(lun)

     call io_close(lun)

end subroutine read_LAAMPS_files



!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
function one_frame(lun)
 implicit none
 real(8), dimension(Natoms, 3) :: one_frame
 integer, intent(in) :: lun
 integer ::  Natoms_file


 read(lun,*) !ITEM: TIMESTEP
 read(lun,*) !timestep
 read(lun,*) !ITEM: NUMBER OF ATOMS
 read(lun,*) Natoms_file
 if (.not.(Natoms .eq. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in velocity file ( ", Natoms_file," ) does", &
                 " not match expected total number of atoms (", &
                 Natoms, ")"
    stop
 endif
 read(lun,*) !ITEM: BOX BOUNDS
 read(lun,*) box(1, :) !bounding box
 read(lun,*) box(2, :)
 read(lun,*) box(3, :)
 read(lun,*) !comment line

 !read velocities
 do ia = 1, Natoms
    read(lun,*) (one_frame(ia, ix), ix=1,3)
 enddo

end function one_frame

!---------------------------------------------------------------------
!-----------------  Print SED ---------------------------------------
!---------------------------------------------------------------------
subroutine print_SED
 implicit none


 do ik = 1, Nk

     call io_open(lunout, filename=trim(fileheader)//"_"//trim(str(ik))//"_SED.dat")

     do i = 1, NPointsOut
         write(lunout, '(f12.2)', advance='no') freqs_smoothed(i)

         do ie = 1, Neig-1
             write(lunout, '(f16.10)', advance='no') all_SED_smoothed(ik, ie, i)
         enddo

         write(lunout, '(f16.10)', advance='yes') all_SED_smoothed(ik, Neig, i)
     enddo
     call io_close(lunout)
enddo !ie = 1, Neig


!call io_open(lunout, filename=trim(fileheader)//"_frequencies.dat")

 !do i = 1, Neig
!    write(lunout, '(f16.10)') freqs(i)
 !enddo
 !call io_close(lunout)

end subroutine print_SED

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end module InputOutput
