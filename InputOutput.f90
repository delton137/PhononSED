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


 call io_assign(luninp)
 open(luninp, file="PhononSED.inp", status='old', action='read')
 read(luninp, *) model
 read(luninp, *) fileheader
 read(luninp, *) Nunitcells
 read(luninp, *) Nk
 read(luninp, *) fvel
 read(luninp, *) feig
 read(luninp, *) Ntimesteps
 read(luninp, *) timestep
 read(luninp, *) READALL
 read(luninp, *) NPointsOut
 read(luninp, *) MaxFreqOut
 call io_close(luninp)

 call io_assign(lunvel)
 open(lunvel, file=fvel, status='old', action='read')

 call io_assign(luneig)
 open(luneig, file=feig, status='old', action='read')

 if (model == 'RDX') then
    write(*,*) "Model is RDX"
    AtomsPerMolecule = 21
    MoleculesPerUnitCell = 8
    AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !168
    Natoms = Nunitcells*AtomsPerUnitCell

    allocate(MassPrefac(AtomsPerUnitCell))
    allocate(freqs(Nk))
    allocate(eig_vecs(Nk, Natoms, 3))

    !build array of masses
    do i = 1, MoleculesPerUnitCell
        idx = (i-1)*AtomsPerMolecule + 1
        MassPrefac(idx+0:idx+2)   = MC ! Carbon
        MassPrefac(idx+3:idx+8)   = MN ! Nitrogen
        MassPrefac(idx+9:idx+14)  = MO ! Oxygen
        MassPrefac(idx+15:idx+20) = MH ! Hydrogen
    enddo
    MassPrefac = sqrt(MassPrefac/real(Nunitcells))
 endif

 !------------- read length of velocities file ----------------------
 if (READALL) then
    write(*,*) "Reading size of file ..."

    Nlines = 0
    do
        read(lunvel, *, iostat=EoF)
        if (.not.(EoF .eq. 0)) then
            exit
        else
            Nlines = Nlines + 1
        endif
    enddo
    Ntimesteps = floor(real(Nlines)/real(Natoms+9)) - 1
    write(*, *) "There are ", Ntimesteps, " timesteps in the file"

endif

 allocate(qdot(Ntimesteps, Natoms, 3))
 allocate(r(Nunitcells, 3))
 allocate(SED(Ntimesteps))
 allocate(oneSED(Ntimesteps))

 !-------------- set up r-vectors (to be filled in) ------

end subroutine read_input_file


!------------------------------------------------------------
!---------------- Read k-point file ------------------------
!------------------------------------------------------------
subroutine read_eigvector_file()
 integer :: Natoms_file, Nk_file

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
 read(luneig, *) Nk_file  !N kpoints
 !if (.not.(Nk .eq. Nk_file)) then
!     write(*,*) "ERROR: N_k in eigenvector file does not match expected &
!                 number of eigenvectors"
!    stop
! endif

 read(luneig, *) !K point at 0.0... in BZ

 do ik = 1, Nk
    read(luneig, *) !Mode    x
    read(luneig, *) freqs(ik)
    do ia = 1, AtomsPerUnitCell
        read(luneig, *) eig_vecs(ik, ia, :)
    enddo
 enddo

end subroutine read_eigvector_file

!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
subroutine read_velocities_file(lun)
 implicit none
 integer, intent(in) :: lun

     allocate(velocities(Ntimesteps, Natoms, 3))
     do t = 1, Ntimesteps
         call read_one_frame(lunvel, t)
     enddo

end subroutine read_velocities_file


!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
subroutine read_one_frame(lun, t)
 implicit none
 integer, intent(in) :: lun, t
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
    read(lun,*) (velocities(t, ia, ix), ix=1,3)
 enddo


end subroutine read_one_frame

!---------------------------------------------------------------------
!-----------------  Print SED ---------------------------------------
!---------------------------------------------------------------------
subroutine print_SED
 implicit none

 call io_open(lunout, filename=trim(fileheader)//"_SED.dat")

 do i = 1, NPointsOut
    write(lunout, '(f12.2)', advance='no') freqs_smoothed(i)

    do j = 1, Nk-1
        write(lunout, '(f16.10)', advance='no') all_SED_smoothed(j, i)
    enddo

    write(lunout, '(f16.10)', advance='yes') all_SED_smoothed(Nk, i)
enddo

 call io_close(lunout)
 call io_open(lunout, filename=trim(fileheader)//"_frequencies.dat")

 do i = 1, Nk
    write(lunout, '(f16.10)') freqs(i)
 enddo

 call io_close(lunout)

end subroutine print_SED

end module InputOutput
