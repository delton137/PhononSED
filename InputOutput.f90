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
 integer :: luninp

 call io_assign(luninp)
 open(luninp, file="PhononSED.inp", status='old', action='read')
 read(luninp, *) model
 read(luninp, *) Nunitcells
 read(luninp, *) fvel
 read(luninp, *) feig
 read(luninp, *) Ntimesteps
 read(luninp, *) READALL
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
    Nk = 504

    allocate(MassPrefac(Natoms))
    allocate(freqs(Nk))
    allocate(eig_vecs(Nk, Natoms, 3))

    !build array of masses
    do i = 1, MoleculesPerUnitCell
        idx = (i-1)*AtomsPerMolecule
        MassPrefac(idx+0:idx+2)   = MC ! Carbon
        MassPrefac(idx+3:idx+8)   = MN ! Nitrogen
        MassPrefac(idx+9:idx+14)  = MO ! Oxygen
        MassPrefac(idx+15:idx+20) = MH ! Hydrogen
    enddo
    MassPrefac = sqrt(MassPrefac/Nunitcells)
 endif


end subroutine read_input_file

!------------------------------------------------------------
!---------------- Read k-point file ------------------------
!------------------------------------------------------------
subroutine read_eigvector_file()
 integer :: Natoms_file, EOF=0, Nk_file

 read(luneig, *) Natoms_file
 if (.not.(Natoms .eq. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in eigenvector file does not match expected \\
                 number of atoms per unit cell"
    stop
 endif
 do ia = 1, Natoms
     read(luneig, *)
 enddo

 ik = 1
 read(luneig, *) !int
 read(luneig, *) Nk_file  !N kpoints
 if (.not.(Nk .eq. Nk_file)) then
     write(*,*) "ERROR: N_k in eigenvector file does not match expected \\
                 number of eigenvectors"
    stop
 endif
 read(luneig, *) !K point at 0.0... in BZ

 do while (EOF .eq. 0)
    read(luneig, *) !Mode    x
    read(luneig, *) freqs(ik)
    do ia = 1, Natoms
        read(luneig, '(3xf10.6)', iostat = EOF)  eig_vecs(ik, ia, : )
    enddo
    ik = ik + 1
 enddo

end subroutine read_eigvector_file

!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
subroutine read_velocities_file(lun)
 implicit none
 integer, intent(in) :: lun

 if (READALL .eqv. .false.) then
     allocate(velocities(Ntimesteps, Natoms, 3))
     do t = 1, Ntimesteps
         call read_one_frame(lunvel, t)
     enddo
 endif

 if (READALL .eqv. .true.) then
     !! to be filled in
 endif

end subroutine read_velocities_file


!-----------------------------------------------------------------------
!----------------- Read one frame of .xyz velocities file -------------
!-----------------------------------------------------------------------
subroutine read_one_frame(lun, t)
 implicit none
 integer, intent(in) :: lun, t
 integer ::  Natoms_file


 read(lun,*) !comment line
 read(lun,*) Natoms_file
 if (.not.(Natoms .eq. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in eigenvector file does not match expected \\
                 number of atoms per unit cell"
    stop
 endif

 read(lun,*) box(1, :) !bounding box
 read(lun,*) box(2, :)
 read(lun,*) box(3, :)
 read(lun,*) !comment line

    !read velocities
 do ia = 1, Natoms
    read(lun,*) (velocities(t, ia, ix), ix=1,3)
 enddo

end subroutine read_one_frame



end module InputOutput
