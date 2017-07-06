!-------------------------------------------------------------------------------------
! Global variable storage area
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
module main_vars
 implicit none
 integer :: Natoms, Nunitcells, AtomsPerMolecule, MoleculesPerUnitCell, AtomsPerUnitCell
 integer :: Nk, lunvel, luneig, ik, ia, i, j, idx, ix, Ntimesteps, t
 real(8), dimension(3,3) :: box
 logical :: READALL
 real(8), dimension(:), allocatable :: MassPrefac, freqs
 real(8), dimension(:,:,:), allocatable :: eig_vecs, velocities
 double complex, dimension(:,:,:), allocatable :: qdot
 real(8) :: MC = 12.011000, MN = 14.007200, MO = 15.999430, MH = 1.0080000
 character(len=200) :: model, fvel, feig


end module main_vars
