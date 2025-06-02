! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @dir tblite/scf/mixer
!> Routines for implementing electronic mixing

!> @file tblite/scf/mixer.f90
!> Provides an electronic mixer implementation

!> Routines to set up electronic mixers
module tblite_scf_mixer
   use mctc_env, only : wp, sp, dp
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, &
   & not_used, orbital_resolved
   use tblite_scf_gambits_broyden, only : gambits_broyden_type
   use tblite_scf_gambits_diis, only : gambits_diis_type
   use tblite_scf_mixer_input, only : mixer_input
   use tblite_scf_mixer_broyden, only : broyden_mixer, new_broyden
   use tblite_scf_mixer_type, only : mixer_type, mixers_type
   use iso_c_binding
   implicit none
   private

   public :: new_mixer, get_mixer_dimension

   interface
      type(c_ptr) function c_new_broyden(ndim, memory, alpha, nao, prec) bind(C,name="SetupBroyden")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
         integer(c_int), value :: prec
      end function c_new_broyden

      type(c_ptr) function c_new_diis(ndim, memory, alpha, overlap, nao, runmode, io_prec, prec) bind(C,name="SetupDIIS")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         real(c_double), intent(in) :: overlap(*)
         integer(c_int), value :: nao
         integer(c_int), value :: runmode
         integer(c_int), value :: io_prec
         integer(c_int), value :: prec
      end function c_new_diis
   end interface

contains

!> Create a new instance of the mixer
   subroutine new_mixer(self, input, ndim, nao, nspin, overlap, info)
      !> Instance of the mixer on exit
      class(mixers_type), intent(out) :: self
      !> Mixer input parameters
      type(mixer_input), intent(in) :: input
      !> Dimensions of the Broyden mixers
      integer, intent(in) :: ndim
      !> Number of atomic orbitals
      integer, intent(in) :: nao
      !> Number of spin channels
      integer, intent(in) :: nspin
      !> Overlap matrix for DIIS
      real(wp), intent(in) :: overlap(:,:)
      !> Information on wavefunction data used to construct Hamiltonian
      type(scf_info), intent(inout) :: info

      integer :: prec, i

      write(*,*) "In new_mixer, with type"
      if (wp == sp) then
         prec = 0
      else if (wp == dp) then
         prec = 1
      end if

      select case(input%type)
       case(0)
         write(*,*) "In case 0"
         block
            type(broyden_mixer), allocatable :: mixer
            allocate(mixer)
            call new_broyden(mixer, nspin*ndim, input%broyden)
            mixer%info = info
            allocate(self%mixer(1),source=mixer)
         end block
         allocate(self%type(1))
         self%type(1) = 0
         write(*,*) "End of case 0"
       case(1)
         block
            type(gambits_broyden_type), allocatable :: mixer
            allocate(mixer)
            mixer%ptr = c_new_broyden(nspin*ndim, input%gambits_broyden%memory, input%gambits_broyden%damp, nao, prec)
            mixer%info = info
            allocate(self%mixer(1),source=mixer)
         end block
         allocate(self%type(1))
         self%type(1) = 1

       case(2)
         block
            type(gambits_diis_type), allocatable :: mixer(:)
            allocate(mixer(nspin))
            do i=1,nspin
               mixer(i)%ptr = c_new_diis(nao**2, input%gambits_diis%memory, input%gambits_diis%damp, overlap, nao, &
               & input%gambits_diis%runmode, prec, input%gambits_diis%prec)
               call mixer(i)%diis_info(info)
               mixer(i)%nspin = nspin
               mixer(i)%channel = i
            end do
            allocate(self%mixer,source=mixer)
         end block
         allocate(self%type(nspin))
         self%type(1) = 2
      end select

   end subroutine new_mixer

   function get_mixer_dimension(mol, bas, info) result(ndim)
      use mctc_io, only : structure_type
      use tblite_basis, only : basis_type
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      type(scf_info), intent(in) :: info
      integer :: ndim

      ndim = 0

      select case(info%charge)
       case(atom_resolved)
         ndim = ndim + mol%nat
       case(shell_resolved)
         ndim = ndim + bas%nsh
      end select

      select case(info%dipole)
       case(atom_resolved)
         ndim = ndim + 3*mol%nat
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         ndim = ndim + 6*mol%nat
      end select
   end function get_mixer_dimension

end module tblite_scf_mixer
