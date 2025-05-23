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
!> Proxy module for electronic mixing routines

!> Provides an electronic mixer implementation
module tblite_scf_mixer
   use mctc_env, only : wp
   use tblite_scf_mixer_broyden, only : broyden_mixer, broyden_input, new_broyden
   use tblite_basis, only : basis_type
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_scf_info, only : scf_info
   use tblite_wavefunction, only : wavefunction_type
   implicit none
   private

   public :: mixer_type, new_mixer, get_mixer_dimension
   public :: get_mixer, set_mixer, diff_mixer


   !> Input for selecting electronic mixer
   type, public :: mixer_input
      !> Input for Broyden mixer
      type(broyden_input), allocatable :: broyden
   end type mixer_input

contains

!> Create a new instance of the mixer
subroutine new_mixer(self, memory, ndim, damp)
   !> Instance of the mixer on exit
   class(mixer_type), allocatable, intent(out) :: self
   integer, intent(in) :: memory
   integer, intent(in) :: ndim
   real(wp), intent(in) :: damp

   block
      type(broyden_mixer), allocatable :: mixer
      allocate(mixer)
      call new_broyden(mixer, ndim, broyden_input(memory, damp))
      call move_alloc(mixer, self)
   end block
end subroutine new_mixer

function get_mixer_dimension(mol, bas, info) result(ndim)
   use mctc_io, only : structure_type
   use tblite_scf_info, only : atom_resolved, shell_resolved
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

subroutine set_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%set(wfn%qat)
   case(shell_resolved)
      call mixer%set(wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%set(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%set(wfn%qpat)
   end select
end subroutine set_mixer

subroutine diff_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   class(mixer_type), intent(inout) :: mixer
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%diff(wfn%qat)
   case(shell_resolved)
      call mixer%diff(wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%diff(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%diff(wfn%qpat)
   end select
end subroutine diff_mixer

subroutine get_mixer(mixer, bas, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   use tblite_scf_iterator, only : get_qat_from_qsh
   class(mixer_type), intent(inout) :: mixer
   type(basis_type), intent(in) :: bas
   type(wavefunction_type), intent(inout) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%get(wfn%qat)
   case(shell_resolved)
      call mixer%get(wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%get(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%get(wfn%qpat)
   end select
end subroutine get_mixer


end module tblite_scf_mixer
