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

!> @file tblite/scf/mixer/type.f90
!> Provides an electronic mixer implementation

!> Base class for electronic mixing
module tblite_scf_mixer_type
   use mctc_env, only : error_type, wp
   use tblite_basis_type, only : basis_type
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, orbital_resolved
   use tblite_scf_utils, only : get_qat_from_qsh
   use tblite_wavefunction, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for electronic mixing
   type, public, abstract :: mixer_type
      type(scf_info) :: info
   contains
      !> Apply mixing to the density
      procedure(next), deferred :: next
      !> Set new density
      generic :: set => set_1d, set_2d, set_3d
      !> Set new density from 1D array
      procedure(set_1d), deferred :: set_1d
      !> Set new density from 2D array
      procedure :: set_2d
      !> Set new density from 3D array
      procedure :: set_3d
      !> Set difference between new and old density
      generic :: diff => diff_1d, diff_2d, diff_3d
      !> Set difference between new and old density from 1D array
      procedure(diff_1d), deferred :: diff_1d
      !> Set difference between new and old density from 2D array
      procedure :: diff_2d
      !> Set difference between new and old density from 3D array
      procedure :: diff_3d
      !> Get density
      generic :: get => get_1d, get_2d, get_3d
      !> Get density as 1D array
      procedure(get_1d), deferred :: get_1d
      !> Get density as 2D array
      procedure :: get_2d
      !> Get density as 3D array
      procedure :: get_3d
      !> Get error metric from mixing
      procedure(get_error), deferred :: get_error
      !> Destroy mixer
      procedure :: cleanup
   end type mixer_type

   abstract interface
      !> Apply mixing to the density
      subroutine next(self, iscf, wfn, error)
         import :: mixer_type, wavefunction_type, error_type
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Iteration counter
         integer, intent(in) :: iscf
         !> Tight-binding wavefunction data
         type(wavefunction_type), intent(inout) :: wfn
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine next

      !> Set new density from 1D array
      subroutine set_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(in) :: qvec(:)
      end subroutine set_1d

      !> Set difference between new and old density from 1D array
      subroutine diff_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(in) :: qvec(:)
      end subroutine diff_1d

      !> Get density as 1D array
      subroutine get_1d(self, qvec)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(inout) :: self
         !> Density vector
         real(wp), intent(out) :: qvec(:)
      end subroutine get_1d

      !> Get error metric from mixing
      pure function get_error(self,iscf) result(error)
         import :: mixer_type, wp
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Iteration counter
         integer, intent(in) :: iscf
         !> Error metric
         real(wp) :: error
      end function get_error
   end interface


   type, public :: mixers_type
      !> List of mixers
      class(mixer_type), allocatable :: mixer(:)
      !> List of mixer types
      integer, allocatable :: type(:)
   contains
      !> Apply mixing to the density
      procedure :: next_mixer
      !> Get error metric from mixing
      procedure :: get_error_mixer
      !> Destroy mixer
      procedure :: cleanup_mixer
      !> Set new density
      procedure :: set_mixer
      !> Set difference between new and old density
      procedure :: diff_mixer
      !> Get density
      procedure :: get_mixer
   end type mixers_type

contains

!> Set new density from 2D array
subroutine set_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_2d

!> Set new density from 3D array
subroutine set_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%set(qptr)
end subroutine set_3d

!> Set difference between new and old density from 2D array
subroutine diff_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_2d

!> Set difference between new and old density from 3D array
subroutine diff_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(in), target :: qvec(:, :, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%diff(qptr)
end subroutine diff_3d

!> Get density as 2D array
subroutine get_2d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(out), target :: qvec(:, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_2d

!> Get density as 3D array
subroutine get_3d(self, qvec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
   !> Density vector
   real(wp), contiguous, intent(out), target :: qvec(:, :, :)

   real(wp), pointer :: qptr(:)

   qptr(1:size(qvec)) => qvec
   call self%get(qptr)
end subroutine get_3d

!> Cleanup
subroutine cleanup(self)
   !> Instance of the electronic mixer
   class(mixer_type), intent(inout) :: self
end subroutine cleanup

subroutine next_mixer(self, iscf, wfn, error)
   !> Instance of the electronic mixer
   class(mixers_type), intent(inout) :: self
   !> Iteration counter
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: channel

   do channel = 1, size(self%mixer)
      call self%mixer(channel)%next(iscf, wfn, error)
   end do
end subroutine next_mixer


pure function get_error_mixer(self, iscf) result(error)
   !> Instance of the electronic mixer
   class(mixers_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   integer :: channel
   real(wp) :: error
   real(wp) :: perr(size(self%type))

   do channel = 1, size(self%mixer)
      perr(channel) = self%mixer(channel)%get_error(iscf)
   end do

   error = maxval(abs(perr))
end function get_error_mixer

subroutine cleanup_mixer(self)
   !> Instance of the electronic mixer
   class(mixers_type), intent(inout) :: self

   integer :: i

   do i = 1, size(self%mixer)
      call self%mixer(i)%cleanup()
   end do
end subroutine cleanup_mixer

subroutine set_mixer(self, wfn)
   !> Instance of the electronic mixer
   class(mixers_type), intent(inout) :: self
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn

   integer :: channel

   do channel = 1, size(self%mixer)
      select case(self%mixer(channel)%info%charge)
         case(atom_resolved)
         call self%mixer(channel)%set(wfn%qat)
         case(shell_resolved)
         call self%mixer(channel)%set(wfn%qsh)
      end select

      select case(self%mixer(channel)%info%dipole)
         case(atom_resolved)
         call self%mixer(channel)%set(wfn%dpat)
      end select

      select case(self%mixer(channel)%info%quadrupole)
         case(atom_resolved)
         call self%mixer(channel)%set(wfn%qpat)
      end select

      select case(self%mixer(channel)%info%density)
         case(orbital_resolved)
         call self%mixer(channel)%set(wfn%density(:,:,channel))
      end select

      select case(self%mixer(channel)%info%fock)
         case(orbital_resolved)
         call self%mixer(channel)%set(wfn%coeff(:,:,channel))
      end select
   end do
end subroutine set_mixer

subroutine diff_mixer(self, wfn)
   !> Instance of the electronic mixer
   class(mixers_type), intent(inout) :: self
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn

   integer :: channel

   do channel = 1, size(self%mixer)
      select case(self%mixer(channel)%info%charge)
         case(atom_resolved)
         call self%mixer(channel)%diff(wfn%qat)
         case(shell_resolved)
         call self%mixer(channel)%diff(wfn%qsh)
      end select

      select case(self%mixer(channel)%info%dipole)
         case(atom_resolved)
         call self%mixer(channel)%diff(wfn%dpat)
      end select

      select case(self%mixer(channel)%info%quadrupole)
         case(atom_resolved)
         call self%mixer(channel)%diff(wfn%qpat)
      end select

      select case(self%mixer(channel)%info%density)
         case(orbital_resolved)
         call self%mixer(channel)%diff(wfn%density(:,:,channel))
      end select

      select case(self%mixer(channel)%info%fock)
         case(orbital_resolved)
         call self%mixer(channel)%diff(wfn%coeff(:,:,channel))
      end select
   end do
end subroutine diff_mixer

subroutine get_mixer(self, bas, wfn)
   !> Instance of the electronic mixer
   class(mixers_type), intent(inout) :: self
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn

   integer :: channel

   do channel = 1, size(self%mixer)
      select case(self%mixer(channel)%info%charge)
         case(atom_resolved)
         call self%mixer(channel)%get(wfn%qat)
         case(shell_resolved)
         call self%mixer(channel)%get(wfn%qsh)
         call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      end select

      select case(self%mixer(channel)%info%dipole)
         case(atom_resolved)
         call self%mixer(channel)%get(wfn%dpat)
      end select

      select case(self%mixer(channel)%info%quadrupole)
         case(atom_resolved)
         call self%mixer(channel)%get(wfn%qpat)
      end select

      select case(self%mixer(channel)%info%density)
         case(orbital_resolved)
         call self%mixer(channel)%get(wfn%density(:,:,channel))
      end select

      select case(self%mixer(channel)%info%fock)
         case(orbital_resolved)
         call self%mixer(channel)%get(wfn%coeff(:,:,channel))
      end select
   end do
end subroutine get_mixer

end module tblite_scf_mixer_type
