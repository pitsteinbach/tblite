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

!> @file tblite/scf/mixer/gambits_broyden.f90
!> Provides an electronic mixer implementation

!> Implementing Broyden mixing via GAMBITS
module tblite_scf_gambits_broyden
   use mctc_env, only : error_type, dp, sp
   use gambits, only : gambits_context_type
   use gambits_api_broyden
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_wavefunction, only : wavefunction_type
   use iso_c_binding, only : c_ptr, c_size_t
   implicit none
   private

   public :: gambits_broyden_type, new_gambits_broyden

   !> Electronic mixer using modified Broyden scheme via GAMBITS
   type, extends(mixer_type) :: gambits_broyden_type
      type(c_ptr) :: ptr
      type(gambits_context_type) :: ctx

   contains
      !> Set new object to mix
      procedure :: set_1d => set_broyden_dp, set_broyden_sp
      !> Get mixed object
      procedure :: get_1d => get_broyden_dp, get_broyden_sp
      !> Set difference between two consecutive objects to mix
      procedure :: diff_1d => diff_broyden_dp, diff_broyden_sp
      !> Perform mixing
      procedure :: next => next_broyden_dp, next_broyden_sp
      !> Get error
      procedure :: get_error => get_error_dp, get_error_sp
      !> Destroy mixer pointer
      procedure :: cleanup
      !> Update context
      procedure :: update_ctx
   end type gambits_broyden_type

contains

!> Create a new instance of the GAMBITS Broyden mixer
subroutine new_gambits_broyden(self, ndim, memory, alpha, nao, prec, error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(out) :: self
   !> Number of dimensions
   integer, intent(in) :: ndim
   !> Memory size for the mixer
   integer, intent(in) :: memory
   !> Damping parameter
   real(dp), intent(in) :: alpha
   !> Number of atomic orbitals
   integer, intent(in) :: nao
   !> Precision (0: single, 1: double)
   integer, intent(in) :: prec
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call self%ctx%setup(int(1, kind=c_size_t))
   self%ptr = c_new_broyden(self%ctx%ptr, ndim, memory, alpha, nao, prec)
   call self%update_ctx(self%ctx, error)
end subroutine new_gambits_broyden

!> Set the vector to mix
subroutine set_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call set_mixer_data_dp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine set_broyden_dp

subroutine set_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call set_mixer_data_sp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine set_broyden_sp

!> Get the differences of the mixed vector
subroutine diff_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call diff_mixer_data_dp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine diff_broyden_dp

subroutine diff_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call diff_mixer_data_sp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine diff_broyden_sp

!> Get the mixed vector
subroutine get_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(out) :: qvec(:)

   call get_mixer_data_dp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine get_broyden_dp

subroutine get_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(out) :: qvec(:)

   call get_mixer_data_sp(self%ctx%ptr, self%ptr, qvec, size(qvec))
end subroutine get_broyden_sp

subroutine next_broyden_sp(self, iscf, wfn, error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(sp), allocatable,target :: density(:,:)
   real(sp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   dptr(1:size(density)) => density
   call next_broyden_data_sp(self%ctx%ptr, self%ptr, iscf, dptr)
   call self%update_ctx(self%ctx, error)
end subroutine next_broyden_sp

subroutine next_broyden_dp(self, iscf, wfn, error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), allocatable,target :: density(:,:)
   real(dp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   dptr(1:size(density)) => density
   call next_broyden_data_dp(self%ctx%ptr, self%ptr, iscf, dptr)
   call self%update_ctx(self%ctx, error)
end subroutine next_broyden_dp

!> Get the density error
pure function get_error_dp(self, iscf) result(error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(dp) :: error

   error = get_broyden_error_dp(self%ctx%ptr, self%ptr, iscf, error)
end function get_error_dp

pure function get_error_sp(self, iscf) result(error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(sp) :: error

   error = get_broyden_error_sp(self%ctx%ptr, self%ptr, iscf, error)
end function get_error_sp

subroutine cleanup(self)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self

   call destroy_mixer(self%ctx%ptr, self%ptr)
   call self%ctx%delete()
end subroutine cleanup

subroutine update_ctx(self, ctx, error)
   !> Mixer object
   class(gambits_broyden_type), intent(inout) :: self
   !> GAMBITS context
   type(gambits_context_type), intent(in) :: ctx
   !> Error handling
   type(error_type), allocatable, optional, intent(out) :: error

   call ctx%get_message(self%msg)
   if (present(error) .and. ctx%failed()) then
      allocate(error)
      error%stat = 1
      call ctx%get_error(error%message)
   end if
end subroutine update_ctx

end module tblite_scf_gambits_broyden
