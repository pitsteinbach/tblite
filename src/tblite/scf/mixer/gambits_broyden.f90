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
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_wavefunction, only : wavefunction_type
   use iso_c_binding
   implicit none
   private

   public :: gambits_broyden_type, new_gambits_broyden

   !> Electronic mixer using modified Broyden scheme via GAMBITS
   type, extends(mixer_type) :: gambits_broyden_type
      type(c_ptr) :: ptr

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
   end type gambits_broyden_type

   interface
         type(c_ptr) function c_new_broyden(ndim, memory, alpha, nao, prec) bind(C,name="SetupBroyden")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
         integer(c_int), value :: prec
      end function c_new_broyden

      subroutine set_mixer_data_dp(mixer,target,size) bind(C,name="SetDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine set_mixer_data_dp

      subroutine set_mixer_data_sp(mixer,target,size) bind(C,name="SetDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine set_mixer_data_sp

      subroutine get_mixer_data_dp(mixer,target,size) bind(C,name="GetDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(inout) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine get_mixer_data_dp

      subroutine get_mixer_data_sp(mixer,target,size) bind(C,name="GetDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(inout) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine get_mixer_data_sp

      subroutine diff_mixer_data_dp(mixer,target,size) bind(C,name="DiffDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine diff_mixer_data_dp

      subroutine diff_mixer_data_sp(mixer,target,size) bind(C,name="DiffDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine diff_mixer_data_sp

      subroutine next_broyden_data_dp(mixer,iter,density) bind(C,name="NextDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_double), intent(in) :: density(*)
      end subroutine next_broyden_data_dp

      subroutine next_broyden_data_sp(mixer,iter,density) bind(C,name="NextSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_float), intent(in) :: density(*)
      end subroutine next_broyden_data_sp

      pure double precision function get_broyden_error_dp(mixer,iter,dummy) bind(C,name="GetErrorDP")
         use iso_c_binding
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_double), value :: dummy
      end function get_broyden_error_dp

      pure real function get_broyden_error_sp(mixer,iter,dummy) bind(C,name="GetErrorSP")
         use iso_c_binding
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_float), value :: dummy
      end function get_broyden_error_sp

      subroutine destroy_mixer(mixer) bind(C,name="Destroy")
         use iso_c_binding
         type(c_ptr), value :: mixer
      end subroutine destroy_mixer

   end interface

contains

!> Create a new instance of the GAMBITS Broyden mixer
subroutine new_gambits_broyden(self, ndim, memory, alpha, nao, prec)
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

   self%ptr = c_new_broyden(ndim, memory, alpha, nao, prec)

end subroutine new_gambits_broyden

!> Set the vector to mix
subroutine set_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call set_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine set_broyden_dp

subroutine set_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call set_mixer_data_sp(self%ptr, qvec, size(qvec))
end subroutine set_broyden_sp

!> Get the differences of the mixed vector
subroutine diff_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call diff_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine diff_broyden_dp

subroutine diff_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call diff_mixer_data_sp(self%ptr, qvec, size(qvec))
end subroutine diff_broyden_sp

!> Get the mixed vector
subroutine get_broyden_dp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(out) :: qvec(:)

   call get_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine get_broyden_dp

subroutine get_broyden_sp(self, qvec)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(out) :: qvec(:)

   call get_mixer_data_sp(self%ptr, qvec, size(qvec))
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
   call next_broyden_data_sp(self%ptr, iscf, dptr)
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
   call next_broyden_data_dp(self%ptr, iscf, dptr)
end subroutine next_broyden_dp

!> Get the density error
pure function get_error_dp(self, iscf) result(error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(dp) :: error

   error = 0.0_dp
   error = get_broyden_error_dp(self%ptr, iscf, error)
end function get_error_dp

pure function get_error_sp(self, iscf) result(error)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(sp) :: error

   error = 0.0_sp
   error = get_broyden_error_sp(self%ptr, iscf, error)
end function get_error_sp

subroutine cleanup(self)
   !> Instance of the GAMBITS Broyden mixer
   class(gambits_broyden_type), intent(inout) :: self

   call destroy_mixer(self%ptr)
end subroutine cleanup

end module tblite_scf_gambits_broyden
