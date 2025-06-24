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

!> @file tblite/scf/mixer/gambits_diis.f90
!> Provides an electronic mixer implementation

!> Implementing DIIS (direct inversion in the iterative subspace) mixing
!> via GAMBITS
module tblite_scf_gambits_diis
   use gambits_api_context, only : gambits_context_type
   use mctc_env, only : error_type, wp, dp, sp
   use tblite_scf_info, only : scf_info, not_used, orbital_resolved
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_wavefunction, only : wavefunction_type
   use iso_c_binding
   implicit none
   private

   public :: gambits_diis_type, new_gambits_diis

   !> Electronic mixer using DIIS scheme via GAMBITS
   type, extends(mixer_type) :: gambits_diis_type

      type(c_ptr) :: ptr
      integer :: nspin
      integer :: channel

   contains
      !> Set information on wavefunction data to be used for mixing
      procedure :: diis_info
      !> Set new object to mix
      procedure :: set_1d => set_diis_dp, set_diis_sp
      !> Get mixed object
      procedure :: get_1d => get_diis_dp, get_diis_sp
      !> Set difference between two consecutive objects to mix
      procedure :: diff_1d => diff_diis_dp, diff_diis_sp
      !> Perform mixing
      procedure :: next => next_diis_dp, next_diis_sp
      !> Get error
      procedure :: get_error => get_error_dp, get_error_sp
      !> Get timing information
      procedure :: get_timings
      !> Destroy mixer pointer
      procedure :: cleanup
      !> Update context
      procedure :: update_ctx
   end type gambits_diis_type

   interface
      type(c_ptr) function c_new_diis(ctx, ndim, memory, alpha, overlap, nao, runmode, io_prec, prec) bind(C,name="SetupDIIS")
         use iso_c_binding
         type(c_ptr), value :: ctx
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         real(c_double), intent(in) :: overlap(*)
         integer(c_int), value :: nao
         integer(c_int), value :: runmode
         integer(c_int), value :: io_prec
         integer(c_int), value :: prec
      end function c_new_diis

      subroutine set_mixer_data_dp(ctx, mixer, target, size) bind(C,name="SetDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine set_mixer_data_dp

      subroutine set_mixer_data_sp(ctx, mixer, target, size) bind(C,name="SetDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine set_mixer_data_sp

      subroutine get_mixer_data_dp(ctx, mixer, target, size) bind(C,name="GetDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(inout) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine get_mixer_data_dp

      subroutine get_mixer_data_sp(ctx, mixer, target, size) bind(C,name="GetDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(inout) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine get_mixer_data_sp

      subroutine diff_mixer_data_dp(ctx, mixer, target, size) bind(C,name="DiffDataDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine diff_mixer_data_dp

      subroutine diff_mixer_data_sp(ctx, mixer, target, size) bind(C,name="DiffDataSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine diff_mixer_data_sp

      subroutine next_diis_data_dp(ctx, mixer, iter, density) bind(C,name="NextDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_double), intent(in) :: density(*)
      end subroutine next_diis_data_dp

      subroutine next_diis_data_sp(ctx, mixer, iter, density) bind(C,name="NextSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: ctx
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_float), intent(in) :: density(*)
      end subroutine next_diis_data_sp

      double precision function get_diis_error_dp(ctx, mixer, iter, dummy) bind(C,name="GetErrorDP")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_double), value :: dummy
      end function get_diis_error_dp

      real function get_diis_error_sp(ctx, mixer, iter, dummy) bind(C,name="GetErrorSP")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_float), value :: dummy
      end function get_diis_error_sp

      subroutine destroy_mixer(ctx, mixer) bind(C,name="Destroy")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: mixer
      end subroutine destroy_mixer

      subroutine get_timing(ctx, mixer) bind(C,name="GetTiming")
         use iso_c_binding
         type(c_ptr), value :: ctx
         type(c_ptr), value :: mixer
      end subroutine get_timing

   end interface

contains

subroutine new_gambits_diis(self, ctx, ndim, memory, alpha, overlap, nao, runmode, io_prec, prec, error)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(out) :: self
   !> GAMBITS context
   type(gambits_context_type), intent(in) :: ctx
   !> Number of dimensions for the mixer
   integer, intent(in) :: ndim
   !> Memory size for the mixer
   integer, intent(in) :: memory
   !> Damping parameter for the mixer
   real(wp), intent(in) :: alpha
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:,:)
   !> Number of atomic orbitals
   integer, intent(in) :: nao
   !> Runmode (size-dependent: 0, cpu: 1, gpu: 2)
   integer, intent(in) :: runmode
   !> IO precision (FP32: 0, FP64: 1)
   integer, intent(in) :: io_prec
   !> Working precision (FP32: 0, FP64: 1)
   integer, intent(in) :: prec
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   self%ptr = c_new_diis(ctx%ptr, ndim, memory, alpha, overlap, nao, runmode, io_prec, prec)
   call self%update_ctx(ctx, error)
end subroutine new_gambits_diis

subroutine diis_info(self, info)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Information on wavefunction data used to construct Hamiltonian
   type(scf_info), intent(in) :: info

   self%info = info
   self%info%charge = not_used
   self%info%dipole = not_used
   self%info%quadrupole = not_used
   self%info%density = not_used
   self%info%fock = orbital_resolved
end subroutine diis_info

!> Set the vector to mix
subroutine set_diis_dp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call set_mixer_data_dp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine set_diis_dp

subroutine set_diis_sp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call set_mixer_data_sp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine set_diis_sp

!> Get the differences of the mixed vector
subroutine diff_diis_dp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call diff_mixer_data_dp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine diff_diis_dp

subroutine diff_diis_sp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call diff_mixer_data_sp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine diff_diis_sp

!> Get the mixed vector
subroutine get_diis_dp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(out) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call get_mixer_data_dp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine get_diis_dp

subroutine get_diis_sp(self, qvec, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(out) :: qvec(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   call get_mixer_data_sp(ctx%ptr, self%ptr, qvec, size(qvec))
   call self%update_ctx(ctx, error)
end subroutine get_diis_sp

subroutine next_diis_dp(self, iscf, wfn, error, ctx)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   real(dp), target, allocatable :: density(:, :)
   real(dp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   density = wfn%density(:,:,self%channel)
   dptr(1:size(density)) => density
   call next_diis_data_dp(ctx%ptr, self%ptr, iscf, dptr)
   call self%update_ctx(ctx, error)
end subroutine next_diis_dp

subroutine next_diis_sp(self, iscf, wfn, error, ctx)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   real(sp), target, allocatable :: density(:, :)
   real(sp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   density = wfn%density(:,:,self%channel)
   dptr(1:size(density)) => density
   call next_diis_data_sp(ctx%ptr, self%ptr, iscf, dptr)
   call self%update_ctx(ctx, error)
end subroutine next_diis_sp

!> Get the density error
function get_error_dp(self, iscf, error, ctx) result(err)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Current iteration
   integer, intent(in) :: iscf
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   real(dp) :: err

   err = 0.0_dp
   err = get_diis_error_dp(ctx%ptr, self%ptr, iscf, err)
   err = err * (self%nspin*1.0)**2
   call self%update_ctx(ctx, error)
end function get_error_dp

!> Get the density error
function get_error_sp(self, iscf, error, ctx) result(err)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Current iteration
   integer, intent(in) :: iscf
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(in), optional :: ctx

   real(sp) :: err

   err = 0.0_sp
   err = get_diis_error_sp(ctx%ptr, self%ptr, iscf, err)
   err = err * (self%nspin*1.0)**2
   call self%update_ctx(ctx, error)
end function get_error_sp

subroutine cleanup(self, error, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> GAMBITS context
   type(gambits_context_type), intent(inout), optional :: ctx

   call destroy_mixer(ctx%ptr, self%ptr)
   call self%update_ctx(ctx, error)
end subroutine cleanup

subroutine get_timings(self, ctx)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> GAMBITS context
   type(gambits_context_type), intent(inout), optional :: ctx

   call get_timing(ctx%ptr, self%ptr)
   call ctx%get_message(self%msg)
end subroutine get_timings

subroutine update_ctx(self, ctx, error)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self
   !> GAMBITS context
   type(gambits_context_type), intent(in) :: ctx
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call ctx%get_message(self%msg)
   if (ctx%failed()) then
      allocate(error)
      error%stat = 1
      call ctx%get_error(error%message)
   end if
end subroutine update_ctx

end module tblite_scf_gambits_diis
