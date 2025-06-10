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
      !> Destroy mixer pointer
      procedure :: cleanup
   end type gambits_diis_type

   interface
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

      subroutine set_error_vec_dp(mixer,error, xerr, yerr) bind(C,name="SetErrorDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: error(*)
         integer(c_int), intent(in), value :: xerr
         integer(c_int), intent(in), value :: yerr
      end subroutine set_error_vec_dp

      subroutine set_error_vec_sp(mixer,error, xerr, yerr) bind(C,name="SetErrorSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: error(*)
         integer(c_int), intent(in), value :: xerr
         integer(c_int), intent(in), value :: yerr
      end subroutine set_error_vec_sp

      subroutine next_diis_data_dp(mixer,iter,density) bind(C,name="NextDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_double), intent(in) :: density(*)
      end subroutine next_diis_data_dp

      subroutine next_diis_data_sp(mixer,iter,density) bind(C,name="NextSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         integer(c_int), value, intent(in) :: iter
         real(c_float), intent(in) :: density(*)
      end subroutine next_diis_data_sp

      pure double precision function get_diis_error_dp(mixer,iter,dummy) bind(C,name="GetErrorDP")
         use iso_c_binding
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_double), value :: dummy
      end function get_diis_error_dp

      pure real function get_diis_error_sp(mixer,iter,dummy) bind(C,name="GetErrorSP")
         use iso_c_binding
         type(c_ptr), value :: mixer
         integer(c_int), value :: iter
         real(c_float), value :: dummy
      end function get_diis_error_sp

      subroutine destroy_mixer(mixer) bind(C,name="Destroy")
         use iso_c_binding
         type(c_ptr), value :: mixer
      end subroutine destroy_mixer

   end interface

contains

subroutine new_gambits_diis(self, ndim, memory, alpha, overlap, nao, runmode, io_prec, prec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(out) :: self
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

   self%ptr = c_new_diis(ndim, memory, alpha, overlap, nao, runmode, io_prec, prec)

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
subroutine set_diis_dp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call set_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine set_diis_dp

subroutine set_diis_sp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call set_mixer_data_sp(self%ptr, qvec, size(qvec))
end subroutine set_diis_sp


!> Get the differences of the mixed vector
subroutine diff_diis_dp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(in) :: qvec(:)

   call diff_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine diff_diis_dp

subroutine diff_diis_sp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(in) :: qvec(:)

   call diff_mixer_data_sp(self%ptr, qvec, size(qvec))
end subroutine diff_diis_sp

!> Get the mixed vector
subroutine get_diis_dp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(dp), intent(out) :: qvec(:)

   call get_mixer_data_dp(self%ptr, qvec, size(qvec))
end subroutine get_diis_dp

subroutine get_diis_sp(self, qvec)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(inout) :: self
   !> Density vector
   real(sp), intent(out) :: qvec(:)

   call get_mixer_data_sp(self%ptr, qvec, size(qvec))
end subroutine get_diis_sp

subroutine next_diis_dp(self, iscf, wfn, error)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), target, allocatable :: density(:, :)
   real(dp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   density = wfn%density(:,:,self%channel)
   dptr(1:size(density)) => density
   call next_diis_data_dp(self%ptr, iscf, dptr)
end subroutine next_diis_dp

subroutine next_diis_sp(self, iscf, wfn, error)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self
   !> SCF Iteration
   integer, intent(in) :: iscf
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(sp), target, allocatable :: density(:, :)
   real(sp), pointer :: dptr(:)

   allocate(density(size(wfn%density,1),size(wfn%density,2)))
   density = wfn%density(:,:,self%channel)
   dptr(1:size(density)) => density
   call next_diis_data_sp(self%ptr, iscf, dptr)
end subroutine next_diis_sp

!> Get the density error
pure function get_error_dp(self, iscf) result(error)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(dp) :: error

   error = 0.0_dp
   error = get_diis_error_dp(self%ptr, iscf, error)
   error = error * (self%nspin*1.0)**2
end function get_error_dp

!> Get the density error
pure function get_error_sp(self, iscf) result(error)
   !> Instance of the GAMBITS DIIS mixer
   class(gambits_diis_type), intent(in) :: self
   !> Current iteration
   integer, intent(in) :: iscf

   real(sp) :: error

   error = 0.0_sp
   error = get_diis_error_sp(self%ptr, iscf, error)
   error = error * (self%nspin*1.0)**2
end function get_error_sp

subroutine cleanup(self)
   !> Mixer object
   class(gambits_diis_type), intent(inout) :: self

   call destroy_mixer(self%ptr)
end subroutine cleanup

end module tblite_scf_gambits_diis
