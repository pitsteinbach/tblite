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

!> @file tblite/scf/mixers/diis.f90
!> Implementing DIIS (direct inversion in the iterative subspace) mixing

module tblite_scf_mixer_diis
   use tblite_scf_mixer
   use iso_c_binding
   implicit none

   !> DIIS mixer
   type, extends(mixer_type) :: diis_type

   contains
      !> Set new object to mix
      procedure :: set => set_diis_dp, set_diis_sp
      !> Set difference between two consecutive objects to mix
      procedure :: diff => diff_diis_dp, diff_diis_sp
      !> Get mixed object
      procedure :: get => get_diis_dp, get_diis_sp
      !> Construct the error vector
      procedure :: construct_error => construct_error_dp, construct_error_sp
      !> Set the error vector manually
      procedure :: set_error => set_error_dp, set_error_sp
   end type diis_type

   interface
      type(c_ptr) function c_new_diis(ndim, memory, alpha, nao) bind(C,name="SetupDIIS")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
      end function c_new_diis

      subroutine construct_error_vec_dp(mixer,fock,density,overlap) bind(C,name="ConstructErrorDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: fock(*)
         real(c_double), intent(in) :: density(*)
         real(c_double), intent(in) :: overlap(*)
      end subroutine construct_error_vec_dp

      subroutine construct_error_vec_sp(mixer,fock,density,overlap) bind(C,name="ConstructErrorSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: fock(*)
         real(c_float), intent(in) :: density(*)
         real(c_float), intent(in) :: overlap(*)
      end subroutine construct_error_vec_sp

      subroutine set_mixer_data_dp(mixer,fock) bind(C,name="SetQDP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_double), intent(in) :: fock(*)
      end subroutine set_mixer_data_dp

      subroutine set_mixer_data_sp(mixer,fock) bind(C,name="SetQSP")
         use iso_c_binding
         type(c_ptr), value, intent(in) :: mixer
         real(c_float), intent(in) :: fock(*)
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


   end interface

contains
   !> Create a new instance of the DIIS mixer
   subroutine new_diis(self, mol, calc, info)
      use tblite_xtb_calculator, only : xtb_calculator
      use mctc_io, only : structure_type
      use tblite_scf_info, only : scf_info
      !> DIIS object
      class(diis_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
      !> Info data
      type(scf_info) :: info

      self%ndim = calc%bas%nao**2
      self%memory = calc%mixer_mem
      self%ptr = c_new_diis(self%ndim, self%memory, calc%mixer_damping, calc%bas%nao)
   end subroutine new_diis

   !> Calculate the error vector used for mixing
   subroutine construct_error_dp(self, hmat, pmat, smat)
      use mctc_env, only : dp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(dp), intent(in) :: hmat(:,:)
      !> Density matrix
      real(dp), intent(in) :: pmat(:,:)
      !> Overlap matrix
      real(dp), intent(in) :: smat(:,:)

      call construct_error_vec_dp(self%ptr, hmat, pmat, smat)
   end subroutine construct_error_dp

   subroutine construct_error_sp(self, hmat, pmat, smat)
      use mctc_env, only : sp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(sp), intent(in) :: hmat(:,:)
      !> Density matrix
      real(sp), intent(in) :: pmat(:,:)
      !> Overlap matrix
      real(sp), intent(in) :: smat(:,:)

      call construct_error_vec_sp(self%ptr, hmat, pmat, smat)
   end subroutine construct_error_sp

   !> Set an error vector manually to be used for mixing
   subroutine set_error_dp(self, error)
      use mctc_env, only : dp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Error matrix
      real(dp), intent(in) :: error(:,:)

      call set_error_vec_dp(self%ptr, error, size(error(:,1)), size(error(1,:)))
   end subroutine set_error_dp

   subroutine set_error_sp(self, error)
      use mctc_env, only : sp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Error matrix
      real(sp), intent(in) :: error(:,:)

      call set_error_vec_sp(self%ptr, error, size(error(:,1)), size(error(1,:)))
   end subroutine set_error_sp

   !> Set the vector to mix
   subroutine set_diis_dp(self, hmat)
      use mctc_env, only : dp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(dp), intent(in) :: hmat(:,:)
      
      call set_mixer_data_dp(self%ptr, hmat)
   end subroutine set_diis_dp

   subroutine set_diis_sp(self, hmat)
      use mctc_env, only : sp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(sp), intent(in) :: hmat(:,:)
  
      call set_mixer_data_sp(self%ptr, hmat)
   end subroutine set_diis_sp

   !> Get the differences of the mixed vector
   subroutine diff_diis_dp(self, hmat)
      use mctc_env, only : dp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(dp), intent(in) :: hmat(:,:)

      call diff_mixer_data_dp(self%ptr, hmat, size(hmat))
   end subroutine diff_diis_dp

   subroutine diff_diis_sp(self, hmat)
      use mctc_env, only : sp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(sp), intent(in) :: hmat(:,:)

      call diff_mixer_data_sp(self%ptr, hmat, size(hmat))
   end subroutine diff_diis_sp

   !> Get the mixed vector
   subroutine get_diis_dp(self, hmat)
      use mctc_env, only : dp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(dp), intent(out) :: hmat(:,:)

      call get_mixer_data_dp(self%ptr, hmat, size(hmat))
   end subroutine get_diis_dp

   subroutine get_diis_sp(self, hmat)
      use mctc_env, only : sp
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Fock matrix
      real(sp), intent(out) :: hmat(:,:)

      call get_mixer_data_sp(self%ptr, hmat, size(hmat))
   end subroutine get_diis_sp

end module tblite_scf_mixer_diis
