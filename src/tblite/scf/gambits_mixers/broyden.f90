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

!> @file tblite/scf/mixers/broyden.f90
!> Implementing Broyden mixing

!> Provides an electronic mixer implementation
module tblite_scf_gambits_broyden
   use tblite_scf_gambits_mixer
   use tblite_basis_type, only : basis_type
   use tblite_xtb_calculator, only : xtb_calculator
   use mctc_io, only : structure_type
   use tblite_wavefunction, only : wavefunction_type
   use tblite_scf_info, only : scf_info
   use iso_c_binding
   implicit none

!> Broyden mixer
   type, extends(gambits_mixer_type) :: gambits_broyden_type

   contains
      !> Set new object to mix
      procedure :: set => set_broyden_dp, set_broyden_sp
      !> Set difference between two consecutive objects to mix
      procedure :: diff => diff_broyden_dp, diff_broyden_sp
      !> Get mixed object
      procedure :: get => get_broyden_dp, get_broyden_sp
   end type gambits_broyden_type

   interface
      type(c_ptr) function c_new_broyden(ndim, memory, alpha, nao) bind(C,name="SetupBroyden")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
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
   end interface

contains

   !> Create a new instance of the Broyden mixer
   subroutine new_broyden(self, mol, calc, wfn, info)
      use tblite_scf_mixer, only : get_mixer_dimension
      !> Broyden object
      class(gambits_broyden_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
      !> Tight-binding wavefunction data
      type(wavefunction_type), intent(in) :: wfn
      !> Info data
      type(scf_info) :: info

      self%ndim = wfn%nspin * get_mixer_dimension(mol,calc%bas,info)
      self%memory = calc%mixer_mem
      self%ptr = c_new_broyden(self%ndim, self%memory, calc%mixer_damping, calc%bas%nao)
   end subroutine new_broyden

   !> Set the vector to mix
   subroutine set_broyden_dp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : dp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(dp), intent(in) :: qat(:)
      !> Shell charges
      real(dp), intent(in) :: qsh(:)
      !> Dipole moments
      real(dp), intent(in) :: dpat(:,:)
      !> Quadrupole moments
      real(dp), intent(in) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info

         select case(info%charge)
         case(atom_resolved)
            call set_mixer_data_dp(self%ptr, qat, size(qat))
         case(shell_resolved)
            call set_mixer_data_dp(self%ptr, qsh, size(qsh))
         end select

         select case(info%dipole)
         case(atom_resolved)
            call set_mixer_data_dp(self%ptr, dpat, size(dpat))
         end select

         select case(info%quadrupole)
         case(atom_resolved)
            call set_mixer_data_dp(self%ptr, qpat, size(qpat))
         end select
   end subroutine set_broyden_dp

   subroutine set_broyden_sp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : sp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(sp), intent(in) :: qat(:)
      !> Shell charges
      real(sp), intent(in) :: qsh(:)
      !> Dipole moments
      real(sp), intent(in) :: dpat(:,:)
      !> Quadrupole moments
      real(sp), intent(in) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info

      select case(info%charge)
       case(atom_resolved)
         call set_mixer_data_sp(self%ptr, qat, size(qat))
       case(shell_resolved)
         call set_mixer_data_sp(self%ptr, qsh, size(qsh))
      end select

      select case(info%dipole)
       case(atom_resolved)
         call set_mixer_data_sp(self%ptr, dpat, size(dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call set_mixer_data_sp(self%ptr, qpat, size(qpat))
      end select
   end subroutine set_broyden_sp


!> Get the differences of the mixed vector
   subroutine diff_broyden_dp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : dp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(dp), intent(in) :: qat(:)
      !> Shell charges
      real(dp), intent(in) :: qsh(:)
      !> Dipole moments
      real(dp), intent(in) :: dpat(:,:)
      !> Quadrupole moments
      real(dp), intent(in) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info

      select case(info%charge)
       case(atom_resolved)
         call diff_mixer_data_dp(self%ptr, qat, size(qat))
       case(shell_resolved)
         call diff_mixer_data_dp(self%ptr, qsh, size(qsh))
      end select

      select case(info%dipole)
       case(atom_resolved)
         call diff_mixer_data_dp(self%ptr, dpat, size(dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call diff_mixer_data_dp(self%ptr, qpat, size(qpat))
      end select
   end subroutine diff_broyden_dp

   subroutine diff_broyden_sp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : sp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(sp), intent(in) :: qat(:)
      !> Shell charges
      real(sp), intent(in) :: qsh(:)
      !> Dipole moments
      real(sp), intent(in) :: dpat(:,:)
      !> Quadrupole moments
      real(sp), intent(in) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info

      select case(info%charge)
       case(atom_resolved)
         call diff_mixer_data_sp(self%ptr, qat, size(qat))
       case(shell_resolved)
         call diff_mixer_data_sp(self%ptr, qsh, size(qsh))
      end select

      select case(info%dipole)
       case(atom_resolved)
         call diff_mixer_data_sp(self%ptr, dpat, size(dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call diff_mixer_data_sp(self%ptr, qpat, size(qpat))
      end select
   end subroutine diff_broyden_sp

!> Get the mixed vector
   subroutine get_broyden_dp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : dp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(dp), intent(out) :: qat(:)
      !> Shell charges
      real(dp), intent(out) :: qsh(:)
      !> Dipole moments
      real(dp), intent(out) :: dpat(:,:)
      !> Quadrupole moments
      real(dp), intent(out) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info
      
      select case(info%charge)
         case(atom_resolved)
            call get_mixer_data_dp(self%ptr, qat, size(qat))
         case(shell_resolved)
            call get_mixer_data_dp(self%ptr, qsh, size(qsh))
      end select

      select case(info%dipole)
         case(atom_resolved)
            call get_mixer_data_dp(self%ptr, dpat, size(dpat))
      end select

      select case(info%quadrupole)
         case(atom_resolved)
            call get_mixer_data_dp(self%ptr, qpat, size(qpat))
      end select
   end subroutine get_broyden_dp

   subroutine get_broyden_sp(self, qat, qsh, dpat, qpat, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      use mctc_env, only : sp
      !> Instance of the Broyden mixer
      class(gambits_broyden_type), intent(inout) :: self
      !> Atom charges
      real(sp), intent(out) :: qat(:)
      !> Shell charges
      real(sp), intent(out) :: qsh(:)
      !> Dipole moments
      real(sp), intent(out) :: dpat(:,:)
      !> Quadrupole moments
      real(sp), intent(out) :: qpat(:,:)
      !> Info data
      type(scf_info) :: info

      select case(info%charge)
         case(atom_resolved)
            call get_mixer_data_sp(self%ptr, qat, size(qat))
         case(shell_resolved)
            call get_mixer_data_sp(self%ptr, qsh, size(qsh))
      end select

      select case(info%dipole)
         case(atom_resolved)
            call get_mixer_data_sp(self%ptr, dpat, size(dpat))
      end select

      select case(info%quadrupole)
         case(atom_resolved)
            call get_mixer_data_sp(self%ptr, qpat, size(qpat))
      end select
   end subroutine get_broyden_sp

end module tblite_scf_gambits_broyden
