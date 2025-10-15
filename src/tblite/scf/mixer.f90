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
   use gambits, only : gambits_context_type
   use mctc_env, only : wp, sp, dp, error_type
   use mctc_io, only : structure_type
   use tblite_basis, only : basis_type
   use tblite_context, only : context_type
   use tblite_scf_gambits_broyden, only : gambits_broyden_type, new_gambits_broyden
   use tblite_scf_gambits_diis, only : gambits_diis_type, new_gambits_diis
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, &
   & not_used, orbital_resolved
   use tblite_scf_mixer_broyden, only : broyden_mixer, new_broyden
   use tblite_scf_mixer_input, only : mixer_input, mixer_kind
   use tblite_scf_mixer_type, only : mixers_type
   use iso_c_binding, only : c_size_t
   implicit none
   private

   public :: new_mixer, get_mixer_dimension

contains

!> Create a new instance of the mixer
subroutine new_mixer(self, ctx, input, ndim, nao, nspin, overlap, info, error)
   !> Instance of the mixer on exit
   class(mixers_type), intent(out) :: self
   !> Calculation context
   class(context_type), intent(inout) :: ctx
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
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: prec, i

   if (wp == sp) then
      prec = 0
   else if (wp == dp) then
      prec = 1
   end if

   select case(input%kind)
      case(mixer_kind%broyden)
      block
         type(broyden_mixer), allocatable :: mixer
         allocate(mixer)
         call new_broyden(mixer, nspin*ndim, input)
         mixer%info = info
         allocate(self%mixer(1),source=mixer)
      end block
      allocate(self%kind(1))
      self%kind = mixer_kind%broyden

      case(mixer_kind%gambits_broyden)
      block
         type(gambits_broyden_type), allocatable :: mixer
         allocate(mixer)
         call new_gambits_broyden(mixer, nspin*ndim, input%memory(input%kind), input%damp, nao, prec, error)
         mixer%info = info
         allocate(self%mixer(1),source=mixer)
      end block
      if (len(self%mixer(1)%msg) > 0) call ctx%message(self%mixer(1)%msg)
      if (allocated(self%mixer(1)%msg)) deallocate(self%mixer(1)%msg) 
      
      allocate(self%kind(1))
      self%kind = mixer_kind%gambits_broyden

      case(mixer_kind%gambits_diis)
      block
         type(gambits_diis_type), allocatable :: mixer(:)
         allocate(mixer(nspin))
         do i=1,nspin
            call new_gambits_diis(mixer(i), nao**2, input%memory(input%kind), input%damp, overlap, nao, &
               & input%runmode, prec, input%prec, error)
            if (allocated(error)) return
            if (len(mixer(i)%msg) > 0) call ctx%message(mixer(i)%msg)
            if (allocated(mixer(i)%msg)) deallocate(mixer(i)%msg) 
            call mixer(i)%diis_info(info)
            mixer(i)%nspin = nspin
            mixer(i)%channel = i
         end do
         allocate(self%mixer,source=mixer)
      end block
      
      allocate(self%kind(nspin))
      self%kind = mixer_kind%gambits_diis
   end select

end subroutine new_mixer

function get_mixer_dimension(mol, bas, info) result(ndim)
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