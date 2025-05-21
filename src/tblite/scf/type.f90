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

!> @file tblite/scf/mixers.f90
!> Mixers base class

!> Provides an electronic mixer implementation

module tblite_scf_mixers
    use mctc_env, only: wp
    use tblite_scf_mixer_broyden, only: broyden_type, new_broyden
    use tblite_scf_mixer_diis, only: diis_type, new_diis
 
    use tblite_xtb_calculator, only : xtb_calculator
    use mctc_io, only : structure_type
    use tblite_wavefunction_type, only : wavefunction_type
    use tblite_scf_info, only : scf_info
    use iso_c_binding
 
    type, public :: mixers_type
       !> Broyden mixer object
       type(broyden_type), allocatable :: broyden
       !> DIIS mixer object
       type(diis_type), allocatable :: diis
       !> Pointer to the current mixer
       type(c_ptr) :: currptr
 
    contains
       !> Create mixers
       procedure :: setup
       !> Destroy mixers
       procedure :: destroy
    end type
 
    interface
       subroutine destroy_mixer(mixer) bind(C,name="Destroy")
          use iso_c_binding
          type(c_ptr), value :: mixer
       end subroutine destroy_mixer
    end interface
 
 contains
    subroutine setup(self, mixer_type, mol, calc, wfn, info)
       !> Mixers object
       class(mixers_type), intent(inout) :: self
       !> Type(s) of self-consistent iteration mixing (0: Broyden, 1: DIIS)
       integer, intent(in) :: mixer_type
       !> Molecular structure data
       type(structure_type), intent(in) :: mol
       !> Single-point calculator
       type(xtb_calculator), intent(in) :: calc
       !> Tight-binding wavefunction data
       type(wavefunction_type), intent(in) :: wfn
       !> Info data
       type(scf_info) :: info
 
       select case (mixer_type)
        case(0)
         allocate(self%broyden)
         call new_broyden(self%broyden, mol, calc, wfn, info)
         self%currptr=self%broyden%ptr
        case(1)
         allocate(self%diis)
          call new_diis(self%diis, mol, calc, info)
          self%currptr=self%diis%ptr
       end select
 
 end subroutine setup
 
 
 subroutine destroy(self)
    !> Mixer object
    class(mixers_type), intent(inout) :: self
 
    if (allocated(self%broyden)) then
       call destroy_mixer(self%broyden%ptr)
     endif
 
    if (allocated(self%diis)) then
       call destroy_mixer(self%diis%ptr)
    endif
 
 end subroutine destroy
 
 end module tblite_scf_mixers
 