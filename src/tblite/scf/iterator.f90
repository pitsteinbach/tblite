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

!> @file tblite/scf/iterator.f90
!> Provides the implementation of the actual self-consistent field iteractions

!> Iterator for evaluating the Hamiltonian self-consistently
module tblite_scf_iterator
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache, container_list
   use tblite_disp, only : dispersion_type
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
   & get_mulliken_atomic_multipoles
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_scf_mixer_input, only : mixer_type
   use tblite_scf_mixer_type, only : mixers_type
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type, add_pot_to_h1
   use tblite_scf_utils, only : get_density, get_electronic_energy, reduce, get_qat_from_qsh
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: next_scf

contains

!> Evaluate self-consistent iteration for the density-dependent Hamiltonian
subroutine next_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, ccache, dcache, icache, energies, error)
   !> Current iteration count
   integer, intent(inout) :: iscf
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Convergence accelerator
   class(mixers_type), intent(inout) :: mixer
   !> Information on wavefunction data used to construct Hamiltonian
   type(scf_info), intent(in) :: info
   !> Container for coulombic interactions
   type(tb_coulomb), intent(in), optional :: coulomb
   !> Container for dispersion interactions
   class(dispersion_type), intent(in), optional :: dispersion
   !> Container for general interactions
   type(container_list), intent(in), optional :: interactions

   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Restart data for coulombic interactions
   type(container_cache), intent(inout), optional :: ccache
   !> Restart data for dispersion interactions
   type(container_cache), intent(inout), optional :: dcache
   !> Restart data for interaction containers
   type(container_cache), intent(inout), optional :: icache

   !> Self-consistent energy
   real(wp), intent(inout) :: energies(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: eao(:)
   real(wp) :: ts

   if (iscf > 0 .and. (mixer%type(1) == mixer_type%broyden .or. mixer%type(1) == mixer_type%gambits_broyden)) then
      call mixer%next_mixer(iscf, wfn, error)
      if (allocated(error)) return
      call mixer%get_mixer(bas, wfn)
   end if

   iscf = iscf + 1
   call pot%reset
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_potential(mol, ccache, wfn, pot)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_potential(mol, dcache, wfn, pot)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_potential(mol, icache, wfn, pot)
   end if
   call add_pot_to_h1(bas, ints, pot, wfn%coeff)

   call mixer%set_mixer(wfn)

   if (mixer%type(1) == mixer_type%gambits_diis .and. iscf > 1) then
      call mixer%next_mixer(iscf, wfn, error)
      if (allocated(error)) return
      call mixer%get_mixer(bas, wfn)
   end if

   call get_density(wfn, solver, ints, ts, error)
   if (allocated(error)) return

   call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
   call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)

   call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)
   call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
      & wfn%qpat)

   call mixer%diff_mixer(wfn)

   allocate(eao(bas%nao), source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, eao)

   energies(:) = ts / size(energies)
   call reduce(energies, eao, bas%ao2at)
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_energy(mol, ccache, wfn, energies)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if
end subroutine next_scf

end module tblite_scf_iterator
