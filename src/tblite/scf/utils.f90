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

!> @file tblite/scf/utils.f90
!> Provides the implementation of the actual self-consistent field iteractions

!> Helper modules for the self-consistent field iteractions
module tblite_scf_utils
   use mctc_env, only : wp, error_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_type, only : wavefunction_type, get_density_matrix
   use tblite_wavefunction_fermi, only : get_fermi_filling
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: get_electronic_energy, reduce
   public :: next_density, get_qat_from_qsh

contains

subroutine get_electronic_energy(h0, density, energies)
   real(wp), intent(in) :: h0(:, :)
   real(wp), intent(in) :: density(:, :, :)
   real(wp), intent(inout) :: energies(:)

   integer :: iao, jao, spin

   !$omp parallel do collapse(3) schedule(runtime) default(none) &
   !$omp reduction(+:energies) shared(h0, density) private(spin, iao, jao)
   do spin = 1, size(density, 3)
      do iao = 1, size(density, 2)
         do jao = 1, size(density, 1)
            energies(iao) = energies(iao) + h0(jao, iao) * density(jao, iao, spin)
         end do
      end do
   end do
end subroutine get_electronic_energy

subroutine reduce(reduced, full, map)
   real(wp), intent(inout) :: reduced(:)
   real(wp), intent(in) :: full(:)
   integer, intent(in) :: map(:)

   integer :: ix

   do ix = 1, size(map)
      reduced(map(ix)) = reduced(map(ix)) + full(ix)
   end do
end subroutine reduce

subroutine get_qat_from_qsh(bas, qsh, qat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: qsh(:, :)
   real(wp), intent(out) :: qat(:, :)

   integer :: ish, ispin

   qat(:, :) = 0.0_wp
   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:qat) shared(bas, qsh) private(ish)
   do ispin = 1, size(qsh, 2)
      do ish = 1, size(qsh, 1)
         qat(bas%sh2at(ish), ispin) = qat(bas%sh2at(ish), ispin) + qsh(ish, ispin)
      end do
   end do
end subroutine get_qat_from_qsh

subroutine next_density(wfn, solver, ints, ts, error)
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Electronic entropy
   real(wp), intent(out) :: ts
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi, stmp(2)
   real(wp), allocatable :: focc(:)
   integer :: spin

   call solver%get_density(wfn%coeff, ints%overlap, wfn%emo, wfn%focc, wfn%density, error)
   do spin = 1, 2
      call get_electronic_entropy(wfn%focc(:, spin), wfn%kt, stmp(spin))
   end do
   ts = sum(stmp)
end subroutine next_density

subroutine get_electronic_entropy(occ, kt, s)
   real(wp), intent(in) :: occ(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: s

   s = sum(log(occ ** occ * (1 - occ) ** (1 - occ))) * kt
end subroutine get_electronic_entropy

end module tblite_scf_utils
