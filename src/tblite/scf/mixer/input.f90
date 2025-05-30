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

!> @file tblite/scf/mixer/input.f90
!> Provides an electronic mixer implementation

!> Module for processing mixer input parameters
module tblite_scf_mixer_input
   use mctc_env, only : wp
   implicit none
   private

   !> Configuration for the Broyden mixer
   type, public :: broyden_input
      !> Number of steps to keep in memory
      integer :: memory = 250
      !> Damping parameter
      real(wp) :: damp = 0.4_wp
   end type broyden_input

   !> Configuration for the GAMBITS Broyden mixer
   type, public :: gambits_broyden_input
      !> Number of steps to keep in memory
      integer :: memory = 250
      !> Damping parameter
      real(wp) :: damp = 0.4_wp
      !> Number of atomic orbitals
      integer :: nao
   end type gambits_broyden_input

   !> Configuration for the GAMBITS DIIS mixer
   type, public :: gambits_diis_input
      !> Number of steps to keep in memory
      integer :: memory = 5
      !> Damping parameter
      real(wp) :: damp = 0.4_wp
      !> Size-dependent (0), CPU (1), or GPU (2) runmode
      integer :: runmode = 0
      !> Working precision (FP32: 0, FP64: 1)
      integer :: prec = 1
      !> Number of atomic orbitals
      integer :: nao
   end type gambits_diis_input

   !> Input parameters for electronic mixer
   type, public :: mixer_input
      !> Mixer type
      integer :: type
      !> Input for Broyden mixer
      type(broyden_input) :: broyden
      !> Input for GAMBITS Broyden mixer
      type(gambits_broyden_input) :: gambits_broyden
      !> Input for GAMBITS DIIS mixer
      type(gambits_diis_input) :: gambits_diis
   end type mixer_input

end module tblite_scf_mixer_input
