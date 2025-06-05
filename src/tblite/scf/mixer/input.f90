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

   public :: mixer_type, mixer_precision, mixer_runmode

type enum_mixers
   !> Native broyden mixer
   integer :: broyden = 1
   !> Gambits broyden mixer
   integer :: gambits_broyden = 2
   !> Gambits DIIS mixer
   integer :: gambits_diis = 3
end type enum_mixers

type enum_precision
   !> Single precision (FP32)
   integer :: single = 0
   !> Double precision (FP64)
   integer :: double = 1
end type enum_precision

type enum_runmode
   !> Default runmode (size-dependent)
   integer :: default = 0
   !> CPU runmode
   integer :: cpu = 1
   !> GPU runmode
   integer :: gpu = 2
end type enum_runmode

   type(enum_mixers), parameter :: mixer_type = enum_mixers()
   type(enum_precision), parameter :: mixer_precision = enum_precision()
   type(enum_runmode), parameter :: mixer_runmode = enum_runmode()

   !> Input parameters for electronic mixer
   type, public :: mixer_input
      !> Mixer type
      integer :: type = mixer_type%broyden
      !> Damping parameter
      real(wp) :: damp = 0.4_wp
      !> Number of steps to keep in memory
      integer, dimension(3) :: memory = [250, 250, 5] ! Defaults in order of enumerator
      !> Number of atomic orbitals
      integer :: nao
      !> DIIS Working precision (FP32: 0, FP64: 1)
      integer :: prec = mixer_precision%single
      !> DIIS runmode (size-dependent: 0, cpu: 1, gpu: 2)
      integer :: runmode = mixer_runmode%default
   end type mixer_input

end module tblite_scf_mixer_input
