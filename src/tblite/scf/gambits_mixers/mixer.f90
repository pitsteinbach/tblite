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
module tblite_scf_gambits_mixer
   use iso_c_binding
   implicit none



   !> Electronic mixer
   type, public, abstract :: gambits_mixer_type
      integer :: ndim
      integer :: memory
      integer :: start
      type(c_ptr) :: ptr
   contains
      !> Apply mixing to the density
      procedure :: next
   end type gambits_mixer_type

   interface
   subroutine next_mixer(mixer,iter) bind(C,name="Next")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: mixer
      integer(c_int), value, intent(in) :: iter
   end subroutine next_mixer
   end interface
 
   contains

   subroutine next(self, iter)
      !> Mixer object
      class(gambits_mixer_type), intent(inout) :: self
      !> SCF Iteration
      integer, intent(in) :: iter

      call next_mixer(self%ptr, iter)
   end subroutine next

end module tblite_scf_gambits_mixer
