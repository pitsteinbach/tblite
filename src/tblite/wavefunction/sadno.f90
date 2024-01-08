module tblite_wavefunction_guess_sadno
   use mctc_env, only : wp
   use mctc_io, only : structure_type, new
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_lapack_solver, only : lapack_solver
   use tblite_context, only : context_type, context_terminal
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_wavefunction_guess, only : sad_guess
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_calculator, only : new_xtb_calculator
   use tblite_param, only : param_record
   use mctc_env, only : error_type
   implicit none
   private
   public :: sadno_guess
   real(wp), parameter :: kt = 3.166808578545117e-06_wp
contains
    
   subroutine sadno_guess(mol, calc, wfn, method, param)
      type(structure_type), intent(in) :: mol
      type(xtb_calculator), intent(in) :: calc
      type(wavefunction_type), intent(inout) :: wfn
      type(xtb_calculator) :: calc_atom
      type(structure_type) :: atom
      type(wavefunction_type) :: wfn_atom
      character(len=*), optional :: method
      type(param_record), optional :: param
      type(error_type), allocatable :: error
      real(kind=wp) :: xyz(3, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      real(kind=wp), allocatable :: atomic_densities(:,:,:)
      integer :: num(1)
      integer, allocatable :: nao_id(:), nao_at(:)
      integer :: maxnao, i, ish, id, offset
      real(kind = wp) :: en
      type(context_type) :: context
      !set up calc object for each id type of atom in an analogue fashion to the molecular setting 
      !iterate over all store atomic ref densities
      allocate(nao_at(mol%nat), source=0)
      associate(bas => calc%bas)
      do ish = 1, bas%nsh
         nao_at(bas%sh2at(ish)) = nao_at(bas%sh2at(ish)) + bas%nao_sh(ish)
      end do
      id = 0
      allocate(nao_id(mol%nid), source=0)
      do i = 1, mol%nat
         if (mol%id(i) .ne. id) then
            nao_id(mol%id(i)) = nao_at(i)
            id = mol%id(i)
         end if
      end do
      end associate
      maxnao = maxval(nao_id)
      allocate(atomic_densities(mol%nid, maxnao, maxnao), source = 0.0_wp)
      context%solver = lapack_solver(1)
      context%terminal = context_terminal(.false.)
      context%verbosity = 1
      en = 0.0_wp
      do i = 1, mol%nid
         num(1) = mol%num(i) 
         call new(atom, num, xyz, charge=0.0_wp, uhf=0)
         call new_wavefunction(wfn_atom, 1, calc%bas%nsh_id(i), nao_id(i), 1, 300 * kt)
         if (present(method)) then
            select case(method)
            case("gfn2")
               call new_gfn2_calculator(calc_atom, atom)
            case("gfn1")
               call new_gfn1_calculator(calc_atom, atom)
            case("ipea1")
               call new_ipea1_calculator(calc_atom, atom)
            end select
         else if (present(param)) then
            call new_xtb_calculator(calc_atom, atom, param, error)
         end if
         call sad_guess(atom, calc_atom, wfn_atom)
         call xtb_singlepoint(context, atom, calc_atom, wfn_atom, 1.0_wp, en)
         write(*,*) en
         atomic_densities(i, : , :) = wfn_atom%density(:, :, 1)
      enddo
      offset = 0
      do i = 1, mol%nat
         wfn%density(offset+1:offset+nao_at(i), offset+1:offset+nao_at(i), 1) = atomic_densities(mol%id(i),1:nao_at(i), 1:nao_at(i)) 
         offset = offset + nao_at(i)
      end do
      write(*,*) wfn%density(1, 1, 1)
   end subroutine
end module