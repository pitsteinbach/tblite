module tblite_post_processing_bond_orders
    use mctc_env, only : wp
    use tblite_post_processing_type, only : post_processing_type
    use tblite_wavefunction_type, only : wavefunction_type, get_density_matrix
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_container, only : container_cache
    use tblite_results, only : results_type
    use tblite_context, only : context_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_double_dictionary, only : double_dictionary_type
    use tblite_wavefunction_spin, only : updown_to_magnet
    use tblite_timer, only : timer_type, format_time
    use tblite_output_format, only : format_string
    use tblite_blas, only : gemm
    implicit none
    private

    public :: new_wbo, wiberg_bond_orders

    type, extends(post_processing_type) :: wiberg_bond_orders
    contains
        procedure :: compute
        procedure :: print_timer
    end type

    character(len=24), parameter :: label = "Mayer-Wiberg bond orders"
    type(timer_type) :: timer
contains

subroutine new_wbo(new_wbo_type)
    type(wiberg_bond_orders), intent(inout) :: new_wbo_type
        new_wbo_type%label = label
end subroutine

subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
        class(wiberg_bond_orders),intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction strcuture data
        type(wavefunction_type), intent(in) :: wfn
        !> integral container for dipole and quadrupole integrals for CAMMs
        type(integral_type) :: integrals
        !> Single-point calculator conatiner
        type(xtb_calculator), intent(in) :: calc
        !> Context container for writing to stdout
        type(context_type), intent(inout) :: ctx
        type(container_cache), intent(inout) :: cache_list(:)
        type(double_dictionary_type), intent(inout) :: dict
        real(kind=wp), allocatable :: wbo(:, :, :)
        real(kind=wp), allocatable :: focc_(:, :)
        integer :: prlevel, nspin, i, j
        real(kind=wp) :: nel_
        real(wp), allocatable :: pmat(:, :, :)

        call timer%push("total")
        nspin = size(wfn%density, dim=3)
        
        if ((nspin == 1) .and. (wfn%nel(1) /= wfn%nel(2))) then
            allocate(wbo(mol%nat, mol%nat, 2), source=0.0_wp)
            allocate(focc_(calc%bas%nao, 2), source=0.0_wp)
            allocate(pmat(calc%bas%nao, calc%bas%nao, 2))
            do j = 1,2
                nel_ = wfn%nel(j)
                do i = 1, size(wfn%focc)
                    if (nel_ > 1.0_wp) then
                        focc_(i,j) = 1.0_wp
                        nel_ = nel_ - 1.0_wp
                    else
                        focc_(i,j) = nel_
                        exit
                    end if
                end do
            end do
            do j = 1, 2
                call get_density_matrix(focc_(:, j), wfn%coeff(:, :, nspin), pmat(:, :, j))
            end do
            call get_mayer_bond_orders(calc%bas, integrals%overlap, pmat, wbo)
        else
            allocate(wbo(mol%nat, mol%nat, nspin), source=0.0_wp)
            call get_mayer_bond_orders(calc%bas, integrals%overlap, wfn%density, wbo)
        end if
        
        call dict%add_entry("wbo", wbo(:, :, :))
        call timer%pop()
end subroutine

subroutine get_mayer_bond_orders(bas, smat, pmat, mbo)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Wiberg/Mayer bond orders
   real(wp), intent(out) :: mbo(:, :, :)

   integer :: iao, jao, iat, jat, spin
   real(wp) :: pao
   real(wp), allocatable :: psmat(:, :)

   allocate(psmat(bas%nao, bas%nao))

   mbo(:, :, :) = 0.0_wp
   do spin = 1, size(pmat, 3)
      call gemm(pmat(:, :, spin), smat, psmat)
      !$omp parallel do default(none) collapse(2) &
      !$omp shared(bas, psmat, mbo, spin) private(iao, jao, iat, jat, pao)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            iat = bas%ao2at(iao)
            jat = bas%ao2at(jao)
            pao = merge(psmat(iao, jao) * psmat(jao, iao), 0.0_wp, iat /= jat)
            !$omp atomic
            mbo(jat, iat, spin) = mbo(jat, iat, spin) + pao
         end do
      end do
   end do

   call updown_to_magnet(mbo)
end subroutine get_mayer_bond_orders

subroutine print_timer(self, prlevel, ctx)
    !> Instance of the interaction container
    class(wiberg_bond_orders), intent(in) :: self
    integer :: prlevel
    type(context_type) :: ctx
    real(wp) :: ttime, stime
    integer :: it
    character(len=*), parameter :: labels(*) = [character(len=20):: &
         & ]

    

    if (prlevel > 2) then
        call ctx%message(label//" timing details:")
        ttime = timer%get("total")
        call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))        
        do it = 1, size(labels)
            stime = timer%get(labels(it))
            if (stime <= epsilon(0.0_wp)) cycle
            call ctx%message(" - "//labels(it)//format_time(stime) &
                & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
        end do
        call ctx%message("")
    end if
end subroutine

end module