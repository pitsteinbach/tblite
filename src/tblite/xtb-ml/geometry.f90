module tblite_xtbml_geometry_based
    implicit none
    use mctc_env, only : wp
    use tblite_xtbml_feature_type, only : xtbml_feature_type
    private
    integer :: n_features
    
type :: enum_geometry_features
    !> atomic coordination number 
    integer :: CN = 1
    !> extended coordination number
    integer :: delta_CN = 2
end type

type, public, extends(xtbml_feature_type) :: xtbml_geometry_features_type
    character(len=*) :: label = "geometry-based features"
    real(wp), allocatable ::  cn_atom(:)
    real(wp), allocatable ::  delta_cn(:, :)
    type(enum_geometry_features) :: labels = enum_geometry_features()

contains
    procedure :: compute_features
    procedure :: compute_extended
    procedure :: get_n_features
    procedure :: get_feature_labels
end type

contains

subroutine compute_features(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
    use tblite_ncoord_exp, only : new_exp_ncoord, exp_ncoord_type
    class(xtbml_type), intent(inout) :: self
    !> Molecular structure data
    type(structure_type), intent(in) :: mol
    !> Wavefunction strcuture data
    type(wavefunction_type), intent(in) :: wfn
    !> Integral container
    type(integral_type) :: integrals
    !> Single-point calculator
    type(xtb_calculator), intent(in) :: calc
    !> Container
    type(container_cache), intent(inout) :: cache
    !> Context type
    type(context_type),intent(inout) :: ctx
    !> Print Level
    integer, intent(in) :: prlevel
    n_features = n_features + 1
    allocate(self%cn_atom(mol%nat))
    
    call new_exp_ncoord(ncoord_exp, mol)
    call ncoord_exp%get_cn(mol, self%cn_atom) 
    
end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
    class(xtbml_type), intent(inout) :: self
    !> Molecular structure data
    type(structure_type), intent(in) :: mol
    !> Wavefunction strcuture data
    type(wavefunction_type), intent(in) :: wfn
    !> Integral container
    type(integral_type) :: integrals
    !> Single-point calculator
    type(xtb_calculator), intent(in) :: calc
    !> Container
    type(container_cache), intent(inout) :: cache
    !> Context type
    type(context_type),intent(inout) :: ctx
    !> Print Level
    integer, intent(in) :: prlevel
    if (allocated(calc%array)) then 
        self%a = calc%a_array
    else 
        self%a = [1.0]
    end if
    if (.not.allocated(self%inv_cn_a)) then
        call self%populate_inv_cn_array(mol%nat, mol%at, mol%xyz)
    endif
    call get_delta_cn(mol%nat, n, self%cn_atom, mol%id, mol%xyz, self%delta_cn)

end subroutine

subroutine get_delta_cn(nat, n_a, cn, at, xyz, delta_cn)
    use tblite_timer, only : timer_type, format_time
    integer, intent(in) :: nat, at(nat), n_a
    real(wp), intent(in) :: cn(nat), xyz(3, nat)
    real(wp), intent(out) :: delta_cn(nat, n_a)
    real(wp):: delta_cn_tmp(nat, n_a), stime
    integer :: i, j, k
    type(timer_type) :: timer
    
    delta_cn = 0.0_wp
    !$omp parallel do default(none) collapse(2)&
    !$omp shared(nat, self%inv_cn_a, n_a)&
    !$omp private(delta_cn,i,j,k)
    do k = 1, n_a
       do i = 1, nat
          do j = 1, nat
             if (i == j) cycle
             !$omp atomic
             delta_cn(i, k) = delta_cn(i, k) + cn(j)/self%inv_cn_a(i, j, k)
 
          end do
       end do
    end do
    !$omp end parallel do
 end subroutine





end module tblite_xtbml_geometry_based   