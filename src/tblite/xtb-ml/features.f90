module tblite_xtbml_feature_type
    implicit none
    use mctc_env, only : wp
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_res ults, only : results_type
   use tblite_integral_type, only : integral_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   private

type, public, abstract :: xtbml_feature_type
    real(wp), allocatable :: rcov(:)
    real(wp), allocatable :: inv_cn_a(:, :, :)
    real(wp), allocatable :: a(:)
    real(wp) :: k1 = 16.0_wp

contains
    procedure(compute_features), deferred :: compute_features
    procedure(compute_extended), deferred :: compute_extended
    procedure(get_n_features), deferred :: get_n_features
    procedure, private :: populate_inv_cn_array
end type xtbml_feature_type

abstract interface
    subroutine compute_features(self, mol, wfn, integrals, calc, cache, prlevel, ctx)
        import :: wp, wavefunction_type, structure_type, wavefunction_type, integral_type, xtb_calculator,&
        & container_cache, context_type, xtbml_feature_type
        class(xtbml_feature_type), intent(inout) :: self
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
    end subroutine
    subroutine get_n_features(self,n)
        import :: xtbml_feature_type
        class(xtbml_feature_type), intent(inout) :: self
        !> Number of features in feature type
        integer :: n
    end subroutine
end interface

contains

subroutine get_rcov(mol)
    use tblite_data_covrad, only : get_covalent_rad
    type(structure_type), intent(in) :: mol
    if (allocated(rcov)) then
       deallocate (rcov)
    end if
    allocate (rcov(mol%nid), source=0.0_wp)
    rcov(:) = get_covalent_rad(mol%num)
end subroutine

subroutine populate_inv_cn_array(nat, at, xyz)
    integer, intent(in):: nat, at(nat)
    real(wp), intent(in) :: xyz(:, :)
    real(wp) :: result
    integer :: i, j, k
    n_a = size(self%a)
    if (allocated(self%rcov)) then
        deallocate (slef%rcov)
    end if
    allocate (self%inv_rcov(nat), source=0.0_wp)
    call get_rcov(mol)
    if (allocated(self%inv_cn_a)) then
       deallocate (slef%inv_cn_a)
    end if
    allocate (self%inv_cn_a(nat, nat, n_a), source=0.0_wp)
    !$omp parallel do default(none) collapse(2)&
    !$omp shared(self%a,nat, at,xyz,self%inv_cn_a,n_a)&
    !$omp private(result,i,j,k)
    do k = 1, n_a
       do i = 1, nat
          do j = 1, nat
             !if (i == j) cycle
             call inv_cn(nat, i, j, at, xyz, self%a(k), result)
             self%inv_cn_a(i, j, k) = result
          end do
       end do
    end do
    !$omp end parallel do
 
 end subroutine populate_inv_cn_array

subroutine inv_cn(nat, a, b, at, xyz, dampening_fact, result)
   
    integer, intent(in) :: a, b, nat
    integer, intent(in) :: at(:)
    real(wp), intent(in)  :: xyz(3, nat), dampening_fact
    real(wp), intent(out) :: result
    real(wp) :: rab(3), r, rco, r2
 
    result = 0.0_wp
 
    rab = xyz(1:3, a) - xyz(1:3, b)
    r2 = sum(rab**2)
    r = sqrt(r2)
 
    rco = dampening_fact*(self%rcov(at(a)) + self%rcov(at(b)))
 
    result = 1.0_wp/exp_count(k1, r, rco)
 
end subroutine

pure elemental function exp_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 1.0_wp/(1.0_wp + exp(-k*(r0/r - 1.0_wp)))
end function exp_count

end module  