module tblite_xtbml_potenetial_features
  use mctc_env, only : wp
  use mctc_io_convert, only : autoev
  use tblite_xtbml_feature_type, only : xtbml_feature_type
  use tblite_wavefunction_type, only : wavefunction_type
  use mctc_io, only : structure_type
  use tblite_integral_type, only : integral_type
  use tblite_basis_type, only : basis_type 
  use tblite_container, only : container_cache
  use tblite_context , only : context_type
  use tblite_double_dictionary, only : double_dictionary_type   
  use tblite_xtbml_convolution, only : xtbml_convolution_type
  use tblite_xtbml_atomic_frontier, only : atomic_frontier_orbitals
  use tblite_container, only : container_list, container_type
  use tblite_xtb_calculator, only : xtb_calculator
  use tblite_repulsion, only : tb_repulsion
  use tblite_xtb_coulomb, only : tb_coulomb
  use tblite_classical_halogen, only : halogen_correction
  use tblite_scf_potential, only : potential_type, new_potential
  use tblite_scf_iterator, only : reduce
  implicit none
  private
  character(len=*), parameter :: label = "potential-based features"
   type, public, extends(xtbml_feature_type) :: xtbml_potential_features_type
      
    contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure :: setup
    end type

contains

subroutine setup(self)
  class(xtbml_potential_features_type) :: self
  self%label = label
  if (allocated(self%dict)) deallocate(self%dict)
  allocate(self%dict)
  !if (allocated(self%dict_ext)) deallocate(self%dict_ext)
  !allocate(self%dict_ext)
end subroutine

subroutine compute_features(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx)
  class(xtbml_potential_features_type), intent(inout) :: self
  !> Molecular structure data
  type(structure_type), intent(in) :: mol
  !> Wavefunction strcuture data
  type(wavefunction_type), intent(in) :: wfn
  !> Integral container
  type(integral_type) :: integrals
  !> Single-point calculator
  type(xtb_calculator), intent(in) :: calc
  !> List of containers 
  type(container_cache), intent(inout) :: cache_list(:)
  type(container_cache), allocatable :: cache
  !> Context type
  type(context_type),intent(inout) :: ctx
  !> Print Level
  integer, intent(in) :: prlevel
  !> Container for cached variables
  class(container_type), allocatable :: cont
  !> potential type
  type(potential_type) :: pot
  real(kind=wp), allocatable :: vat(:)
  integer :: i
  self%label = label
  
  allocate(vat(mol%nat), source=0.0_wp)
  call new_potential(pot, mol, calc%bas, wfn%nspin)
  if (allocated(calc%coulomb)) then
    cache = cache_list(2)
    associate(cont => calc%coulomb)
    call cont%update(mol, cache)

    if (allocated(cont%es2)) then
        call cont%es2%update(mol, cache)
        call cont%es2%get_potential(mol, cache, wfn, pot)
    end if
    call reduce(vat, pot%vsh(:, 1), calc%bas%sh2at)
    call self%dict%add_entry("V_IES", vat)
    pot%vsh = 0.0_wp
    if (allocated(cont%es3)) then
        call cont%es3%update(mol, cache)
        call cont%es3%get_potential(mol, cache, wfn, pot)
    end if
    vat = 0.0_wp
    if (cont%es3%shell_resolved) then
      call reduce(vat, pot%vsh(:, 1), calc%bas%sh2at)
      call self%dict%add_entry("V_IXC", vat)
    else
      call self%dict%add_entry("V_IXC", pot%vat)
    end if
    deallocate(cache)
    end associate
  end if

end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx, convolution)
  use tblite_output_format, only : format_string
  class(xtbml_potential_features_type), intent(inout) :: self
  !> Molecular structure data
  type(structure_type), intent(in) :: mol
  !> Wavefunction strcuture data
  type(wavefunction_type), intent(in) :: wfn
  !> Integral container
  type(integral_type) :: integrals
  !> Single-point calculator
  type(xtb_calculator), intent(in) :: calc
  type(container_cache), intent(inout) :: cache_list(:)
  !> Context type
  type(context_type),intent(inout) :: ctx
  !> Print Level
  integer, intent(in) :: prlevel
  !> Convolution container
  type(xtbml_convolution_type) :: convolution
end subroutine

end module