module tblite_xtbml_energy_features
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
  use tblite_scf_iterator, only : get_electronic_energy, reduce
  use tblite_disp_d3, only : d3_dispersion, new_d3_dispersion
  use tblite_classical_halogen, only : halogen_correction
  use tblite_disp_d4, only : d4_dispersion, new_d4_dispersion
   implicit none
   private
  character(len=*), parameter :: label = "energy-based features"
   type, public, extends(xtbml_feature_type) :: xtbml_energy_features_type
      
    contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure :: setup
    end type

contains

subroutine setup(self)
  class(xtbml_energy_features_type) :: self
  self%label = label
  if (allocated(self%dict)) deallocate(self%dict)
  allocate(self%dict)
  !if (allocated(self%dict_ext)) deallocate(self%dict_ext)
  !allocate(self%dict_ext)
end subroutine

subroutine compute_features(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx)
  class(xtbml_energy_features_type), intent(inout) :: self
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
  type(d3_dispersion), allocatable :: d3
  type(d4_dispersion), allocatable :: d4
  class(container_type), allocatable :: cont
  real(wp), allocatable :: tmp_energy(:), e_ao(:), e_disp_tot(:), e_disp_ATM(:), tot_energy(:)
  integer :: i
  self%label = label
  
  allocate(e_ao(calc%bas%nao), source=0.0_wp)
  allocate(tmp_energy(mol%nat), tot_energy(mol%nat),source=0.0_wp)
  call get_electronic_energy(integrals%hamiltonian, wfn%density, e_ao)
  call reduce(tmp_energy, e_ao, calc%bas%ao2at)
  deallocate(e_ao)
  call self%dict%add_entry("E_EHT", tmp_energy)
  tot_energy = tmp_energy
  tmp_energy = 0.0_wp
  if (allocated(calc%repulsion)) then
    associate(cont => calc%repulsion)
      cache = cache_list(1)
      call cont%update(mol, cache)
      call cont%get_engrad(mol, cache, tmp_energy)
      call self%dict%add_entry("E_rep", tmp_energy)
      deallocate(cache)
    end associate
  end if
  tot_energy = tot_energy + tmp_energy
  tmp_energy = 0.0_wp
  if (allocated(calc%coulomb)) then
    cache = cache_list(2)
    associate(cont => calc%coulomb)
    call cont%update(mol, cache)
    if (allocated(cont%es2)) then
        call cont%es2%update(mol, cache)
        call cont%es2%get_energy(mol, cache, wfn, tmp_energy)
    end if
    if (allocated(cont%es3)) then
        call cont%es3%update(mol, cache)
        call cont%es3%get_energy(mol, cache, wfn, tmp_energy)
    end if
    tot_energy = tot_energy + tmp_energy
    call self%dict%add_entry("E_ies_ixc", tmp_energy)
    if (allocated(cont%aes2)) then
        tmp_energy = 0.0_wp
        call cont%aes2%get_AXC(mol, wfn, tmp_energy)
        call self%dict%add_entry("E_AXC", tmp_energy)
        tot_energy = tot_energy + tmp_energy
        tmp_energy = 0.0_wp 
        call cont%aes2%get_energy_aes_xtb(mol, cache, wfn, tmp_energy)
        call self%dict%add_entry("E_AES", tmp_energy)
        tot_energy = tot_energy + tmp_energy
    end if
    deallocate(cache)
    end associate
  end if
  tmp_energy = 0.0_wp
  if (allocated(calc%halogen)) then 
    cache = cache_list(3)
    associate(cont => calc%halogen)
      call cont%update(mol, cache)
      call cont%get_engrad(mol, cache, tmp_energy)
      call self%dict%add_entry("E_HX", tmp_energy)
    end associate
    deallocate(cache)
  end if
  tot_energy = tot_energy + tmp_energy
  tmp_energy = 0.0_wp
  if (allocated(calc%dispersion)) then
    cache = cache_list(4)
    associate(cont => calc%dispersion)
    select type(cont)
    type is (d3_dispersion)
        allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
        call cont%update(mol, cache)
        call cont%get_engrad(mol, cache, e_disp_tot)
        
        allocate(d3)
        call new_d3_dispersion(d3, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, a2=cont%param%a2, s9=cont%param%s9)
        call d3%update(mol, cache)
        call d3%get_engrad(mol, cache, e_disp_ATM)
        call self%dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
        call self%dict%add_entry("E_disp3", e_disp_ATM)
        tmp_energy = e_disp_tot
        deallocate(d3, e_disp_ATM, e_disp_tot)
    type is (d4_dispersion)
        allocate(e_disp_tot(mol%nat), e_disp_ATM(mol%nat), source=0.0_wp)
        call cont%update(mol, cache)
        call cont%get_engrad(mol, cache, e_disp_tot)
        call cont%get_energy(mol, cache, wfn, e_disp_tot)
            
        allocate(d4)
        call new_d4_dispersion(d4, mol, s6=0.0_wp, s8=0.0_wp, a1=cont%param%a1, a2=cont%param%a2, s9=cont%param%s9)
        call d4%update(mol, cache)
        call d4%get_engrad(mol, cache, e_disp_ATM)
        call self%dict%add_entry("E_disp2", e_disp_tot-e_disp_ATM)
        call self%dict%add_entry("E_disp3", e_disp_ATM)
        tmp_energy = e_disp_tot
        deallocate(d4, e_disp_tot, e_disp_ATM)
    end select
    end associate
    deallocate(cache)
  end if
  tot_energy = tot_energy + tmp_energy
  tmp_energy = 0.0_wp
  if (allocated(calc%interactions)) then
    cache = cache_list(5)
    associate(cont => calc%dispersion)
        call cont%update(mol, cache)
        call cont%get_engrad(mol, cache, tmp_energy)
        call cont%get_energy(mol, cache, wfn, tmp_energy)
        call self%dict%add_entry(cont%info(0, ""), tmp_energy)
    end associate
    deallocate(cache)
  end if
  tot_energy = tot_energy + tmp_energy
  call self%dict%add_entry("E_tot", tot_energy)
  call self%dict%add_entry("w_tot", tot_energy/sum(tot_energy))
  deallocate(tot_energy)

end subroutine

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, prlevel, ctx, convolution)
  use tblite_output_format, only : format_string
  class(xtbml_energy_features_type), intent(inout) :: self
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