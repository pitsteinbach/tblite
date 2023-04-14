module xtbml_base

    use mctc_env, only : wp
    use tblite_wavefunction, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_container, only : container_cache
    use tblite_results, only: results_type
    use mctc_io_convert, only : autoev
    
    use xtbml_class, only: xtbml_type

    type, public, extends(xtbml_type) :: xtbml_base_type
        
        
contains
    procedure :: get_xtbml
    procedure :: pack_res
end type xtbml_base_type
contains
    subroutine get_xtbml(self,mol,wfn,integrals,erep,calc,ccache,dcache,prlevel,res)
        use xtbml_functions
        implicit none
        class(xtbml_base_type),intent(inout) :: self
        type(structure_type), intent(in) :: mol
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(wavefunction_type), intent(in) :: wfn
        type(container_cache),intent(inout) :: ccache,dcache
        type(container_cache) :: dcache2
        type(results_type),intent(inout) :: res
        real(wp), INTENT(IN) ::  erep(mol%nat)
        integer, intent(in) :: prlevel
        real(wp) :: e_gfn2_tot
        integer :: ml_out
        
        self%n_features = 34
        self%a = 1.0_wp
        allocate(self%feature_labels(self%n_features))
        self%feature_labels = [ character(len=20) :: "CN","delta_CN",&
        &"q_s","q_p","q_d",&
        &"dipm_s","dipm_p","dipm_d",&
        &"qm_s","qm_p","qm_d",&
        &"p_A","delta_p_A","dipm_A","delta_dipm_A","qm_A","delta_qm",&
        &"delta_dipm_e","delta_qm_e","delta_dipm_Z","delta_qm_Z",&
        &"response","gap","chem.pot","HOAO","LUAO","E_repulsion","E_EHT",&
        &"E_disp_2","E_disp_3","E_ies_ixc","E_aes","E_axc","E_tot"]
        !get individual coulombic energy contributions in an atomwise vector
        call self%get_geometry_density_based(mol,wfn,integrals,calc)
        call self%get_energy_based(mol,wfn,calc,integrals,ccache,dcache,erep,e_gfn2_tot)

        call atomic_frontier_orbitals(mol%nat,calc%bas%nao,wfn%focc(:,1),wfn%emo(:,1)*autoev,calc%bas%ao2at,wfn%coeff(:,:,1),&
        integrals%overlap(:,:),self%response,self%egap,self%chempot,self%ehoao,self%eluao)

        
        call self%pack_res(mol%nat,calc%bas%nsh,calc%bas%nsh_at,e_gfn2_tot,res)
        allocate(res%w_xtbml(mol%nat),source=0.0_wp)
        res%w_xtbml = self%w_tot
        if (prlevel > 1) then
            ml_out = 42
            open(file='ml_feature_tblite.csv', newunit=ml_out)
            call self%print_out(ml_out,mol%nat,mol%num,mol%id,res)
        endif

    end subroutine get_xtbml

    subroutine pack_res(self,nat,nsh_tot,at2nsh,e_tot,res)
        use xtbml_functions , only: pack_shellwise
        implicit none
        integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
        real(wp), intent(in) :: e_tot
        type(results_type),intent(inout) :: res
        class(xtbml_base_type), intent(inout) :: self
        allocate(res%xtbml_labels(self%n_features),source=self%feature_labels)
        allocate(res%ml_features(nat,self%n_features),source=0.0_wp)
        res%ml_features(:,1) = self%cn_atom(:)
        res%ml_features(:,2) = self%delta_cn(:)
        call pack_shellwise(self%mulliken_shell,res,3,at2nsh,nat)
        call pack_shellwise(self%dipm_shell,res,6,at2nsh,nat)
        call pack_shellwise(self%qm_shell,res,9,at2nsh,nat)
        res%ml_features(:,12) = self%partial_charge_atom(:)
        res%ml_features(:,13) = self%delta_partial_charge(:)
        res%ml_features(:,14) = self%dipm_atom(:)
        res%ml_features(:,15) = self%delta_dipm(:)
        res%ml_features(:,16) = self%qm_atom(:)
        res%ml_features(:,17) = self%delta_qm(:)
        res%ml_features(:,18) = self%delta_dipm_e(:)
        res%ml_features(:,19) = self%delta_qm_e(:)
        res%ml_features(:,20) = self%delta_dipm_Z(:)
        res%ml_features(:,21) = self%delta_qm_Z(:)
        res%ml_features(:,22) = self%response(:)
        res%ml_features(:,23) = self%egap(:)
        res%ml_features(:,24) = self%chempot(:)
        res%ml_features(:,25) = self%ehoao(:)
        res%ml_features(:,26) = self%eluao(:)
        res%ml_features(:,27) = self%e_rep_atom(:)
        res%ml_features(:,28) = self%e_EHT(:)
        res%ml_features(:,29) = self%e_disp_2(:)
        res%ml_features(:,30) = self%e_disp_3(:)
        res%ml_features(:,31) = self%e_ies_ixc(:)
        res%ml_features(:,32) = self%e_aes(:)
        res%ml_features(:,33) = self%e_axc(:)
        res%ml_features(:,34) = e_tot
        
    end subroutine

end module xtbml_base
