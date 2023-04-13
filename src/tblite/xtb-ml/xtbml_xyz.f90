module xtbml_xyz

    use mctc_env, only : wp
    use tblite_wavefunction, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_container, only : container_cache
    use tblite_scf_iterator, only : get_electronic_energy,reduce
    use tblite_results, only: results_type
    use mctc_io_convert, only : autoev
    
    use xtbml_class, only: xtbml_type


type, public, extends(xtbml_type) :: xtbml_xyz_type
        
contains
    procedure :: get_xtbml
    procedure :: pack_res
end type xtbml_xyz_type
contains
    subroutine get_xtbml(self,mol,wfn,integrals,erep,calc,ccache,dcache,prlevel,res)
        use xtbml_functions
        implicit none
        class(xtbml_xyz_type),intent(inout) :: self
        type(structure_type), intent(in) :: mol
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(wavefunction_type), intent(in) :: wfn
        type(container_cache),intent(inout) :: ccache,dcache
        type(results_type),intent(inout) :: res
        real(wp), INTENT(IN) ::  erep(mol%nat)
        integer, intent(in) :: prlevel
        real(wp) :: e_gfn2_tot
        integer :: ml_out
        
        self%n_features = 97
        self%a = 1.0_wp
        allocate(self%feature_labels(self%n_features))
        self%feature_labels = [ character(len=20) :: "CN","delta_CN",&
        &"q_s","q_p","q_d",&
        &"dipm_s","dipm_p","dipm_d",&
        &"dipm_s_x","dipm_s_y","dipm_s_z",&
        &"dipm_p_x","dipm_p_y","dipm_p_z",&
        &"dipm_d_x","dipm_d_y","dipm_d_z",&
        &"qm_s","qm_p","qm_d",&
        &"qm_s_xx","qm_s_xy","qm_s_yy","qm_s_xz","qm_s_yz","qm_s_zz",&
        &"qm_p_xx","qm_p_xy","qm_p_yy","qm_p_xz","qm_p_yz","qm_p_zz",&
        &"qm_d_xx","qm_d_xy","qm_d_yy","qm_d_xz","qm_d_yz","qm_d_zz",&
        &"p_A","delta_p_A","dipm_A",&
        &"dipm_A_x","dipm_A_y","dipm_A_z",&
        &"delta_dipm_A",&
        &"delta_dipm_A_x","delta_dipm_A_y","delta_dipm_A_z",&
        &"qm_A",&
        &"qm_A_xx","qm_A_xy","qm_A_yy","qm_A_xz","qm_A_yz","qm_A_zz",&
        &"delta_qm",&
        &"delta_qm_A_xx","delta_qm_A_xy","delta_qm_A_yy","delta_qm_A_xz","delta_qm_A_yz","delta_qm_A_zz",&
        &"delta_dipm_e",&
        &"delta_dipm_e_x","delta_dipm_e_y","delta_dipm_e_z",&
        &"delta_qm_e",&
        &"delta_qm_e_xx","delta_qm_e_xy","delta_qm_e_yy","delta_qm_e_xz","delta_qm_e_yz","delta_qm_e_zz",&
        &"delta_dipm_Z",&
        &"delta_dipm_Z_x","delta_dipm_Z_y","delta_dipm_Z_z",&
        &"delta_qm_Z",&
        &"delta_qm_Z_xx","delta_qm_Z_xy","delta_qm_Z_yy","delta_qm_Z_xz","delta_qm_Z_yz","delta_qm_Z_zz",&
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
        use xtbml_functions, only : pack_mult_xyz, pack_mult_xyz_shell,pack_shellwise
        implicit none
        integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
        real(wp), intent(in) :: e_tot
        type(results_type),intent(inout) :: res
        class(xtbml_xyz_type), intent(inout) :: self
        integer :: i, nsh
        res%n_features = self%n_features
        allocate(res%xtbml_labels(res%n_features))
        res%xtbml_labels = self%feature_labels
        allocate(res%ml_features(nat,self%n_features),source=0.0_wp)
        res%ml_features(:,1) = self%cn_atom(:)
        res%ml_features(:,2) = self%delta_cn(:)
        call pack_shellwise(self%mulliken_shell,res,3,at2nsh,nat) !3-5
        call pack_shellwise(self%dipm_shell,res,6,at2nsh,nat) ! 6-8
        call pack_mult_xyz_shell(self%dipm_shell_xyz,res,9,nat,at2nsh) !packs xyz for s to d shell 9-17
        call pack_shellwise(self%qm_shell,res,18,at2nsh,nat) ! 18-20
        call pack_mult_xyz_shell(self%qm_shell_xyz,res,21,nat,at2nsh) !21-38
        write(*,*) self%qm_shell_xyz(:,1:3)
        res%ml_features(:,39) = self%partial_charge_atom(:)
        res%ml_features(:,40) = self%delta_partial_charge(:)
        res%ml_features(:,41) = self%dipm_atom(:)
        call pack_mult_xyz(self%dipm_atom_xyz,res,42,nat) !42-44
        res%ml_features(:,45) = self%delta_dipm(:)
        call pack_mult_xyz(self%delta_dipm_xyz,res,46,nat) !46-48
        res%ml_features(:,49) = self%qm_atom(:)
        call pack_mult_xyz(self%qm_atom_xyz,res,50,nat) !50-55
        res%ml_features(:,56) = self%delta_qm(:)
        call pack_mult_xyz(self%delta_qm_xyz,res,57,nat) !57-62
        res%ml_features(:,63) = self%delta_dipm_e(:)
        call pack_mult_xyz(self%delta_dipm_e_xyz,res,64,nat) !64-66
        res%ml_features(:,67) = self%delta_qm_e(:)
        call pack_mult_xyz(self%delta_qm_e_xyz,res,68,nat) !66-71
        res%ml_features(:,74) = self%delta_dipm_Z(:)
        call pack_mult_xyz(self%delta_dipm_Z_xyz,res,75,nat) !73-75
        res%ml_features(:,78) = self%delta_qm_Z(:)
        call pack_mult_xyz(self%delta_qm_e_xyz,res,79,nat) !77-82
        res%ml_features(:,85) = self%response(:)
        res%ml_features(:,86) = self%egap(:)
        res%ml_features(:,87) = self%chempot(:)
        res%ml_features(:,88) = self%ehoao(:)
        res%ml_features(:,89) = self%eluao(:)
        res%ml_features(:,90) = self%e_rep_atom(:)
        res%ml_features(:,91) = self%e_EHT(:)
        res%ml_features(:,92) = self%e_disp_2(:)
        res%ml_features(:,93) = self%e_disp_3(:)
        res%ml_features(:,94) = self%e_ies_ixc(:)
        res%ml_features(:,95) = self%e_aes(:)
        res%ml_features(:,96) = self%e_axc(:)
        res%ml_features(:,97) = e_tot
        
    end subroutine


end module xtbml_xyz
