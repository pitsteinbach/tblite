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
    integer,parameter :: n_delta = 12
    type, public, extends(xtbml_type) :: xtbml_base_type
        
     

contains
    procedure :: get_xtbml
    procedure :: pack_res
end type xtbml_base_type
contains
subroutine get_xtbml(self, mol, wfn, integrals, calc, ccache, dcache,rcache, prlevel, ctx, res)
   use tblite_xtbml_functions
   use tblite_timer, only : timer_type
   implicit none

   class(xtbml_base_type), intent(inout) :: self
   type(structure_type), intent(in) :: mol
   type(integral_type) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   type(container_cache), intent(inout) :: ccache, dcache, rcache
   type(results_type), intent(inout) :: res
   type(context_type), intent(inout) :: ctx

   integer, intent(in) :: prlevel
   
   real(wp) :: e_gfn2_tot
   integer :: ml_out
   logical :: print_afo
   character(len=30), allocatable :: tmp_labels(:)
   
   self%n_features = 40 + ((size(self%a) - 1)*n_delta)
   allocate (self%feature_labels(self%n_features))
   tmp_labels = [character(len=30) :: "CN",&
      &"p_s", "p_p", "p_d",&
      &"dipm_s", "dipm_p", "dipm_d",&
      &"qm_s", "qm_p", "qm_d",&
      &"q_A", "dipm_A", "qm_A",&
      &"response", "gap", "chem_pot", "HOAO_a", "LUAO_a", "HOAO_b", "LUAO_b",&
      &"E_repulsion", "E_EHT",&
      &"E_disp_2", "E_disp_3", "E_ies_ixc", "E_aes", "E_axc", "E_tot"]

   allocate (self%delta_labels(n_delta))
   self%delta_labels = [character(len=30) :: "delta_CN",&
      &"delta_q_A", "delta_dipm_A", "delta_qm_A",&
      &"delta_dipm_e", "delta_qm_e", "delta_dipm_Z", "delta_qm_Z",&
      &"delta_gap", "delta_chem_pot", "delta_HOAO", "delta_LUAO"]
   allocate(self%geometry_features)
   call self%geometry_features%compute_features(mol, wfn , integrals, calc, dcache, prlevel, ctx)
   call self%geometry_features%compute_extended(mol, wfn , integrals, calc, dcache, prlevel, ctx)

   !get individual coulombic energy contributions in an atomwise vector
   call self%get_geometry_density_based(mol, wfn, integrals, calc)
   call self%get_energy_based(mol, wfn, calc, integrals, ccache, dcache, erep, e_gfn2_tot)
   print_afo = .false.
   if (prlevel > 1) then
      print_afo = .true.
   end if

   call self%get_frontier(mol,calc%bas,wfn,integrals%overlap(:,:),print_afo,ctx)

   call self%compute_extended(mol, wfn, calc)
   call self%get_extended_frontier(mol, wfn)
   call self%pack_res(mol%nat, calc%bas%nsh, calc%bas%nsh_at, e_gfn2_tot, tmp_labels, res)

   if (prlevel > 1) then
      ml_out = 42
      open (file='ml_feature_tblite.csv', newunit=ml_out)
      call self%print_out(ml_out, mol%nat, mol%num, mol%id, res)
   end if
   deallocate (self%delta_labels)
   if (debug) then
      call self%print_timer(ctx)
   end if
end subroutine get_xtbml

    subroutine pack_res(self,nat,nsh_tot,at2nsh,e_tot,labels,res)
        use xtbml_functions , only: pack_shellwise
        use tblite_output_format, only: format_string
        implicit none
        integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
        real(wp), intent(in) :: e_tot
        type(results_type),intent(inout) :: res
        class(xtbml_base_type), intent(inout) :: self
        integer :: k, n_other, i , offset
        character(len=30),intent(in) :: labels(:)

        allocate(res%ml_features(nat,self%n_features),source=0.0_wp)
        res%ml_features(:,1) = self%cn_atom(:)
        call pack_shellwise(self%mulliken_shell,res,2,at2nsh,nat)
        call pack_shellwise(self%dipm_shell,res,5,at2nsh,nat)
        call pack_shellwise(self%qm_shell,res,8,at2nsh,nat)
        res%ml_features(:,11) = self%partial_charge_atom(:)
        res%ml_features(:,12) = self%dipm_atom(:)
        res%ml_features(:,13) = self%qm_atom(:)
        res%ml_features(:,14) = self%response(:)
        res%ml_features(:,15) = self%egap(:)
        res%ml_features(:,16) = self%chempot(:)
        res%ml_features(:,17) = self%ehoao_a(:)
        res%ml_features(:,18) = self%eluao_a(:)
        res%ml_features(:,19) = self%ehoao_b(:)
        res%ml_features(:,20) = self%eluao_b(:)
        res%ml_features(:,21) = self%e_rep_atom(:)
        res%ml_features(:,22) = self%e_EHT(:)
        res%ml_features(:,23) = self%e_disp_2(:)
        res%ml_features(:,24) = self%e_disp_3(:)
        res%ml_features(:,25) = self%e_ies_ixc(:)
        res%ml_features(:,26) = self%e_aes(:)
        res%ml_features(:,27) = self%e_axc(:)
        res%ml_features(:,28) = e_tot
        n_other = 28
        do i =1, n_other
            self%feature_labels(i) = trim(labels(i))
        enddo
            
        do k =1 , size(self%a)
            offset = n_other + (k-1)*n_delta
            res%ml_features(:,offset+1) = self%delta_cn(:,k)
            res%ml_features(:,offset+2) = self%delta_partial_charge(:,k)
            res%ml_features(:,offset+3) = self%delta_dipm(:,k)
            res%ml_features(:,offset+4) = self%delta_qm(:,k)
            res%ml_features(:,offset+5) = self%delta_dipm_e(:,k)
            res%ml_features(:,offset+6) = self%delta_qm_e(:,k)
            res%ml_features(:,offset+7) = self%delta_dipm_Z(:,k)
            res%ml_features(:,offset+8) = self%delta_qm_Z(:,k)
            res%ml_features(:,offset+9) = self%delta_egap(:,k)
            res%ml_features(:,offset+10) = self%delta_chempot(:,k)
            res%ml_features(:,offset+11) = self%delta_ehoao(:,k)
            res%ml_features(:,offset+12) = self%delta_eluao(:,k)
            do i = 1, n_delta
                self%feature_labels(offset+i) = trim(self%delta_labels(i)) //'_' //adjustl(format_string(self%a(k),'(f12.2)'))
            enddo
        enddo

        res%n_features = self%n_features
        call move_alloc(self%feature_labels,res%xtbml_labels)
        call move_alloc(self%w_tot,res%w_xtbml)
        

    end subroutine

end module xtbml_base
