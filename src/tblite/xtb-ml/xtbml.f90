module xtbml_feature_calc

    use mctc_env, only : wp
    use tblite_wavefunction, only : wavefunction_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_container, only : container_cache
    use tblite_scf_iterator, only : get_electronic_energy,reduce
    use tblite_data_covrad, only : get_covalent_rad
    use tblite_wavefunction_mulliken, only: get_mulliken_shell_multipoles
    use tblite_results, only: results_type
    use mctc_io_convert, only : autoev
    use tblite_ncoord_exp, only:  new_exp_ncoord,exp_ncoord_type
    !use tblite_ncoord_type, only: ncoord_type
    !use xtb_mctc_convert, only : evtoau
    !use xtb_type_data
    use xtb_type_ml

    real(wp),allocatable :: rcov(:)
    real(wp),parameter :: k1 = 16.0_wp
    real(wp), PARAMETER :: dampening_fact = 1.0_wp
    integer,parameter :: nfeatures = 34


contains
    subroutine do_ML_print(mol,wfn,integrals,erep,calc,ccache,dcache,res)
        implicit none
        type(ml_features) :: ml
        type(structure_type), intent(in) :: mol
        type(integral_type) :: integrals
        !> Single-point calculator
        type(xtb_calculator), intent(in) :: calc
        type(wavefunction_type), intent(inout) :: wfn
        type(container_cache) :: ccache,dcache
        type(container_cache) :: dcache2
        type(exp_ncoord_type) :: ncoord_exp
        type(results_type) :: res
        real(wp), INTENT(IN) ::  erep(mol%nat)
        real(wp) :: e_gfn2_tot
        real(wp), allocatable :: e_ao(:),e_disp(:),e_1e(:)
        real(wp) :: dipm_shellwise(3,calc%bas%nsh,wfn%nspin), qp_shellwise(6,calc%bas%nsh,wfn%nspin), dipm_atom(3,mol%nat), qp_atom(6,mol%nat)
        real(wp) :: delta_dipm(3,mol%nat),delta_qp(6,mol%nat),z(mol%nat)
        reaL(wp) :: delta_dipm_only_p(3,mol%nat), delta_qp_only_p(6,mol%nat),delta_dipm_only_Z (3,mol%nat),delta_qp_only_Z(6,mol%nat), mull_charge_atomic(mol%nat)
        integer :: ml_out
        !allocate ml type
        call ml%allocate(mol%nat,size(wfn%qsh))
        
        !get individual coulombic energy contributions in an atomwise vector
        !call calc%coulomb%aes2%update(mol,ccache)
        call calc%coulomb%aes2%get_energy(mol,ccache,wfn,ml%e_aes)
        !write(*,*) ml%e_aes
        !call calc%coulomb%es2%update(mol,ccache)
        call calc%coulomb%es2%get_energy(mol,ccache,wfn,ml%e_ies_ixc)
        !call calc%coulomb%es3%update(mol,ccache)
        call calc%coulomb%es3%get_energy(mol,ccache,wfn,ml%e_ies_ixc)

        call calc%dispersion_2body%update(mol,dcache2)
        call calc%dispersion_2body%get_engrad(mol,dcache2,ml%e_disp_3)
        
        allocate (e_disp(mol%nat),source=0.0_wp)
        call calc%dispersion%update(mol,dcache)
        call calc%dispersion%get_energy(mol,dcache,wfn,e_disp)

        ml%e_disp_2(:) = e_disp(:) - ml%e_disp_3(:)

        ml%e_rep_atom = erep
        
        !Compute E_EHT
        allocate(e_ao(calc%bas%nao),source=0.0_wp)
        call get_electronic_energy(integrals%hamiltonian,wfn%density,e_ao)
        call reduce(ml%e_EHT,e_ao,calc%bas%ao2at)
        !Compute partition weights based on total energy expression
        call get_total_xtb_weights(mol%nat,ml%e_EHT,ml%e_rep_atom,ml%e_disp_2,ml%e_disp_3,&
                                ml%e_ies_ixc,ml%e_aes,ml%e_axc,ml%w_tot,e_gfn2_tot)
        
        !compute delta CN
        allocate(rcov(mol%nid),source=0.0_wp)
        rcov(:) = get_covalent_rad(mol%num)
        !write(*,*) rcov
        call new_exp_ncoord(ncoord_exp,mol)
        call ncoord_exp%get_cn(mol,ml%cn_atom)
        !write(*,*) mol%xyz
        call get_delta_cn(mol%nat,ml%cn_atom,mol%id,mol%xyz,ml%delta_cn)

        !shellwise mulliken charges
        
        call mulliken_shellwise(calc%bas%nao,calc%bas%nsh,calc%bas%ao2sh,wfn%density(:,:,wfn%nspin),integrals%overlap,ml%shell_mulliken)
        !call sum_up_spin(wfn%qsh,ml%shell_mulliken)
        !write(*,*) ml%shell_mulliken(:)
        call sum_up_spin(wfn%qat,ml%atomic_partialcharge)
        !write(*,*) wfn%qat(:,:)
        
        ml%atomic_partialcharge = wfn%qat(:,1)
        !delta partial charge
        call get_delta_partial(mol%nat,ml%atomic_partialcharge,mol%id,mol%xyz,&
        ml%cn_atom,ml%delta_partial_charge)

        !multipole moments shellwise und then atomwise
        call get_mulliken_shell_multipoles(calc%bas,integrals%dipole, wfn%density, &
        & dipm_shellwise)
        call get_mulliken_shell_multipoles(calc%bas,integrals%quadrupole, wfn%density, &
        & qp_shellwise)
        
        !delta multipole moments
        call get_delta_mm(mol%nat,ml%atomic_partialcharge,wfn%dpat,wfn%qpat,mol%id,mol%xyz,&
        ml%cn_atom,delta_dipm,delta_qp)

        call comp_norm(calc%bas%nsh,dipm_shellwise,qp_shellwise,ml%shell_dipm,ml%shell_qm)
        call comp_norm(mol%nat,wfn%dpat,wfn%qpat,ml%dipm_atom,ml%qm_atom)
        call comp_norm(mol%nat,delta_dipm,delta_qp,ml%delta_dipm,ml%delta_qm)
        
        call mol_set_nuclear_charge(mol%nat,mol%num,mol%id,z)
        call get_delta_mm_Z(mol%nat,z,wfn%dpat,wfn%qpat,mol%id,mol%xyz,ml%cn_atom,delta_dipm_only_Z,delta_qp_only_Z)
        
        !extended CAMMs only mulliken charges
        call sum_up_mulliken(mol%nat,calc%bas%nsh,calc%bas%ao2at,calc%bas%sh2at,ml%shell_mulliken,mull_charge_atomic)
        call get_delta_mm_p(mol%nat,mull_charge_atomic,wfn%dpat,wfn%qpat,mol%id,mol%xyz,ml%cn_atom,&
        delta_dipm_only_p,delta_qp_only_p)

        call comp_norm(mol%nat,delta_dipm_only_p,delta_qp_only_p,ml%delta_dipm_only_p,ml%delta_qp_only_p)
        call comp_norm(mol%nat,delta_dipm_only_Z,delta_qp_only_Z,ml%delta_dipm_only_Z,ml%delta_qp_only_Z)
        
        call atomic_frontier_orbitals(mol%nat,calc%bas%nao,wfn%focc(:,1),wfn%emo(:,1)*autoev,calc%bas%ao2at,wfn%coeff(:,:,1),&
        integrals%overlap(:,:),ml%response,ml%egap,ml%chempot,ml%ehoao,ml%eluao)

        ml_out = 42
        open(file='ml_feature_tblite.csv', newunit=ml_out)
        call pack_res(mol%nat,calc%bas%nsh,calc%bas%nsh_at,ml,e_gfn2_tot,res)
        allocate(res%w_xtbml(mol%nat),source=0.0_wp)
        res%w_xtbml = ml%w_tot
        call print_ML_output(ml_out,mol%nat,calc%bas%nao,calc%bas%nsh,mol%num,mol%id,calc%bas%sh2at,calc%bas%ao2at,res)
              
        !calc%bas%nsh_at

    end subroutine do_ML_print

    

    subroutine get_total_xtb_weights(nat,e_1e,e_rep,e_diff_atom,e_disp_3,e_ies_ixc,e_aes,e_axc,w_tot_xtb,e_tot)
        integer, INTENT(IN) :: nat
        real(wp), INTENT(IN) :: e_1e(nat),e_diff_atom(nat),e_disp_3(nat),e_ies_ixc(nat),e_aes(nat), e_axc(nat),e_rep(nat)
        real(wp), INTENT(OUT) :: w_tot_xtb(nat), e_tot
        real(wp) :: sum_energy, sum_atom(nat)
        integer :: i

        sum_atom = e_1e + e_rep + e_diff_atom + e_disp_3 + e_ies_ixc + e_aes + e_axc
        sum_energy = sum(sum_atom(:))

        w_tot_xtb(:) = sum_atom / sum_energy 

        e_tot = sum_energy

    end subroutine

    subroutine get_delta_cn(nat,cn,at,xyz,delta_cn)
        implicit none 
        integer, INTENT(IN) :: nat, at(nat)
        real(wp), INTENT(IN) :: cn(nat), xyz(3,nat)
        real(wp), INTENT(OUT) :: delta_cn(nat)
        integer :: i,j
        real(wp) :: result
        result = 0.0_wp
        delta_cn = 0.0_wp

        do i = 1, nat
           ! write(*,*) xyz(1:3,i)
            do j = 1, nat
                if (i == j) cycle 
                call inv_cn(nat,i,j,at,xyz,result)
                !write(*,*) result
                delta_cn(i) = delta_cn(i) + cn(j) / result
            enddo
        enddo
    end subroutine

    subroutine inv_cn(nat,a,b,at,xyz,result)
        implicit none
        integer, INTENT(IN) :: a, b,nat
        integer, INTENT(IN) :: at(:)
        real(wp), intent(in)  :: xyz(3,nat) 
        real(wp), INTENT(OUT) :: result
        real(wp) :: rab(3), r, rco, den, tmp, r2
     
        result = 0.0_wp
        
        rab = xyz(1:3,a) - xyz(1:3,b)
        r2 = sum( rab**2 )
        r = sqrt(r2)
        
        rco=dampening_fact*(rcov(at(a)) + rcov(at(b)))
        
        result = 1.0_wp / exp_count(k1,r,rco)
     
    end subroutine

    pure elemental function exp_count(k,r,r0) result(count)
        real(wp), intent(in) :: k
        real(wp), intent(in) :: r
        real(wp), intent(in) :: r0
        real(wp) :: count
        count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
    end function exp_count

    subroutine sum_up_spin(q_2,q_1)
        real(wp),intent(in):: q_2(:,:)
        real(wp),intent(out) :: q_1(:)
        integer :: i
        do i = 1, size(q_2,1)
            q_1(i) = sum(q_2(i,:))
        end do
    end subroutine sum_up_spin

    subroutine mulliken_shellwise(nao,nshell,ao2shell,p,s,charges_shell)
        implicit none
        integer, intent(in) :: nao,nshell,ao2shell(:)
        real(wp), intent(in) :: s(:, :)
        real(wp), intent(in) :: p(:, :)
        real(wp), intent(out) :: charges_shell(nshell)
        integer :: a, mu, nu
        
        charges_shell = 0.0_wp
        
        do mu = 1, nao
            do nu = 1, nao
                charges_shell(ao2shell(mu)) = charges_shell(ao2shell(mu)) + p(mu,nu) * s(nu,mu)
            enddo
        enddo

    end subroutine

    subroutine get_delta_partial(nat,atom_partial,at,xyz,cn,delta_partial)
        implicit none 
        integer, INTENT(IN) :: nat , at(nat)
        real(wp), INTENT(IN) :: atom_partial(nat), xyz(3,nat), cn(nat)
        real(wp), INTENT(OUT) :: delta_partial(nat)
        integer :: a, b
        real(wp) :: result

        delta_partial = 0.0_wp
        
        do a = 1, nat
            do b = 1, nat
                !if (a == b) cycle 
                call inv_cn(nat,a,b,at,xyz,result)
                delta_partial(a) = delta_partial(a) + atom_partial(b) / (result*(cn(b)+1))
            enddo
            
        enddo

    end subroutine

    subroutine sum_up_mm(nat,nshell,aoat2,ash,dipm_shell,qm_shell,dipm_at,qm_at)
        implicit none 
        integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
        real(wp), INTENT(IN) :: dipm_shell(:,:), qm_shell(:,:)
        real(wp), INTENT(OUT) :: dipm_at(3,nat), qm_at(6,nat)
        integer :: i

        dipm_at = 0.0_wp
        qm_at = 0.0_wp
        
        do i = 1, nshell
            dipm_at(:,ash(i)) = dipm_at(:,ash(i)) + dipm_shell(:,i)
            qm_at(:,ash(i)) = qm_at(:,ash(i)) + qm_shell(:,i)
        enddo

    end subroutine

    subroutine get_delta_mm(nat,q,dipm,qp,at,xyz,cn,delta_dipm,delta_qp)
        implicit none 
        integer, INTENT(IN) :: nat, at(nat)
        real(wp), INTENT(IN) :: dipm(3,nat), xyz(3,nat), qp(6,nat), q(nat), cn(nat)
        real(wp), INTENT(OUT) :: delta_dipm(3,nat), delta_qp(6,nat)
        integer :: a, b, i
        real(wp) :: result, r_ab(3), tii, qp_part(6,nat)

        !$acc enter data create(delta_dipm(:, :), delta_qp(:, :),qp_part(:,:))
        !$acc kernels default(present)
        delta_dipm = 0.0_wp
        delta_qp = 0.0_wp
        qp_part = 0.0_wp
        !$acc end kernels
        !$acc enter data copyin( nat, , q(:), dipm(:, :), at(:),&
        !$acc& qpint(:, :, :),xyz(:, :)),cn(:)

        !$acc parallel private(r_ab,result)

        !$acc loop gang vector collapse(2)

        do a = 1, nat
            do b = 1, nat
                !if (a == b) cycle 
                call inv_cn(nat,a,b,at,xyz,result)
                r_ab = xyz(:,a) - xyz(:,b)
                !$acc atomic
                delta_dipm(:,a) = delta_dipm(:,a) + (dipm(:,b) - r_ab(:) * q(b)) / (result*(cn(b)+1))
                !sorting of qp xx,xy,yy,xz,yz,zz
                !$acc atomic
                delta_qp(1,a) = delta_qp(1,a) + ( 1.5_wp*(-1*(r_ab(1)*dipm(1,b) + r_ab(1)*dipm(1,b)) + r_ab(1)*r_ab(1)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(2,a) = delta_qp(2,a) + ( 1.5_wp*(-1*(r_ab(1)*dipm(2,b) + r_ab(2)*dipm(1,b)) + r_ab(1)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(3,a) = delta_qp(3,a) + ( 1.5_wp*(-1*(r_ab(2)*dipm(2,b) + r_ab(2)*dipm(2,b)) + r_ab(2)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(4,a) = delta_qp(4,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(1,b) + r_ab(1)*dipm(3,b)) + r_ab(1)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(5,a) = delta_qp(5,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(2,b) + r_ab(2)*dipm(3,b)) + r_ab(2)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(6,a) = delta_qp(6,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(3,b) + r_ab(3)*dipm(3,b)) + r_ab(3)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                qp_part(:,a) = qp_part(:,a) + qp(:,b) / (result*(cn(b)+1))
            enddo
        enddo
        !$acc end parallel

        !$acc exit data copyout(delta_dipm(:, :), delta_qp(:, :), qp_part(:,:))

        call remove_trac_qp(nat,delta_qp,qp_part)

    end subroutine

    subroutine comp_norm(ndim,dipm,qm,dipm_norm,qm_norm)
        implicit none
        integer, INTENT(IN) :: ndim
        real(wp), INTENT(IN) :: dipm(3,ndim), qm(6,ndim)
        real(wp), INTENT(OUT) :: dipm_norm(ndim),qm_norm(ndim)
        real(wp) :: r(ndim), r2(ndim)
        INTEGER :: i, j

        do i = 1, ndim
            r2(i) = dipm(1,i)**2+dipm(2,i)**2+dipm(3,i)**2
        enddo
        r = sqrt(r2)
        dipm_norm = r

        do i = 1, ndim
            r2(i) = qm(1,i)**2 + 2*qm(2,i)**2 + qm(3,i)**2 + 2*qm(4,i)**2 + 2*qm(5,i)**2 + qm(6,i)**2
        enddo
        r = sqrt(r2)
        qm_norm = r

    end subroutine

    subroutine get_delta_mm_Z(nat,q,dipm,qp,at,xyz,cn,delta_dipm,delta_qp)! all effects due to the electrons are set to 0, only the distribution of positive charges is left
        implicit none 
        integer, INTENT(IN) :: nat, at(nat)
        real(wp), INTENT(IN) :: dipm(3,nat), xyz(3,nat), qp(6,nat), q(nat), cn(nat)
        real(wp), INTENT(OUT) :: delta_dipm(3,nat), delta_qp(6,nat)
        integer :: a, b, i
        real(wp) :: result, r_ab(3), tii, qp_part(6,nat)

        !$acc enter data create(delta_dipm(:, :), delta_qp(:, :),qp_part(:,:))
        !$acc kernels default(present)
        delta_dipm = 0.0_wp
        delta_qp = 0.0_wp
        qp_part = 0.0_wp
        !$acc end kernels
        !$acc enter data copyin( nat, , q(:), dipm(:, :), at(:),&
        !$acc& qpint(:, :, :),xyz(:, :)),cn(:)

        !$acc parallel private(r_ab,result)

        !$acc loop gang vector collapse(2)
        do a = 1, nat
            do b = 1, nat
                !if (a == b) cycle 
                call inv_cn(nat,a,b,at,xyz,result)
                r_ab = xyz(:,a) - xyz(:,b)
                !$acc atomic
                delta_dipm(:,a) = delta_dipm(:,a) + (- r_ab(:) * q(b)) / (result*(cn(b)+1))
                !sorting of qp xx,xy,yy,xz,yz,zz
                !$acc atomic
                delta_qp(1,a) = delta_qp(1,a) + ( 1.5_wp*( r_ab(1)*r_ab(1)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(2,a) = delta_qp(2,a) + ( 1.5_wp*( r_ab(1)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(3,a) = delta_qp(3,a) + ( 1.5_wp*( r_ab(2)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(4,a) = delta_qp(4,a) + ( 1.5_wp*( r_ab(1)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(5,a) = delta_qp(5,a) + ( 1.5_wp*( r_ab(2)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(6,a) = delta_qp(6,a) + ( 1.5_wp*( r_ab(3)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !qp_part(:,a) = qp_part(:,a) + qp(:,b) / (result*(cn(b)+1))
            enddo
            !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
            !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
        enddo

        !$acc end parallel

        !$acc exit data copyout(delta_dipm(:, :), delta_qp(:, :), qp_part(:,:))
        call remove_trac_qp(nat,delta_qp,qp_part)

    end subroutine

    subroutine get_delta_mm_p(nat,q,dipm,qp,at,xyz,cn,delta_dipm,delta_qp) ! the sign of q was changed to respect the charge of the electrons
        implicit none 
        integer, INTENT(IN) :: nat, at(nat)
        real(wp), INTENT(IN) :: dipm(3,nat), xyz(3,nat), qp(6,nat),  q(nat), cn(nat)
        real(wp), INTENT(OUT) :: delta_dipm(3,nat), delta_qp(6,nat)
        integer :: a, b, i
        real(wp) :: result, r_ab(3), tii, qp_part(6,nat)

        !$acc enter data create(delta_dipm(:, :), delta_qp(:, :),qp_part(:,:))
        !$acc kernels default(present)
        delta_dipm = 0.0_wp
        delta_qp = 0.0_wp
        qp_part = 0.0_wp
        !$acc end kernels
        !$acc enter data copyin( nat, , q(:), dipm(:, :), at(:),&
        !$acc& qpint(:, :, :),xyz(:, :)),cn(:)

        !$acc parallel private(r_ab,result)

        !$acc loop gang vector collapse(2)

        do a = 1, nat
            do b = 1, nat
                !if (a == b) cycle 
                call inv_cn(nat,a,b,at,xyz,result)
                r_ab = xyz(:,a) - xyz(:,b)
                !$acc atomic
                delta_dipm(:,a) = delta_dipm(:,a) + (dipm(:,b) + r_ab(:) * q(b)) / (result*(cn(b)+1))
                !sorting of qp xx,xy,yy,xz,yz,zz
                !$acc atomic
                delta_qp(1,a) = delta_qp(1,a) + ( 1.5_wp*(-1*(r_ab(1)*dipm(1,b) + r_ab(1)*dipm(1,b)) - r_ab(1)*r_ab(1)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(2,a) = delta_qp(2,a) + ( 1.5_wp*(-1*(r_ab(1)*dipm(2,b) + r_ab(2)*dipm(1,b)) - r_ab(1)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(3,a) = delta_qp(3,a) + ( 1.5_wp*(-1*(r_ab(2)*dipm(2,b) + r_ab(2)*dipm(2,b)) - r_ab(2)*r_ab(2)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(4,a) = delta_qp(4,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(1,b) + r_ab(1)*dipm(3,b)) - r_ab(1)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(5,a) = delta_qp(5,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(2,b) + r_ab(2)*dipm(3,b)) - r_ab(2)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                !$acc atomic
                delta_qp(6,a) = delta_qp(6,a) + ( 1.5_wp*(-1*(r_ab(3)*dipm(3,b) + r_ab(3)*dipm(3,b)) - r_ab(3)*r_ab(3)*q(b))) / (result*(cn(b)+1))
                qp_part(:,a) = qp_part(:,a) + qp(:,b) / (result*(cn(b)+1))
            enddo
            !delta_dipm(:,a) = dipm(:,a) + delta_dipm(:,a)
            !delta_qp(:,a) = qp(:,a) + delta_qp(:,a)
        enddo
        !$acc end parallel

        !$acc exit data copyout(delta_dipm(:, :), delta_qp(:, :), qp_part(:,:))
        call remove_trac_qp(nat,delta_qp,qp_part)

    end subroutine

    subroutine remove_trac_qp(nat,qp_matrix,qp_part)
        implicit none
        integer, INTENT(IN) :: nat
        real(wp),INTENT(IN) :: qp_part(6,nat)
        real(wp), INTENT(INOUT) :: qp_matrix(6,nat)
        integer :: i
        real(wp) :: tii

        do i = 1, nat 
            tii = qp_matrix(1,i)+qp_matrix(3,i)+qp_matrix(6,i)
            tii = tii/3.0_wp
            !qp_matrix(1:6,i) = 1.50_wp*qp(1:6,i)
            qp_matrix(1,i) = qp_matrix(1,i)-tii
            qp_matrix(3,i) = qp_matrix(3,i)-tii
            qp_matrix(6,i) = qp_matrix(6,i)-tii
            qp_matrix(:,i) = qp_matrix(:,i) + qp_part(:,i)
         enddo
    end subroutine

    subroutine sum_up_mulliken(nat,nshell,aoat2,ash,mull_shell,mull_at)
        implicit none 
        integer, INTENT(IN) :: nat, nshell, aoat2(:), ash(:)
        real(wp), INTENT(IN) :: mull_shell(nshell)
        real(wp), INTENT(OUT) :: mull_at(nat)
        integer :: i

        mull_at = 0.0_wp
        
        do i = 1, nshell
            mull_at(ash(i)) = mull_at(ash(i)) + mull_shell(i)
        enddo

    end subroutine

    subroutine print_ML_output(out,nat,nao,nshell,at,id2at,ash,ao2at,res)
        integer, INTENT(IN) :: out,nat,nao,ao2at(nao), nshell,ash(nshell), at(nat),id2at(nat)
        type(results_type), intent(in) :: res
        integer :: i,j
        

        write(out,'(a)') 'Atom,w_xtb_tot,coordination number,delta coordination number,&
        mulliken charge s-shell,mulliken charge p-shell, mulliken charge d-shell,&
        dipm s-shell, dipm p-shell, dipm d-shell,&
        qm s-shell, qm p-shell, qm d-shell,&
        atomic partial charges,extended partial charges,dipm atomwise,extended dipm, &
        qm atomwise,extended qm,delta dipm only mull,delta qm only mull,delta dipm only Z,delta qm only Z ,&
        response (a.u.),gap (eV),chem.pot (eV),HOAO (eV),LUAO (eV),E_repulsion,E_EHT,&
        E_disp_2,E_disp_3,E_ies_ixc,E_aes,E_axc,E_tot'

        do i = 1, nat
            write(out,'(i2,a)', advance="no") at(id2at(i)),',' 
            write(out,'(f12.8,a)', advance="no") res%w_xtbml(i),','
            do j = 1, nfeatures-1
            write(out,'(f12.8,a)', advance="no") res%ml_features(i,j),',' 
            end do
            write(out,'(f14.8)') res%ml_features(i,nfeatures)
        enddo

        

    end subroutine

    subroutine mol_set_nuclear_charge(nat,at,id,z)
        implicit none
        integer,intent(in)::nat,at(nat),id(nat)
        real(wp), intent(out) :: z(nat)
        integer :: i
        do i = 1, nat
           z(i) = real(at(id(i)),wp) - real(ncore(at(id(i))))
           if (at(i) > 57 .and. at(i) < 72) z(i) = 3.0_wp
        enddo
     contains
     
     pure elemental integer function ncore(at)
       integer,intent(in) :: at
       if(at.le.2)then
          ncore=0
       elseif(at.le.10)then
          ncore=2
       elseif(at.le.18)then
          ncore=10
       elseif(at.le.29)then   !zn
          ncore=18
       elseif(at.le.36)then
          ncore=28
       elseif(at.le.47)then
          ncore=36
       elseif(at.le.54)then
          ncore=46
       elseif(at.le.71)then
          ncore=54
       elseif(at.le.79)then
          ncore=68
       elseif(at.le.86)then
          ncore=78
       endif
     end function ncore
     end subroutine mol_set_nuclear_charge

     subroutine pack_res(nat,nsh_tot,at2nsh,ml,e_tot,res)
        implicit none
        integer, intent(in) :: nat,nsh_tot,at2nsh(nat)
        real(wp), intent(in) :: e_tot
        type(results_type),intent(inout) :: res
        type(ml_features), intent(in) :: ml
        integer :: i, nsh
        allocate(res%ml_features(nat,nfeatures),source=0.0_wp)
        res%ml_features(:,1) = ml%cn_atom(:)
        res%ml_features(:,2) = ml%delta_cn(:)
        nsh = 1
        do i = 1,nat
            res%ml_features(i,3) = ml%shell_mulliken(nsh) !s shell always filled 
            res%ml_features(i,6) = ml%shell_dipm(nsh)
            res%ml_features(i,9) = ml%shell_qm(nsh)
            nsh = nsh + 1
            if (at2nsh(i) < 3) then
            res%ml_features(i,4) = ml%shell_mulliken(nsh)
            res%ml_features(i,7) = ml%shell_dipm(nsh)
            res%ml_features(i,10) = ml%shell_qm(nsh)
            nsh = nsh + 1
            else
            res%ml_features(i,4) = ml%shell_mulliken(nsh)
            res%ml_features(i,7) = ml%shell_dipm(nsh)
            res%ml_features(i,10) = ml%shell_qm(nsh)
            nsh = nsh + 1
            res%ml_features(i,5) = ml%shell_mulliken(nsh)
            res%ml_features(i,8) = ml%shell_dipm(nsh)
            res%ml_features(i,11) = ml%shell_qm(nsh)
            nsh = nsh + 1
            end if
        end do
        res%ml_features(:,12) = ml%atomic_partialcharge(:)
        res%ml_features(:,13) = ml%delta_partial_charge(:)
        res%ml_features(:,14) = ml%dipm_atom(:)
        res%ml_features(:,15) = ml%delta_dipm(:)
        res%ml_features(:,16) = ml%qm_atom(:)
        res%ml_features(:,17) = ml%delta_qm(:)
        res%ml_features(:,18) = ml%delta_dipm_only_p(:)
        res%ml_features(:,19) = ml%delta_qp_only_p(:)
        res%ml_features(:,20) = ml%delta_dipm_only_Z(:)
        res%ml_features(:,21) = ml%delta_qp_only_Z(:)
        res%ml_features(:,22) = ml%response(:)
        res%ml_features(:,23) = ml%hl_gap(:)
        res%ml_features(:,24) = ml%chempot(:)
        res%ml_features(:,25) = ml%ehoao(:)
        res%ml_features(:,26) = ml%eluao(:)
        res%ml_features(:,27) = ml%e_rep_atom(:)
        res%ml_features(:,28) = ml%e_EHT(:)
        res%ml_features(:,29) = ml%e_disp_2(:)
        res%ml_features(:,30) = ml%e_disp_3(:)
        res%ml_features(:,31) = ml%e_ies_ixc(:)
        res%ml_features(:,32) = ml%e_aes(:)
        res%ml_features(:,33) = ml%e_axc(:)
        res%ml_features(:,34) = e_tot
        
    end subroutine

end module xtbml_feature_calc