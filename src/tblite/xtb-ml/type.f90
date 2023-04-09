module xtb_type_ml
    use mctc_env, only : wp
  
    implicit none
 
    public :: ml_features
 
    private
 
    type :: ml_features
        real(wp),ALLOCATABLE ::  w_1e (:)
        real(wp),ALLOCATABLE ::  w_dif  (:)
        real(wp),ALLOCATABLE ::  w_tot  (:)
        real(wp),ALLOCATABLE ::  cn_atom  (:) ! get it from cn in singlepoint
        real(wp),ALLOCATABLE ::  delta_cn  (:)
        real(wp),ALLOCATABLE ::  shell_mulliken (:)
        real(wp),ALLOCATABLE ::  atomic_partialcharge  (:)
        real(wp),ALLOCATABLE ::  delta_partial_charge  (:)
        real(wp),ALLOCATABLE ::  shell_dipm (:)
        !shell xyz
        real(wp),ALLOCATABLE ::  shell_dipm_xyz (:,:)
        !
        real(wp),ALLOCATABLE ::  dipm_atom  (:)
        !dipm xyz
        real(wp),ALLOCATABLE ::  dipm_atom_xyz  (:,:)

        real(wp),ALLOCATABLE ::  delta_dipm  (:)
        real(wp),ALLOCATABLE ::  shell_qm (:)
        real(wp),ALLOCATABLE ::  qm_atom  (:)
        real(wp),ALLOCATABLE ::  delta_qm  (:)
        real(wp),ALLOCATABLE ::  response  (:)
        real(wp),ALLOCATABLE ::  egap  (:)
        real(wp),ALLOCATABLE ::  chempot  (:)
        real(wp),ALLOCATABLE ::  ehoao  (:)
        real(wp),ALLOCATABLE ::  eluao  (:)
        real(wp),ALLOCATABLE ::  e_rep_atom  (:)
        real(wp),ALLOCATABLE ::  e_EHT  (:)
        real(wp),ALLOCATABLE ::  e_disp_2  (:)
        real(wp),ALLOCATABLE ::  e_disp_3  (:)
        real(wp),ALLOCATABLE ::  e_ies_ixc  (:)
        real(wp),ALLOCATABLE ::  e_aes  (:)
        real(wp),ALLOCATABLE ::  delta_dipm_only_p  (:)
        real(wp),ALLOCATABLE ::  delta_qp_only_p  (:)
        real(wp),ALLOCATABLE ::  delta_dipm_only_Z  (:)
        real(wp),ALLOCATABLE ::  delta_qp_only_Z  (:)
        real(wp),ALLOCATABLE ::  hl_gap (:)
        real(wp),ALLOCATABLE ::  fermi_level (:)
        real(wp),ALLOCATABLE ::  e_axc (:)
        !energy based features; extensions
        real(wp),ALLOCATABLE ::  chempot_ext(:)
        real(wp),ALLOCATABLE ::  e_gap_ext(:)
        real(wp),ALLOCATABLE ::  eluao_ext(:)
        real(wp),ALLOCATABLE ::  ehoao_ext(:)
    contains
        procedure :: allocate => allocate_ml
       ! procedure :: deallocate => deallocate_ml
    end type ml_features
 
 contains
 
 subroutine allocate_ml(self,nat,nshell)
    implicit none
    class(ml_features) :: self
    integer,intent(in)  :: nat,nshell
    allocate( self%w_1e (nat),      source = 0.0_wp )
    allocate( self%w_dif (nat),     source = 0.0_wp )
    allocate( self%w_tot (nat),     source = 0.0_wp )
    allocate( self%cn_atom (nat),   source = 0.0_wp )
    allocate( self%delta_cn (nat),  source = 0.0_wp )
    allocate( self%atomic_partialcharge(nat),   source = 0.0_wp )
    allocate( self%delta_partial_charge (nat),  source = 0.0_wp )
    allocate( self%dipm_atom (nat), source = 0.0_wp )

    allocate( self%dipm_atom_xyz (3,nat), source = 0.0_wp )

    allocate( self%delta_dipm (nat),source = 0.0_wp )
    allocate( self%qm_atom (nat),   source = 0.0_wp )
    allocate( self%delta_qm (nat),  source = 0.0_wp )
    allocate( self%response (nat),  source = 0.0_wp )
    allocate( self%egap (nat),      source = 0.0_wp )
    allocate( self%chempot (nat),   source = 0.0_wp )
    allocate( self%ehoao (nat),     source = 0.0_wp )
    allocate( self%eluao (nat),     source = 0.0_wp )
    allocate( self%e_rep_atom (nat),source = 0.0_wp )
    allocate( self%e_EHT (nat),     source = 0.0_wp )
    allocate( self%e_disp_2 (nat),source = 0.0_wp )
    allocate( self%e_disp_3 (nat),  source = 0.0_wp )
    allocate( self%e_ies_ixc (nat), source = 0.0_wp )
    allocate( self%e_aes (nat),     source = 0.0_wp )
    allocate( self%delta_dipm_only_p (nat),     source = 0.0_wp )
    allocate( self%delta_qp_only_p (nat),       source = 0.0_wp )
    allocate( self%delta_dipm_only_Z (nat),     source = 0.0_wp )
    allocate( self%delta_qp_only_Z (nat),       source = 0.0_wp )
    allocate( self%e_axc (nat),     source = 0.0_wp )

    allocate( self%shell_dipm (nshell),     source = 0.0_wp )
    !shell xyz
    allocate( self%shell_dipm_xyz (3,nshell),     source = 0.0_wp )

    allocate( self%shell_mulliken (nshell), source = 0.0_wp )
    allocate( self%shell_qm (nshell),       source = 0.0_wp )
    !Extensions
    allocate( self%chempot_ext(nat),       source = 0.0_wp )
    allocate( self%e_gap_ext(nat),       source = 0.0_wp )
    allocate( self%eluao_ext(nat),       source = 0.0_wp )
    allocate( self%ehoao_ext(nat),       source = 0.0_wp )

    


 end subroutine allocate_ml
 end module xtb_type_ml
 