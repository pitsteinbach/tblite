module tblite_mulliken_kfock
    use mctc_env, only : wp, dp
    use mctc_io, only : structure_type
    use tblite_exchange_type, only : exchange_type
    use tblite_wavefunction_type, only : wavefunction_type
    use tblite_scf_potential, only : potential_type
    use tblite_container_cache, only : container_cache
    use tblite_exchange_cache, only : exchange_cache
    use tblite_basis_type, only : basis_type
    use iso_c_binding
    use tblite_timer, only : timer_type, format_time
    implicit none
    private


    logical :: allowincr = .true.
    logical :: debug = .false.
    public :: new_mulliken_exchange
    
    type, public, extends(exchange_type) :: mulliken_kfock_type
    !> allow single point precision, 0 false , 1 true
    integer :: allowsingle
    !> wether the fock matrix is build in incremental fashion, 0 false, 1 true
    integer :: incremental
    !> fullrange scale for the K scale
    real(wp) :: frscale, exchangescale = -0.5_wp
    real(wp), allocatable :: omega, lrscale
    !> Averaging scheme for hardness:
    !> 0 = arith.
    !> 1 = geom.
    !> 2 = harmonic
    integer :: average
    !> Smoothening exponent
    !> 1 =  Mataga
    !> 2 = Klopman-type
    integer :: expsmooth
    !> Chemical Hardness per shell with duplicates for reappearing atom types
    real(wp), allocatable :: hardness(:)
    
    type(c_ptr) :: indexer, model
    !> Compute gamma
    logical :: compute_gamma
    !> indexer object needed for computing gamma and computing grad 
    integer :: nao
    integer, allocatable :: ao2at(:)
    real(wp), allocatable :: gamma_(:,:)
    type(timer_type) :: timer
    !>
    contains
    !> Update container cache
    procedure :: update
    !> Get information about density dependent quantities used in the energy
    procedure :: variable_info
    !> Evaluate non-selfconsistent part of the interaction
    !procedure :: get_engrad
    !> Evaluate selfconsistent energy of the interaction
    procedure :: get_energy
    !> Evaluate denisty dependent potential
    procedure :: get_potential_w_overlap
    !> Information on container
    procedure :: info
    procedure :: delete
    end type
    
    character(len=*), parameter :: label = "Mulliken approximated exchange"

    interface
        function SetupIndexer(nat, nao, nsh, aonum, sh2at) result(indexer) bind(C, name="SetupIndexer")
            use iso_c_binding
            !> Number of atoms
            integer(c_int), value, intent(in) :: nat
            !> Number of AOs
            integer(c_int), value, intent(in) :: nao
            !> Number of shells
            integer(c_int), value, intent(in) :: nsh
            !> Convert from ao to at
            integer(c_int), intent(in) :: aonum(*)
            !> Convert from sh to at
            integer(c_int), intent(in) :: sh2at(*)
            type(c_ptr) :: indexer
        end function SetupIndexer
        subroutine DeleteIndexer(indexer) bind(C, name="DeleteIndexer")
            use iso_c_binding
            !> Indexer object
            type(c_ptr), value, intent(in) :: indexer
        end subroutine DeleteIndexer
        subroutine ComputeGamma(indexer, xyz, hardness, frscale, gammamat, average, expsmooth) bind(C, name="ComputeGamma")
            use iso_c_binding
            !> Indexer object
            type(c_ptr), value, intent(in) :: indexer
            !> Atomic coordinates
            real(c_double), intent(in) :: xyz(*)
            !> Chemical hardness
            real(c_double), intent(in) :: hardness(*)
            !> Full-range scale
            real(c_double), value, intent(in) :: frscale
            !> Gamma Values
            real(c_double), intent(out) :: gammamat(*)
            !> Averaging scheme for hardness:
            !> 0 = arith.
            !> 1 = geom.
            !> 2 = harmonic
            integer(c_size_t), value, intent(in) :: average
            !> Smoothening exponent
            !> 1 =  Mataga
            !> 2 = Klopman-type
            integer(c_size_t), value, intent(in) :: expsmooth
        end subroutine ComputeGamma
        function SetupExchangeSQM(indexer, allowsingle, incremental) result(exchange) bind(C, name="SetupExchangeSQM")
            use iso_c_binding
            !> Indexer object
            type(c_ptr), value, intent(in) :: indexer
            !> Allow single precision
            integer(c_int), value, intent(in) :: allowsingle
            !> Incremental Fock build
            integer(c_int), value, intent(in) :: incremental
            !> Exchange object
            type(c_ptr) :: exchange
        end function SetupExchangeSQM
        subroutine DeleteExchangeSQM(exchange) bind(C, name="DeleteExchangeSQM")
            use iso_c_binding
            !> Exchange object
            type(c_ptr), value, intent(in) :: exchange
        end subroutine DeleteExchangeSQM
        subroutine ComputeExchangeSQM(exchange, overlap, density, gamma_, fock) bind(C, name="ComputeExchangeSQM")
            use iso_c_binding
            !> Exchange object
            type(c_ptr), value, intent(in) :: exchange
            !> Gamma values
            real(c_double), intent(in) :: gamma_(*)
            !> Density matrix
            real(c_double), intent(in) :: density(*)
            !> Overlap matrix
            real(c_double), intent(in) :: overlap(*)
            !> Fock matrix
            real(c_double), intent(out) :: fock(*)
        end subroutine ComputeExchangeSQM
    end interface

    interface new_mulliken_exchange
        module procedure :: new_range_separated_mulliken_k_fock
    end interface
         
        contains
        
        pure function info(self, verbosity, indent) result(str)
        use tblite_output_format, only : format_string
        !> Instance of the interaction container
        class(mulliken_kfock_type), intent(in) :: self
        !> Verbosity level
        integer, intent(in) :: verbosity
        !> Indentation level
        character(len=*), intent(in) :: indent
        !> Information on the container
        character(len=:), allocatable :: str
        character(len=*), parameter :: nl = new_line('a')
        str = "Mulliken approximated semi-empirical exchange"//nl//indent//"Using the shell resolved chemical hardness"
        if (allocated(self%omega)) then
            str= str // nl//indent//"Range separted exchange is used:"// &
            & nl//indent//" * Full-range scale: "//format_string(self%frscale,'(f5.2)')//nl//indent//&
            & " * Omega: "//format_string(self%omega,'(f5.2)')//nl//indent//" * Long-range scale: "//format_string(self%lrscale,'(f5.2)')
        else
            str = str //nl//indent//"Full range exchange is used"// nl//indent//" * Full-range scale: "//format_string(self%frscale,'(f5.2)')
        end if
    end function info
    
    !> Update cache from container
    subroutine update(self, mol, cache)
        !> Instance of the multipole container
        class(mulliken_kfock_type), intent(in) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Reusable data container
        type(container_cache), intent(inout) :: cache
        
        type(exchange_cache), pointer :: ptr
        
        call taint(cache, ptr)
        
    end subroutine
    
    !> Inspect container cache and reallocate it in case of type mismatch
    subroutine taint(cache, ptr)
        !> Instance of the container cache
        type(container_cache), target, intent(inout) :: cache
        !> Reference to the container cache
        type(exchange_cache), pointer, intent(out) :: ptr
        
        if (allocated(cache%raw)) then
            call view(cache, ptr)
            if (associated(ptr)) return
            deallocate(cache%raw)
        end if
        
        if (.not.allocated(cache%raw)) then
            block
                type(exchange_cache), allocatable :: tmp
                allocate(tmp)
                call move_alloc(tmp, cache%raw)
                end block
            end if
            
            call view(cache, ptr)
        end subroutine taint
        
        !> Return reference to container cache after resolving its type
        subroutine view(cache, ptr)
            !> Instance of the container cache
            type(container_cache), target, intent(inout) :: cache
            !> Reference to the container cache
            type(exchange_cache), pointer, intent(out) :: ptr
            nullify(ptr)
            select type(target => cache%raw)
            type is(exchange_cache)
            ptr => target
        end select
    end subroutine view
    
    
    !> Create a new Mulliken approximated exchange container
    subroutine new_range_separated_mulliken_k_fock(self, mol, hardness, allowsingle, incremental, frscale, &
        & omega, lrscale, average, exp, bas)
        !> Instance of the multipole container
        type(mulliken_kfock_type), intent(out) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Chemical hardness sorted by shells
        real(wp), intent(in) :: hardness(:)
        !> allow single precision , incremental Fock Build
        logical,intent(in) :: allowsingle, incremental
        !> fullrange scale for K
        real(wp), intent(in) :: frscale
        !> omega if range seperated exchange is used
        real(wp), allocatable, intent(in) :: omega
        !> long range sclae for range seperated exchange treatment
        real(wp), allocatable, intent(in) :: lrscale
        !> Averaging scheme for hardness:
        !> 0 = arith.
        !> 1 = geom.
        !> 2 = arith.
        integer :: average
        !> Smoothening exponent
        !> 1 =  Mataga
        !> 2 = Klopman-type
        integer :: exp
        
        type(basis_type) :: bas
        
        call self%timer%push("Gamma")
        if (allowsingle .eqv. .true.) then
            self%allowsingle = 1
        else
            self%allowsingle= 0
        end if
        if (incremental .eqv. .true.) then
            self%incremental = 1
        else
            self%incremental= 0
        end if
        
        self%frscale = frscale
        if (allocated(omega)) then
            allocate(self%omega,self%lrscale)
            self%omega = omega
            self%lrscale = lrscale
        end if
        self%nao = bas%nao
        self%ao2at = bas%ao2at
        self%indexer = SetupIndexer(mol%nat, bas%nao, bas%nsh, bas%nao_sh, bas%sh2at)
        allocate(self%gamma_(bas%nao,bas%nao))
        call ComputeGamma(self%indexer, mol%xyz, hardness, self%frscale * self%exchangescale, self%gamma_, &
        int(average,kind=c_size_t), int(exp, kind=c_size_t))
        call self%timer%pop()
        call self%timer%push("Mulliken K-Fock")
        self%label = label
    
        self%model = SetupExchangeSQM(self%indexer, self%allowsingle, self%incremental)
        call self%timer%pop()
        
    end subroutine new_range_separated_mulliken_k_fock
    
    
    subroutine get_potential_w_overlap(self, mol, cache, wfn, pot, overlap)
        !> Instance of the exchange container
        class(mulliken_kfock_type), intent(in) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Tight-binding wavefunction data
        type(wavefunction_type), intent(in) :: wfn
        !> Density dependent potential shifts
        type(potential_type), intent(inout) :: pot
        !> Reusable data container
        type(container_cache), intent(inout) :: cache
        !> Overlap
        real(wp), intent(in) :: overlap(:,:)
        
        type(exchange_cache), pointer :: ptr
        integer :: spin
        
        call view(cache, ptr)
        
        if (.not.allocated(ptr%gamma_)) then
            allocate(ptr%gamma_(self%nao,self%nao))
            ptr%gamma_ = self%gamma_
        end if        
        !call self%timer%push("Mulliken K-Fock")
        
        !call self%timer%pop()
        if (self%incremental) then
            if (.not.allocated(ptr%ref_D)) then
                !!do full Fock build in double precision
                if (debug) write(*,*) "Full K Build!"
                if (.not.allocated(ptr%prev_F)) allocate(ptr%prev_F(self%nao, self%nao, wfn%nspin),source = 0.0_wp)
                if (.not.allocated(pot%kao)) allocate(pot%kao(self%nao,self%nao, wfn%nspin))
                do spin = 1, wfn%nspin
                    call ComputeExchangeSQM(self%model, overlap, wfn%density(:, : ,spin), ptr%gamma_, ptr%prev_F(: ,: ,spin))
                    pot%kao(:,:,spin) = ptr%prev_F(:,:,spin)
                end do
                
                if (.not.allocated(ptr%ref_D)) allocate(ptr%ref_D(self%nao, self%nao, wfn%nspin), source= wfn%density(:,:,:))
                if (.not.allocated(ptr%ref_F)) allocate(ptr%ref_F(self%nao, self%nao, wfn%nspin), source= ptr%prev_F(:,:,:))
                !if (.not.allocated(ptr%curr_D))allocate(ptr%curr_D(self%nao, self%nao, wfn%nspin), source= wfn%density(:,:,:))
                
            else
                ptr%curr_D = wfn%density(:,:,:) - ptr%ref_D
                ptr%prev_F(:, :, :) = ptr%ref_F(:, :, :)
                do spin = 1, wfn%nspin
                    call ComputeExchangeSQM(self%model, overlap, ptr%curr_D(:, : ,spin), ptr%gamma_, ptr%prev_F(: ,: ,spin))
                    !pot%kao(:, :, spin) = (ptr%prev_F(:, :, spin) + pot%kao(:, :, spin))
                    !ptr%prev_D = wfn%density(:,:,1)
                end do
                pot%kao(:, :, :) = ptr%prev_F(:, :, :)
                
            end if
        else
            if (debug) write(*,*) "Full K Build!"
            if (.not.allocated(ptr%prev_F)) allocate(ptr%prev_F(self%nao, self%nao, wfn%nspin),source = 0.0_wp)
            if (.not.allocated(pot%kao)) allocate(pot%kao(self%nao,self%nao, wfn%nspin))
            
            do spin = 1, wfn%nspin
                call ComputeExchangeSQM(self%model, overlap, wfn%density(:, : ,spin), ptr%gamma_, ptr%prev_F(: ,: ,spin))
        
                pot%kao(:, :, spin) = ptr%prev_F(:, :, spin)
            end do

            !!if (allocated(ptr%curr_D)) allocate(ptr%ref_D(self%nao, self%nao), source= wfn%density(:,:,1))
            !if (.not.allocated(ptr%curr_D))allocate(ptr%curr_D(self%nao, self%nao), source= wfn%density(:,:,1))
        endif
        
    end subroutine
    
    subroutine get_energy(self, mol, cache, wfn, energies)
        use tblite_blas, only : gemm
        !> Instance of the exchange container
        class(mulliken_kfock_type), intent(in) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction data
        type(wavefunction_type), intent(in) :: wfn
        !> Electrostatic energy
        real(wp), intent(inout) :: energies(:)
        !> Reusable data container
        type(container_cache), intent(inout) :: cache
        type(exchange_cache), pointer :: ptr
        integer :: iao, jao, spin

        real(wp) :: exchange(mol%nat)
        exchange=0.0_wp

        call view(cache, ptr)
        !$omp parallel do collapse(3) schedule(runtime) default(none) &
        !$omp reduction(+:exchange) shared(ptr, wfn, self) private(spin, iao, jao)
        do spin = 1, size(wfn%density, 3)
             do iao = 1, size(wfn%density, 2)
                 do jao = 1, size(wfn%density, 1)
                     exchange(self%ao2at(iao)) = exchange(self%ao2at(iao)) + ptr%prev_F(jao, iao, spin) * wfn%density(jao, iao, spin)
                 end do
             end do
         end do

         if (size(wfn%density,3).eq.2) then
            do iao = 1, size(wfn%density, 2)
                 do jao = 1, size(wfn%density, 1)
                     exchange(self%ao2at(iao)) = exchange(self%ao2at(iao)) + ptr%prev_F(jao, iao, 1) * wfn%density(jao, iao, 2) &
                     & + ptr%prev_F(jao, iao, 2) * wfn%density(jao, iao, 1)
                end do
            end do
        end if
        
        energies(:) = energies(:) + exchange(:) * 0.5

    end subroutine get_energy
    
    subroutine get_gradient_w_overlap(self, mol, cache, wfn, gradient, ao_grad, overlap)
        !> Instance of the exchange container
        class(mulliken_kfock_type), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Wavefunction data
        type(wavefunction_type), intent(in) :: wfn
        !> Reusable data container
        type(container_cache), intent(inout) :: cache
        !> Molecular gradient of the exchange energy
        real(wp), contiguous, intent(inout) :: gradient(:, :)
        !> Strain derivatives of the exchange energy
        real(wp), contiguous, intent(inout) :: ao_grad(:, :)
        !> Overlap
        real(wp), intent(in) :: overlap(:,:)
        real(wp), allocatable :: intermediate(:, :)
        real(wp), allocatable :: grad_before(:, :)
        
        type(exchange_cache), pointer :: ptr
        integer :: spin
        
        call view(cache, ptr)
        allocate(intermediate(self%nao, self%nao), source=0.0_wp)
        allocate(grad_before(self%nao, self%nao), source=gradient)
        !do spin = 1, wfn%nspin
            !call KGradSymSQM(0, self%nao, ptr%gamma_, wfn%density(:, :, spin), overlap, ao_grad, intermediate)
        !end do
        
        !!if (allocated(self%omega)) then 
            !!call GammaGradSQM_rs(self%indexer, mol%xyz, self%hardness, self%frscale * self%exchangescale, self%omega, self%lrscale, intermediate, gradient)
        !!else
            !!call GammaGradSQM_fr(self%indexer, mol%xyz, self%hardness, self%frscale * self%exchangescale, intermediate, gradient)
        !!endif
        
        
        
    end subroutine get_gradient_w_overlap
    
pure function variable_info(self) result(info)
    use tblite_scf_info, only : scf_info, atom_resolved, orbital_resolved, not_used
    !> Instance of the electrostatic container
    class(mulliken_kfock_type), intent(in) :: self
    !> Information on the required potential data
    type(scf_info) :: info
    
    info = scf_info(charge=not_used, dipole=not_used, quadrupole=not_used, density=orbital_resolved)
end function variable_info

subroutine delete(self)
    !> Instance of the exchange container
    class(mulliken_kfock_type), intent(inout) :: self
    !> Delete indexer object
    call DeleteIndexer(self%indexer)
    !> Delete exchange object
    call DeleteExchangeSQM(self%model)
    block
        integer :: it
        real(wp) :: stime
        real(wp) :: ttime = 0.0_wp
        character(len=*), parameter :: label(*) = [character(len=20):: &
        & "Mulliken K-Fock", "Gamma"]
        do it = 1, size(label)
           stime = self%timer%get(label(it))
           ttime = ttime + stime
           if (stime <= epsilon(0.0_wp)) cycle
           write(*,*) " - "//label(it)//format_time(stime)
        end do
        write(*,*) "___________________________________"             
        write(*,*) "   "//"Total: "//format_time(ttime)
    end block
    !> Deallocate gamma matrix
    if (allocated(self%gamma_)) deallocate(self%gamma_)
    !> Deallocate omega and lrscale
    if (allocated(self%omega)) deallocate(self%omega)
    if (allocated(self%lrscale)) deallocate(self%lrscale)
end subroutine delete


end module tblite_mulliken_kfock
