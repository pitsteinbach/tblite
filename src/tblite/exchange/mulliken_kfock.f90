module tblite_mulliken_kfock
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_exchange_type, only : exchange_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_scf_potential, only : potential_type
   use tblite_container_cache, only : container_cache
   use tblite_exchange_cache, only : exchange_cache
   use tblite_basis_type, only : basis_type
   use mullk_fockbuild, only : compute_gamma_fr, compute_gamma_rs, KFockSymSQM
   implicit none
   private
   logical :: allowincr = .true.
   public :: new_mulliken_exchange

   type,public, extends(exchange_type) :: mulliken_kfock_type
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
      !> 2 = arith.
      integer :: average
      !> Smoothening exponent
      !> 1 =  Mataga
      !> 2 = Klopman-type
      integer :: expsmooth
      !> Chemical Hardness per shell with duplicates for reappearing atom types
      real(wp), allocatable :: hardness(:)
      !> number of AOS
      integer :: nao
      !> number of shells
      integer :: nsh
      !> Number of AOs in shell
      integer, allocatable :: aonum(:)
      !> Convert from sh to at
      integer, allocatable :: sh2at(:)
      !> Convert from ao to at
      integer, allocatable :: ao2at(:)
      !> Compute gamma
      logical :: compute_gamma
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
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient_w_overlap
      !> Information on container
      procedure :: info
   end type

   character(len=*), parameter :: label = "Mulliken approximated exchange"

   interface new_mulliken_exchange
      module procedure :: new_range_separated_mulliken_k_fock
   end interface


   interface gamma
         procedure :: compute_gamma_fr
         procedure :: compute_gamma_rs
   end interface gamma

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
      & omega, lrscale, average, exp, incr,bas)
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
   logical :: incr
   type(basis_type) :: bas
    allowincr = incr
   self%nao = bas%nao
   allocate(self%aonum(bas%nsh), self%sh2at(bas%nsh),self%ao2at(self%nao))
   self%aonum = bas%nao_sh
   self%sh2at = bas%sh2at
   self%nsh = bas%nsh
   self%ao2at = bas%ao2at

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

   self%average = average
   self%expsmooth = exp
   self%hardness = hardness

   self%label = label

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

   call view(cache, ptr)

   if (.not.allocated(ptr%gamma_)) then
      allocate(ptr%gamma_(self%nao,self%nao))
      if (allocated(self%omega)) then
         call gamma(self%nao, mol%nat, self%nsh, self%aonum, self%sh2at,&
            & self%average, self%expsmooth, mol%xyz, self%hardness, self%frscale * self%exchangescale, self%omega, self%frscale * self%exchangescale, ptr%gamma_)
      else
         call  gamma(self%nao, mol%nat, self%nsh, self%aonum, self%sh2at,self%average, self%expsmooth, mol%xyz, self%hardness, self%frscale * self%exchangescale , ptr%gamma_)
      end if
   end if

   if (allowincr) then
      if (.not.allocated(ptr%ref_D)) then
         !do full Fock build in double precision
         write(*,*) "Full K Build!"
         if (.not.allocated(ptr%prev_F)) allocate(ptr%prev_F(self%nao, self%nao),source = 0.0_wp)
         call KFockSymSQM(0, 0,self%nao, mol%nat, self%nsh, self%aonum, self%sh2at, &
            & ptr%gamma_, wfn%density(:, :, 1), overlap, ptr%prev_F)
         if (.not.allocated(pot%kao)) allocate(pot%kao(self%nao,self%nao, 1))
         pot%kao(:,:,1) = ptr%prev_F

         if (allocated(ptr%curr_D)) allocate(ptr%ref_D(self%nao, self%nao), source= wfn%density(:,:,1))
         if (.not.allocated(ptr%curr_D))allocate(ptr%curr_D(self%nao, self%nao), source= wfn%density(:,:,1))

      else
         ptr%curr_D = wfn%density(:,:,1) - ptr%ref_D
         pot%kao(:, :, 1) = ptr%prev_F
         call KFockSymSQM(1, 0,self%nao, mol%nat, self%nsh, self%aonum, self%sh2at, &
            & ptr%gamma_, ptr%curr_D, overlap, pot%kao(:, :, 1))
         pot%kao(:, :, 1) = (ptr%prev_F + pot%kao(:, :, 1))
         !ptr%prev_D = wfn%density(:,:,1)


      end if
   else
      write(*,*) "Full K Build!"
      if (.not.allocated(ptr%prev_F)) allocate(ptr%prev_F(self%nao, self%nao),source = 0.0_wp)
      call KFockSymSQM(0, 0,self%nao, mol%nat, self%nsh, self%aonum, self%sh2at, &
         & ptr%gamma_, wfn%density(:, :, 1), overlap, ptr%prev_F)
      if (.not.allocated(pot%kao)) allocate(pot%kao(self%nao,self%nao, 1))
      pot%kao(:,:,1) =  ptr%prev_F
      !if (allocated(ptr%curr_D)) allocate(ptr%ref_D(self%nao, self%nao), source= wfn%density(:,:,1))
      !if (.not.allocated(ptr%curr_D))allocate(ptr%curr_D(self%nao, self%nao), source= wfn%density(:,:,1))
   endif

end subroutine

subroutine get_energy(self, mol, cache, wfn ,energies)
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

   call view(cache, ptr)
   !!$omp parallel do collapse(3) schedule(runtime) default(none) &
   !!$omp reduction(+:energies) shared(ptr, wfn) private(spin, iao, jao)
   do spin = 1, size(wfn%density, 3)
      do iao = 1, size(wfn%density, 2)
         do jao = 1, size(wfn%density, 1)
            energies(self%ao2at(iao)) = energies(self%ao2at(iao)) + ptr%prev_F(jao, iao) * wfn%density(jao, iao, spin) * 0.5
         end do
      end do
   end do

end subroutine get_energy

subroutine get_gradient_w_overlap(self, mol, cache, wfn, gradient, sigma, overlap)
   !> Instance of the exchange container
   class(mulliken_kfock_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)
   !> Overlap
   real(wp), intent(in) :: overlap(:,:)

   type(exchange_cache), pointer :: ptr
   integer :: i

   call view(cache, ptr)

   call KGradSymSQM(0, self%nao, mol%nat, self%nsh, self%aonum, self%sh2at, ptr%gamma_, &
      & wfn%density(:, :, 1), overlap, sigma, gradient)



end subroutine get_gradient_w_overlap

pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved, orbital_resolved, not_used
   !> Instance of the electrostatic container
   class(mulliken_kfock_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=not_used ,dipole=not_used, quadrupole=not_used, density=orbital_resolved)
end function variable_info

end module tblite_mulliken_kfock
