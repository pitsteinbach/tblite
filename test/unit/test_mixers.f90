module test_mixers
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mctc_io_structure, only : new
   use tblite_context_type, only : context_type
   use tblite_results, only : results_type
   use tblite_scf_mixer_input, only : mixer_input, mixer_kind, mixer_precision, mixer_runmode
   use tblite_wavefunction , only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private
   real(wp), parameter :: acc = 0.1_wp
   real(wp), parameter :: thr = 1.e-6_wp*acc
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   public :: collect_mixers
contains

subroutine collect_mixers(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("mixers-broyden", test_broyden),&
      new_unittest("test-diis-cpu", test_diis_cpu)]

end subroutine collect_mixers

subroutine test_broyden(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(mixer_input) :: mixer_config
   type(results_type) :: res

   integer, parameter :: nat=12
   real(wp) :: energy_native = 0.0_wp, energy_gambits = 0.0_wp
   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.06880660023529,       -0.46478030005927,        0.00009471732781,&
   &2.45331325533661,       -0.46484444142679,        0.00042084312846,&
   &3.14559996373466,        0.73409173401801,        0.00008724271521,&
   &2.45340597395333,        1.93314591537421,       -0.00086301044874,&
   &1.06895328716368,        1.93321130458200,       -0.00141386731889,&
   &0.37663958202671,        0.73422879200405,       -0.00090929198808,&
   &0.52864276199175,       -1.40035288735680,        0.00049380958906,&
   &2.99337903563419,       -1.40047547903112,        0.00121759015506,&
   &4.22591259698321,        0.73408789322206,        0.00025398801936,&
   &2.99365942822711,        2.86866756976346,       -0.00166131228918,&
   &0.52879830433456,        2.86879139255056,       -0.00224874122149,&
   &-0.70367266962110,        0.73433126635962,       -0.00138296766859&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/6,6,6,6,6,6,1,1,1,1,1,1/)

   call new(mol, num, xyz*aatoau, uhf=0, charge=0.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   
   mixer_config%kind = mixer_kind%broyden
   mixer_config%memory = 250
   mixer_config%nao = calc%bas%nao
   mixer_config%prec = mixer_precision%double
   mixer_config%damp = 0.4_wp

   calc%mixer_info = mixer_config
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy_native, verbosity=0, results=res)

   call new(mol, num, xyz*aatoau, uhf=0, charge=0.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   mixer_config%kind = mixer_kind%gambits_broyden
   calc%mixer_info = mixer_config
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy_gambits, verbosity=0, results=res)

   if (abs(energy_native - energy_gambits) > thr) then
      call test_failed(error, "GAMBITS Broyden mixing does not give the same energy as native Broyden mixing.")
      print '(2es21.14)', energy_gambits, energy_native
   end if

end subroutine test_broyden

subroutine test_diis_cpu(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(mixer_input) :: mixer_config
   type(results_type) :: res

   integer, parameter :: nat=22
   real(wp) :: energy_native = 0.0_wp, energy_gambits = 0.0_wp
   real(wp), parameter :: xyz(3, nat) = reshape((/&
   &1.40704587900135,-1.266053426,-1.93713467409179,&
   &1.8500720142163,-0.46824073,-1.50918243052625,&
   &-0.03362432546857,-1.392692458,-1.74003582842655,&
   &-0.56857010176787,-1.017644449,-2.61263468250045,&
   &-0.44096297533149,-2.843378101,-1.4889973466575,&
   &-0.47991761435962,-0.552309546,-0.55520223211488,&
   &-1.51566046566003,-2.891873561,-1.32273881899144,&
   &-0.18116520826015,-3.451878075,-2.34920432497853,&
   &0.06989722371032,-3.232990003,-0.60872832970057,&
   &-1.56668254604022,0.00552121,-0.52884675232746,&
   &1.99245341935793,-1.73097166,-3.08869240465405,&
   &3.42884245712259,-1.306600699,-3.28712667180899,&
   &3.8772196423657,-0.888431234,-2.38921454082853,&
   &3.4654854727687,-0.564953085,-4.0831179008844,&
   &4.00253375919125,-2.169709391,-3.61210069945494,&
   &1.40187969243713,-2.438261129,-3.89034129099619,&
   &0.40869198564818,-0.491017096,0.47992425165481,&
   &1.15591901840578,-1.165248428,0.48740266863377,&
   &0.00723492497865,0.116922762,1.73426298331318,&
   &0.88822128835955,0.28499002,2.34645659039969,&
   &-0.47231557974936,1.067376345,1.52286683213051,&
   &-0.70199988222212,-0.504859383,2.28058248842892&
   &/),shape=(/3,nat/))
   integer, parameter :: num(nat) = (/7,1,6,1,6,6,1,1,1,8,6,6,1,1,1,8,7,1,6,1,1,1/)

   call new(mol, num, xyz*aatoau, uhf=0, charge=0.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   
   mixer_config%kind = mixer_kind%broyden
   mixer_config%memory = 250
   mixer_config%nao = calc%bas%nao
   mixer_config%prec = mixer_precision%double
   mixer_config%damp = 0.4_wp

   calc%mixer_info = mixer_config
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy_native, verbosity=0, results=res)

   call new(mol, num, xyz*aatoau, uhf=0, charge=0.0_wp)
   call new_gfn2_calculator(calc, mol, error)
   
   mixer_config%kind = mixer_kind%gambits_diis
   mixer_config%memory = 5
   mixer_config%prec = mixer_precision%double
   calc%mixer_info = mixer_config
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy_gambits, verbosity=0, results=res)

   if (abs(energy_native - energy_gambits) > thr) then
      call test_failed(error, "GAMBITS DIIS mixing does not give the same energy as native Broyden mixing.")
      print '(2es21.14)', energy_gambits, energy_native
      write(*,*) "Energy difference:", abs(energy_gambits - energy_native)
   end if

end subroutine test_diis_cpu

end module test_mixers
