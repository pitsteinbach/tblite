! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/xtb/gfn2.f90
!> Provides the parametrization for the GFN2-xTB Hamiltonian

!> Implementation of the GFN2-xTB Hamiltonian to parametrize an xTB calculator.
module tblite_sqmbox_gfn2
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_symbols, only : to_symbol
   use tblite_basis_type, only : basis_type, new_basis, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_coulomb_charge, only : new_effective_coulomb, effective_coulomb, &
      & arithmetic_average, coulomb_kernel
   use tblite_coulomb_multipole, only : new_damped_multipole
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_disp, only : d4_dispersion, new_d4_dispersion
   use tblite_ncoord, only : new_ncoord
   use tblite_param, only : param_record
   use tblite_repulsion, only : new_repulsion
   use tblite_sqmbox_calculator, only : sqmbox_calculator
   use tblite_xtb_h0, only : new_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_gfn2_sqmbox_calculator
   public :: export_gfn2_param


   integer, parameter :: max_elem = 86
   integer, parameter :: max_shell = 3

   !> Use older eV to Eh conversion for consistency
   real(wp), parameter :: evtoau = 1.0_wp / 27.21138505_wp

   !> Dispersion parameters
   real(wp), parameter :: s6 = 1.0_wp, s8 = 2.7_wp, a1 = 0.52_wp, a2 = 5.0_wp, s9 = 5.0_wp

   !> Repulsion parameters
   real(wp), parameter :: rep_kexp = 1.5_wp, rep_kexp_light = 1.0_wp, rep_rexp = 1.0_wp

   real(wp), parameter :: wexp = 0.5_wp
   real(wp), parameter :: enscale = 2.0e-2_wp
   real(wp), parameter :: kdiag(0:4) = [1.85_wp, 2.23_wp, spread(2.23_wp, 1, 3)]

   !> Exponents of repulsion term for GFN2-xTB repulsion
   real(wp), parameter :: rep_alpha(max_elem) = [&
      & 2.213717_wp, 3.604670_wp, 0.475307_wp, 0.939696_wp, 1.373856_wp, &
      & 1.247655_wp, 1.682689_wp, 2.165712_wp, 2.421394_wp, 3.318479_wp, &
      & 0.572728_wp, 0.917975_wp, 0.876623_wp, 1.187323_wp, 1.143343_wp, &
      & 1.214553_wp, 1.577144_wp, 0.896198_wp, 0.482206_wp, 0.683051_wp, &
      & 0.574299_wp, 0.723104_wp, 0.928532_wp, 0.966993_wp, 1.071100_wp, &
      & 1.113422_wp, 1.241717_wp, 1.077516_wp, 0.998768_wp, 1.160262_wp, &
      & 1.122923_wp, 1.222349_wp, 1.249372_wp, 1.230284_wp, 1.296174_wp, &
      & 0.908074_wp, 0.574054_wp, 0.697345_wp, 0.706172_wp, 0.681106_wp, &
      & 0.865552_wp, 1.034519_wp, 1.019565_wp, 1.031669_wp, 1.094599_wp, &
      & 1.092745_wp, 0.678344_wp, 0.936236_wp, 1.024007_wp, 1.139959_wp, &
      & 1.122937_wp, 1.000712_wp, 1.017946_wp, 1.012036_wp, 0.585257_wp, &
      & 0.716259_wp, 0.737643_wp, 0.729950_wp, 0.734624_wp, 0.739299_wp, &
      & 0.743973_wp, 0.748648_wp, 0.753322_wp, 0.757996_wp, 0.762671_wp, &
      & 0.767345_wp, 0.772020_wp, 0.776694_wp, 0.781368_wp, 0.786043_wp, &
      & 0.790717_wp, 0.852852_wp, 0.990234_wp, 1.018805_wp, 1.170412_wp, &
      & 1.221937_wp, 1.197148_wp, 1.204081_wp, 0.919210_wp, 1.137360_wp, &
      & 1.399312_wp, 1.179922_wp, 1.130860_wp, 0.957939_wp, 0.963878_wp, &
      & 0.965577_wp]

   !> Effective nuclear charge for GFN2-xTB repulsion
   real(wp), parameter :: rep_zeff(max_elem) = [&
      &  1.105388_wp,  1.094283_wp,  1.289367_wp,  4.221216_wp,  7.192431_wp, &
      &  4.231078_wp,  5.242592_wp,  5.784415_wp,  7.021486_wp, 11.041068_wp, &
      &  5.244917_wp, 18.083164_wp, 17.867328_wp, 40.001111_wp, 19.683502_wp, &
      & 14.995090_wp, 17.353134_wp,  7.266606_wp, 10.439482_wp, 14.786701_wp, &
      &  8.004267_wp, 12.036336_wp, 15.677873_wp, 19.517914_wp, 18.760605_wp, &
      & 20.360089_wp, 27.127744_wp, 10.533269_wp,  9.913846_wp, 22.099503_wp, &
      & 31.146750_wp, 42.100144_wp, 39.147587_wp, 27.426779_wp, 32.845361_wp, &
      & 17.363803_wp, 44.338211_wp, 34.365525_wp, 17.326237_wp, 24.263093_wp, &
      & 30.562732_wp, 48.312796_wp, 44.779882_wp, 28.070247_wp, 38.035941_wp, &
      & 28.674700_wp,  6.493286_wp, 26.226628_wp, 63.854240_wp, 80.053438_wp, &
      & 77.057560_wp, 48.614745_wp, 63.319176_wp, 51.188398_wp, 67.249039_wp, &
      & 46.984607_wp, 50.927529_wp, 48.676714_wp, 47.669448_wp, 46.662183_wp, &
      & 45.654917_wp, 44.647651_wp, 43.640385_wp, 42.633120_wp, 41.625854_wp, &
      & 40.618588_wp, 39.611322_wp, 38.604057_wp, 37.596791_wp, 36.589525_wp, &
      & 35.582259_wp, 40.186772_wp, 54.666156_wp, 55.899801_wp, 80.410086_wp, &
      & 62.809871_wp, 56.045639_wp, 53.881425_wp, 14.711475_wp, 51.577544_wp, &
      & 58.801614_wp,102.368258_wp,132.896832_wp, 52.301232_wp, 81.771063_wp, &
      &128.133580_wp]

   real(wp), parameter :: gexp = 2.0_wp

   real(wp), parameter :: hubbard_parameter(max_elem) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp, &
      & 0.472633_wp, 0.513586_wp, 0.589187_wp, 0.396299_wp, 0.346651_wp, &
      & 0.271594_wp, 0.477760_wp, 0.344970_wp, 0.202969_wp, 0.564152_wp, &
      & 0.432236_wp, 0.802051_wp, 0.571748_wp, 0.235052_wp, 0.261253_wp, &
      & 0.424373_wp, 0.210481_wp, 0.340000_wp, 0.711958_wp, 0.461440_wp, &
      & 0.952957_wp, 0.586134_wp, 0.368054_wp, 0.711205_wp, 0.509183_wp, &
      & 0.273310_wp, 0.263740_wp, 0.392012_wp, 0.461812_wp, 0.900000_wp, &
      & 0.942294_wp, 0.750000_wp, 0.383124_wp, 0.424164_wp, 0.236569_wp, &
      & 0.245937_wp, 0.597716_wp, 0.662889_wp, 0.660710_wp, 0.658531_wp, &
      & 0.656352_wp, 0.654173_wp, 0.651994_wp, 0.649815_wp, 0.647635_wp, &
      & 0.645456_wp, 0.643277_wp, 0.641098_wp, 0.638919_wp, 0.636740_wp, &
      & 0.634561_wp, 0.662597_wp, 0.449812_wp, 0.685426_wp, 0.224623_wp, &
      & 0.364388_wp, 0.548507_wp, 0.353574_wp, 0.438997_wp, 0.457611_wp, &
      & 0.418841_wp, 0.168152_wp, 0.900000_wp, 1.023267_wp, 0.288848_wp, &
      & 0.303400_wp]

   real(wp), parameter :: shell_hubbard(0:2, max_elem) = 1.0_wp + reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp, &
      & 0.0_wp,-0.0800000_wp,-0.2046716_wp, 0.0_wp,-0.3800000_wp,-0.4921114_wp, &
      & 0.0_wp,-0.4500000_wp,-0.0379088_wp, 0.0_wp,-0.4700000_wp, 0.7405872_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0545811_wp, 0.0_wp,-0.6500000_wp, 0.4046615_wp, &
      & 0.0_wp,-0.6500000_wp,-0.2418493_wp, 0.0_wp,-0.6000000_wp,-0.0611188_wp, &
      & 0.0_wp, 0.0700000_wp, 1.3333066_wp, 0.0_wp, 0.0684343_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5416555_wp,-0.3000000_wp, 0.0_wp,-0.3809089_wp,-0.1500000_wp, &
      & 0.0_wp,-0.4104743_wp,-0.5000000_wp, 0.0_wp, 0.1192113_wp,-0.2500000_wp, &
      & 0.0_wp, 0.5203002_wp, 0.4000000_wp, 0.0_wp,-0.2503223_wp,-0.0700000_wp, &
      & 0.0_wp, 0.9386493_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp, &
      & 0.0_wp,-0.4500000_wp,-0.3349288_wp, 0.0_wp,-0.1100000_wp,-0.4422630_wp, &
      & 0.0_wp,-0.0500000_wp,-0.3562950_wp, 0.0_wp,-0.3000000_wp,-0.4301371_wp, &
      & 0.0_wp,-0.6000000_wp, 0.3956819_wp, 0.0_wp,-0.6500000_wp,-0.3052305_wp, &
      & 0.0_wp,-0.6500000_wp,-0.1881774_wp, 0.0_wp,-0.6000000_wp, 0.0931707_wp, &
      & 0.0_wp,-0.0300000_wp, 0.8024848_wp, 0.0_wp, 0.2388669_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5867460_wp,-0.2800000_wp, 0.0_wp,-0.5090746_wp,-0.0600000_wp, &
      & 0.0_wp,-0.6278501_wp,-0.5500000_wp, 0.0_wp,-0.1555334_wp, 0.0600000_wp, &
      & 0.0_wp,-0.0338735_wp, 0.3000000_wp, 0.0_wp,-0.2302667_wp,-0.2300000_wp, &
      & 0.0_wp, 0.2494305_wp, 0.0000000_wp, 0.0_wp, 2.2247532_wp,-0.2300000_wp, &
      & 0.0_wp,-0.3000000_wp,-0.4699666_wp, 0.0_wp,-0.3000000_wp,-0.5539659_wp, &
      & 0.0_wp,-0.2769230_wp,-0.5462784_wp, 0.0_wp,-0.2538460_wp,-0.5385909_wp, &
      & 0.0_wp,-0.2307691_wp,-0.5309034_wp, 0.0_wp,-0.2076921_wp,-0.5232158_wp, &
      & 0.0_wp,-0.1846151_wp,-0.5155283_wp, 0.0_wp,-0.1615381_wp,-0.5078408_wp, &
      & 0.0_wp,-0.1384612_wp,-0.5001533_wp, 0.0_wp,-0.1153842_wp,-0.4924658_wp, &
      & 0.0_wp,-0.0923072_wp,-0.4847782_wp, 0.0_wp,-0.0692302_wp,-0.4770907_wp, &
      & 0.0_wp,-0.0461533_wp,-0.4694032_wp, 0.0_wp,-0.0230763_wp,-0.4617157_wp, &
      & 0.0_wp, 0.0000007_wp,-0.4540282_wp, 0.0_wp, 0.1000000_wp,-0.4486165_wp, &
      & 0.0_wp, 0.0500000_wp,-0.3394380_wp, 0.0_wp, 0.3700000_wp,-0.3419199_wp, &
      & 0.0_wp,-0.6000000_wp, 0.6586864_wp, 0.0_wp,-0.6500000_wp, 0.1350223_wp, &
      & 0.0_wp,-0.6500000_wp,-0.0977957_wp, 0.0_wp,-0.6000000_wp,-0.0203212_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0614126_wp, 0.0_wp,-0.5375121_wp, 0.0000000_wp, &
      & 0.0_wp,-0.7133401_wp, 0.0000000_wp, 0.0_wp, 0.7838251_wp, 0.0000000_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0000000_wp, 0.0_wp,-0.8109155_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2532073_wp, 0.2500000_wp, 0.0_wp,-0.0302388_wp,-0.2300000_wp],&
      & shape(shell_hubbard))

   real(wp), parameter :: shell_hubbard_derivs(0:4) = [1.0_wp, 0.5_wp, spread(0.25_wp, 1, 3)]

   real(wp), parameter :: p_hubbard_derivs(max_elem) = 0.1_wp * [&
      & 0.800000_wp, 2.000000_wp, 1.303821_wp, 0.574239_wp, 0.946104_wp, &
      & 1.500000_wp,-0.639780_wp,-0.517134_wp, 1.426212_wp, 0.500000_wp, &
      & 1.798727_wp, 2.349164_wp, 1.400000_wp, 1.936289_wp, 0.711291_wp, &
      &-0.501722_wp, 1.495483_wp,-0.315455_wp, 2.033085_wp, 2.006898_wp, &
      & 0.500000_wp, 1.767268_wp, 0.900000_wp, 0.300000_wp, 0.600000_wp, &
      &-0.500000_wp, 0.300000_wp,-0.200000_wp, 0.500000_wp, 2.312896_wp, &
      & 2.334269_wp,-0.064775_wp, 1.106041_wp, 0.913725_wp, 1.300000_wp, &
      & 0.239815_wp, 2.916203_wp, 1.800000_wp, 0.100000_wp, 0.700000_wp, &
      & 0.500000_wp, 0.919928_wp, 0.600000_wp,-0.500000_wp, 0.300000_wp, &
      & 0.800000_wp, 0.200000_wp, 2.073217_wp, 1.900000_wp,-0.178396_wp, &
      & 1.100000_wp, 0.953683_wp, 1.200000_wp,-0.118925_wp, 2.404185_wp, &
      & 2.069097_wp, 0.012793_wp,-0.100000_wp,-0.100002_wp,-0.100004_wp, &
      &-0.100006_wp,-0.100008_wp,-0.100010_wp,-0.100012_wp,-0.100013_wp, &
      &-0.100015_wp,-0.100017_wp,-0.100019_wp,-0.100021_wp,-0.100023_wp, &
      &-0.100025_wp,-0.100000_wp, 0.200000_wp,-0.200000_wp, 0.800000_wp, &
      & 0.800000_wp,-0.100000_wp, 0.600000_wp, 0.850000_wp,-0.116312_wp, &
      &-0.533933_wp, 0.200000_wp,-0.337508_wp, 1.877978_wp, 1.846485_wp, &
      & 0.097834_wp]

   !> Number of shells
   integer, parameter :: nshell(max_elem) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, &
      & 2, 2, 2, 2, 3, 3]

   !> Angular momentum of each shell
   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
      & 0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 2,  0, 1, 2], shape(ang_shell))

   !> Principal quantum number of each shell
   integer, parameter :: principal_quantum_number(max_shell, max_elem) = reshape([&
      & 1, 0, 0,  1, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  3, 4, 4, &
      & 3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4, &
      & 3, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, &
      & 4, 4, 4,  5, 5, 0,  5, 5, 4,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5, &
      & 4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  5, 5, 0,  5, 5, 5, &
      & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 0,  6, 6, 5, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0, &
      & 6, 6, 5,  6, 6, 5], shape(principal_quantum_number))

   !> Number of primitive gaussians per shell
   integer, parameter :: number_of_primitives(max_shell, max_elem) = reshape([&
      & 3, 0, 0,  3, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0,  4, 4, 0, &
      & 4, 4, 0,  4, 4, 0,  4, 4, 3,  4, 4, 0,  4, 4, 3,  4, 4, 3,  4, 4, 3, &
      & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 0,  4, 4, 3,  3, 4, 4, &
      & 3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4, &
      & 3, 4, 4,  4, 4, 0,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3, &
      & 4, 4, 3,  4, 4, 0,  4, 4, 3,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4, &
      & 3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  4, 4, 0,  4, 4, 3, &
      & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  6, 6, 0,  6, 6, 3, &
      & 3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6, &
      & 3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6, &
      & 3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6,  3, 6, 6, &
      & 3, 6, 6,  3, 6, 6,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0, &
      & 6, 6, 3,  6, 6, 3], shape(number_of_primitives))

   !> Reference occupation of the atom
   real(wp), parameter :: reference_occ(0:2, max_elem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  1.0_wp, 3.0_wp, 0.0_wp, &
      & 1.5_wp, 3.5_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 2.0_wp, &
      & 1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp,  1.0_wp, 1.0_wp, 5.0_wp, &
      & 1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp,  1.0_wp, 1.0_wp, 8.0_wp, &
      & 1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(reference_occ))

   !> Exponent of the Slater function for GFN2-xTB basis set
   real(wp), parameter :: slater_exponent(max_shell, max_elem) = reshape([&
      & 1.230000_wp, 0.000000_wp, 0.000000_wp, 1.669667_wp, 1.500000_wp, 0.000000_wp, &
      & 0.750060_wp, 0.557848_wp, 0.000000_wp, 1.034720_wp, 0.949332_wp, 0.000000_wp, &
      & 1.479444_wp, 1.479805_wp, 0.000000_wp, 2.096432_wp, 1.800000_wp, 0.000000_wp, &
      & 2.339881_wp, 2.014332_wp, 0.000000_wp, 2.439742_wp, 2.137023_wp, 0.000000_wp, &
      & 2.416361_wp, 2.308399_wp, 0.000000_wp, 3.084104_wp, 2.312051_wp, 2.815609_wp, &
      & 0.763787_wp, 0.573553_wp, 0.000000_wp, 1.184203_wp, 0.717769_wp, 1.300000_wp, &
      & 1.352531_wp, 1.391201_wp, 1.000000_wp, 1.773917_wp, 1.718996_wp, 1.250000_wp, &
      & 1.816945_wp, 1.903247_wp, 1.167533_wp, 1.981333_wp, 2.025643_wp, 1.702555_wp, &
      & 2.485265_wp, 2.199650_wp, 2.476089_wp, 2.329679_wp, 2.149419_wp, 1.950531_wp, &
      & 0.875961_wp, 0.631694_wp, 0.000000_wp, 1.267130_wp, 0.786247_wp, 1.380000_wp, &
      & 2.440000_wp, 1.358701_wp, 1.019252_wp, 1.849994_wp, 1.469983_wp, 0.957410_wp, &
      & 1.673577_wp, 1.383176_wp, 0.938025_wp, 1.568211_wp, 1.395427_wp, 1.080270_wp, &
      & 1.839250_wp, 1.222190_wp, 1.240215_wp, 1.911049_wp, 1.022393_wp, 1.294467_wp, &
      & 2.326507_wp, 1.464221_wp, 1.298678_wp, 2.430756_wp, 1.469945_wp, 1.317046_wp, &
      & 2.375425_wp, 1.550837_wp, 1.984703_wp, 1.664847_wp, 1.176434_wp, 0.000000_wp, &
      & 1.720919_wp, 1.591570_wp, 1.050000_wp, 1.990429_wp, 1.830340_wp, 1.100000_wp, &
      & 2.026128_wp, 1.949257_wp, 1.040181_wp, 2.230969_wp, 2.150656_wp, 1.317549_wp, &
      & 2.077587_wp, 2.263120_wp, 1.845038_wp, 2.445680_wp, 2.210494_wp, 1.884991_wp, &
      & 1.017267_wp, 0.870130_wp, 0.000000_wp, 1.419028_wp, 0.928932_wp, 1.500000_wp, &
      & 2.670141_wp, 1.633876_wp, 1.165412_wp, 2.238668_wp, 1.702480_wp, 1.129590_wp, &
      & 1.706832_wp, 1.666463_wp, 1.132172_wp, 1.777658_wp, 1.639917_wp, 1.159781_wp, &
      & 1.918066_wp, 1.918167_wp, 1.346082_wp, 2.102697_wp, 1.749643_wp, 1.348322_wp, &
      & 2.458187_wp, 1.811796_wp, 1.398452_wp, 2.353691_wp, 1.828354_wp, 1.333352_wp, &
      & 2.843549_wp, 1.798462_wp, 1.266649_wp, 1.846689_wp, 1.141823_wp, 0.000000_wp, &
      & 1.963283_wp, 1.685138_wp, 1.050000_wp, 2.551510_wp, 1.893784_wp, 1.100000_wp, &
      & 2.307407_wp, 2.179752_wp, 1.256087_wp, 2.434144_wp, 2.182459_wp, 1.373076_wp, &
      & 2.159500_wp, 2.308379_wp, 1.691185_wp, 2.715140_wp, 2.312510_wp, 1.855707_wp, &
      & 1.225688_wp, 0.823818_wp, 0.000000_wp, 1.528102_wp, 0.991572_wp, 1.500000_wp, &
      & 2.875048_wp, 1.731390_wp, 1.303590_wp, 2.870000_wp, 1.725197_wp, 1.309804_wp, &
      & 2.872308_wp, 1.729767_wp, 1.315495_wp, 2.874615_wp, 1.734337_wp, 1.321186_wp, &
      & 2.876923_wp, 1.738907_wp, 1.326877_wp, 2.879231_wp, 1.743478_wp, 1.332567_wp, &
      & 2.881538_wp, 1.748048_wp, 1.338258_wp, 2.883846_wp, 1.752618_wp, 1.343949_wp, &
      & 2.886154_wp, 1.757188_wp, 1.349640_wp, 2.888462_wp, 1.761758_wp, 1.355331_wp, &
      & 2.890769_wp, 1.766328_wp, 1.361022_wp, 2.893077_wp, 1.770899_wp, 1.366713_wp, &
      & 2.895385_wp, 1.775469_wp, 1.372403_wp, 2.897692_wp, 1.780039_wp, 1.378094_wp, &
      & 2.900000_wp, 1.784609_wp, 1.383785_wp, 2.638729_wp, 2.194333_wp, 1.427467_wp, &
      & 2.018969_wp, 1.996498_wp, 1.407714_wp, 2.155885_wp, 1.892022_wp, 1.458186_wp, &
      & 2.262783_wp, 2.187549_wp, 1.636996_wp, 2.509631_wp, 2.173991_wp, 1.597888_wp, &
      & 2.756134_wp, 2.117548_wp, 1.680343_wp, 2.704492_wp, 2.329136_wp, 1.623286_wp, &
      & 3.241287_wp, 2.183171_wp, 2.084484_wp, 2.244504_wp, 1.470848_wp, 0.000000_wp, &
      & 2.294231_wp, 1.731592_wp, 0.000000_wp, 2.960592_wp, 1.953130_wp, 0.000000_wp, &
      & 2.788267_wp, 2.277039_wp, 0.000000_wp, 3.314810_wp, 2.389456_wp, 0.000000_wp, &
      & 2.220421_wp, 2.408112_wp, 1.500000_wp, 3.109394_wp, 2.541934_wp, 1.790000_wp],&
      & shape(slater_exponent))

   real(wp), parameter :: mp_dmp3 = 3.0_wp, mp_dmp5 = 4.0_wp
   real(wp), parameter :: mp_shift = 1.2_wp, mp_kexp = 4.0_wp, mp_rmax = 5.0_wp
   !> Dipole exchange-correlation kernel
   real(wp), parameter :: p_dkernel(*) = 0.01_wp * [&
      & 5.563889_wp,-1.000000_wp,-0.500000_wp,-0.613341_wp,-0.481186_wp, &
      &-0.411674_wp, 3.521273_wp,-4.935670_wp,-8.339183_wp,10.000000_wp, &
      & 0.000000_wp,-0.082005_wp, 2.633341_wp,-0.025750_wp, 2.110225_wp, &
      &-0.151117_wp,-2.536958_wp,-2.077329_wp,-0.103383_wp,-0.236675_wp, &
      &-0.515177_wp,-0.434506_wp,-0.350000_wp, 0.149669_wp,-0.759168_wp, &
      & 0.412929_wp,-0.247938_wp,-1.261887_wp,-0.700000_wp,-0.100000_wp, &
      & 0.267219_wp, 0.108460_wp,-0.201294_wp,-0.288648_wp,-1.088586_wp, &
      &-0.889357_wp,-0.093328_wp,-0.459925_wp,-0.637291_wp,-0.599615_wp, &
      &-0.288729_wp, 0.346327_wp,-0.458416_wp,-0.081922_wp, 0.007016_wp, &
      &-0.310361_wp,-0.800314_wp,-0.105364_wp, 0.951079_wp, 0.085029_wp, &
      &-0.015519_wp,-0.263414_wp,-0.603648_wp,-0.214447_wp,-0.080000_wp, &
      &-0.260000_wp,-0.395198_wp,-0.723806_wp,-0.704819_wp,-0.685832_wp, &
      &-0.666845_wp,-0.647858_wp,-0.628871_wp,-0.609884_wp,-0.590897_wp, &
      &-0.571910_wp,-0.552923_wp,-0.533936_wp,-0.514949_wp,-0.495961_wp, &
      &-0.476974_wp,-0.537685_wp,-0.200343_wp, 0.065886_wp,-0.587636_wp, &
      &-0.510090_wp,-0.673822_wp,-0.423684_wp, 0.393418_wp,-0.250000_wp, &
      & 0.374018_wp, 1.007016_wp,-0.737252_wp,-1.344854_wp,-0.348123_wp, &
      &-0.167597_wp]
   !> Quadrupole exchange-correlation kernel
   real(wp), parameter :: p_qkernel(*) = 0.01_wp * [&
      & 0.027431_wp,-0.337528_wp, 0.020000_wp,-0.058586_wp,-0.058228_wp, &
      & 0.213583_wp, 2.026786_wp,-0.310828_wp,-0.245955_wp,-0.500000_wp, &
      & 0.020000_wp,-0.005516_wp,-0.021887_wp,-0.080000_wp, 0.028679_wp, &
      & 0.442859_wp, 0.122783_wp,-1.083404_wp, 0.025000_wp, 0.010000_wp, &
      &-0.042004_wp, 0.059660_wp, 0.009764_wp, 0.137744_wp, 0.229903_wp, &
      & 0.267734_wp, 0.048237_wp,-0.080000_wp,-0.345631_wp, 0.007658_wp, &
      &-0.003616_wp,-0.003589_wp, 0.014149_wp, 0.085728_wp, 0.216935_wp, &
      &-0.415024_wp, 0.015000_wp, 0.015000_wp, 0.010460_wp,-0.012944_wp, &
      & 0.041491_wp, 0.312549_wp, 0.155242_wp, 0.359228_wp, 0.008570_wp, &
      &-0.040485_wp,-0.020810_wp, 0.012250_wp,-0.002031_wp,-0.008243_wp, &
      &-0.020630_wp,-0.026864_wp, 0.069660_wp,-0.156200_wp, 0.008000_wp, &
      & 0.015000_wp,-0.030000_wp,-0.025000_wp,-0.024615_wp,-0.024231_wp, &
      &-0.023846_wp,-0.023462_wp,-0.023077_wp,-0.022692_wp,-0.022308_wp, &
      &-0.021923_wp,-0.021538_wp,-0.021154_wp,-0.020769_wp,-0.020385_wp, &
      &-0.020000_wp,-0.016478_wp, 0.039599_wp, 1.063309_wp, 0.306870_wp, &
      & 0.759049_wp, 0.322935_wp, 0.098019_wp,-0.020320_wp,-0.032901_wp, &
      &-0.008506_wp,-0.001670_wp, 0.162529_wp, 0.013818_wp, 0.021624_wp, &
      &-0.111556_wp]
   !> Valence coordination number for radii
   real(wp), parameter :: p_vcn(*) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, &
      & 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp]
   !> Cutoff radii for multipole electrostatics
   real(wp), parameter :: p_rad(*) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 4.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp]


   !> Specification of the
   type, public, extends(tb_h0spec) :: gfn2_h0spec
      real(wp) :: kshell(0:2, 0:2)
      real(wp), allocatable :: kpair(:, :)
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure :: get_hscale
      !> Generator for the polynomial parameters for the distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the reference occupation numbers of the atoms
      procedure :: get_reference_occ
   end type gfn2_h0spec

   interface gfn2_h0spec
      module procedure :: new_gfn2_h0spec
   end interface gfn2_h0spec

contains


subroutine new_gfn2_calculator(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call add_basis(calc, mol)
   call add_ncoord(calc, mol)
   call add_hamiltonian(calc, mol)
   call add_repulsion(calc, mol)
   call add_dispersion(calc, mol)
   call add_coulomb(calc, mol)

end subroutine new_gfn2_calculator

subroutine add_basis(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: isp, izp, ish, stat, ng
   integer, allocatable :: nsh_id(:)
   type(cgto_type), allocatable :: cgto(:, :)

   nsh_id = nshell(mol%num)
   allocate(cgto(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nsh_id(isp)
         ng = number_of_primitives(ish, izp)
         call slater_to_gauss(ng, principal_quantum_number(ish, izp), ang_shell(ish, izp), &
            & slater_exponent(ish, izp), cgto(ish, isp), .true., stat)
      end do
   end do

   call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

end subroutine add_basis

subroutine add_ncoord(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_ncoord(calc%ncoord, mol, cn_type="gfn")
end subroutine add_ncoord

subroutine add_hamiltonian(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_hamiltonian(calc%h0, mol, calc%bas, new_gfn2_h0spec(mol))
end subroutine add_hamiltonian

subroutine add_dispersion(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   type(d4_dispersion), allocatable :: tmp,tmp2


   allocate(tmp)
   call new_d4_dispersion(tmp, mol, s6=s6, s8=s8, a1=a1, a2=a2, s9=s9)
   call move_alloc(tmp, calc%dispersion)
   allocate(tmp2)
   call new_d4_dispersion(tmp2, mol, s6=0.0_wp, s8=0.0_wp, a1=a1, a2=a2, s9=s9)
   call move_alloc(tmp2, calc%dispersion_3body)
   
end subroutine add_dispersion

subroutine add_repulsion(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: alpha(:), zeff(:)

   allocate(calc%repulsion)
   alpha = rep_alpha(mol%num)
   zeff = rep_zeff(mol%num)
   call new_repulsion(calc%repulsion, mol, alpha, zeff, rep_kexp, rep_kexp_light, rep_rexp)
end subroutine add_repulsion

subroutine add_coulomb(calc, mol)
   !> Instance of the xTB evaluator
   type(sqmbox_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: hardness(:, :), hubbard_derivs(:, :)
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)
   type(effective_coulomb), allocatable :: es2

   allocate(calc%coulomb)
   allocate(es2)
   call get_shell_hardness(mol, calc%bas, hardness)
   call new_effective_coulomb(es2, mol, gexp, hardness, arithmetic_average, &
      & calc%bas%nsh_id)
   call move_alloc(es2, calc%coulomb%es2)

   allocate(calc%coulomb%es3)
   call get_hubbard_derivs(mol, calc%bas, hubbard_derivs)
   call new_onsite_thirdorder(calc%coulomb%es3, mol, hubbard_derivs, calc%bas%nsh_id)

   allocate(calc%coulomb%aes2)
   dkernel = p_dkernel(mol%num)
   qkernel = p_qkernel(mol%num)
   rad = p_rad(mol%num)
   vcn = p_vcn(mol%num)
   call new_damped_multipole(calc%coulomb%aes2, mol, mp_dmp3, mp_dmp5, dkernel, qkernel, &
      & mp_shift, mp_kexp, mp_rmax, rad, vcn)

end subroutine add_coulomb

subroutine get_shell_hardness(mol, bas, hardness)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell resolved hardness parameters
   real(wp), allocatable, intent(out) :: hardness(:, :)

   integer :: isp, izp, ish, il

   allocate(hardness(maxval(bas%nsh_id), mol%nid))
   hardness(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hardness(ish, isp) = hubbard_parameter(izp) * shell_hubbard(il, izp)
      end do
   end do
end subroutine get_shell_hardness


subroutine get_hubbard_derivs(mol, bas, hubbard_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell resolved Hubbard derivatives
   real(wp), allocatable, intent(out) :: hubbard_derivs(:, :)

   integer :: isp, izp, ish, il

   allocate(hubbard_derivs(maxval(bas%nsh_id), mol%nid))
   hubbard_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hubbard_derivs(ish, isp) = p_hubbard_derivs(izp) * shell_hubbard_derivs(il)
      end do
   end do
end subroutine get_hubbard_derivs




elemental function kshell(k, l)
   integer, intent(in) :: k, l
   real(wp) :: kshell
   kshell = merge(2.0_wp, (kdiag(l)+kdiag(k))/2, k==2.and.any(l==[0,1]).or.l==2.and.any(k==[0,1]))
end function kshell

end module tblite_xtb_gfn2
