#
# Copyright (C) 1997 Arve Kylling
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 1, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# To obtain a copy of the GNU General Public License write to the
# Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
# usa.

#-------------------------------------------------------------------------

use lib "../src";
use Carp;
use Cwd;
use File::Copy;
use File::Spec;
use Getopt::Long;
use Time::HiRes qw (time);
use UVTools;
use strict;

my $found_libnetcdf      = "yes"      eq "yes";
my $found_libgsl         = "yes"         eq "yes";
my $found_avhrr          = "yes"          eq "yes";
my $found_mystic         = "yes"         eq "yes";
my $found_mystic3d       = "no"       eq "yes";
my $found_twomaxrnd3c    = "yes"    eq "yes";
my $found_tenstream      = "no"      eq "yes";
my $found_dts            = "no"            eq "yes"; 
my $found_vroom          = "yes"          eq "yes";
my $found_lidar          = "no"          eq "yes";
my $found_lowtran        = "yes"        eq "yes";
my $found_sos            = "yes"            eq "yes";
my $found_opac           = "no"           eq "yes";
my $found_wcmie          = "no"          eq "yes";
my $found_baum           = "no"           eq "yes";
my $found_baum_v36       = "no"       eq "yes";
my $found_hey            = "no"            eq "yes";
my $found_reptran_coarse = "yes" eq "yes";
my $found_reptran_ext    = "yes"    eq "yes";
my $found_hamster        = "no"        eq "yes"; 
my $found_tests          = "yes"          eq "yes";
my $found_nco            = "no"            eq "yes";
my $found_nccmp          = "no"          eq "yes";

my $global_t_start = time;

my $help                  = 0;
my $verbose               = 0;
my $verbose_uvspec        = 0;
my $dryrun                = 0;
my $answers               = 0;
my $force_answers         = 0;
my $valgrind              = 0;
my $all                   = 1;   # =0 just requested test is done
                                 # =1 all except airmass, fdisort and ftwostr tests are done
                                 # =2 all except airmass tests are done
                                 # =3 all tests including airmass test are done
                                 #    invoked with --test="everything"
                                 # =4 all tests including all optional packages tests

my $mc_photons_multiplier = 1;   # multiply the number of photons in montecarlo tests to increase quality of results
                                 # is set to 100 in case of answers to get better baseline results

my $test                  = "default";

my $global_atol = -1; # Lower limit below which a result is neglected
my $global_rtol = -1; # Accepted relative error between test cases and answers

my $ret = GetOptions(
         "-h"                       , \$help,
         "--help"                   , \$help,
         "--verbose"                , \$verbose,
         "--verbose_uvspec"         , \$verbose_uvspec,
         "--dry"                    , \$dryrun,
         "--answers"                , \$answers,
         "--force_answers"          , \$force_answers,
         "--global_atol=f"          , \$global_atol,
         "--global_rtol=f"          , \$global_rtol,
         "--valgrind"               , \$valgrind,
         "--mc_photons_multiplier=f", \$mc_photons_multiplier,
         "--test=s"                 , \$test);

die "Could not read commandline options" if (!$ret);


if ( $test ne "default" )    { $all = 0; }
if ( $test eq "everything" ) { $all = 3; }

my $fail = 0;       # number of tests where uvspec returned with an error code
my $fail_diff = 0;  # number of tests with differences in output

my $ITER_MC_TESTS = 3; # repeat tests if we found differences in outputs, if none of them produces correct results, we assume it a failed test
if ($answers) { $mc_photons_multiplier = $mc_photons_multiplier * 100; }

my $show_help = $help || $verbose;
if ($show_help) {
  printf stderr "Running various libRadtran tests. This may take some time....\n";
  printf stderr "\n";
  printf stderr "The numbers in the parenthesis behind the name of the tests are:\n";
  printf stderr " 1st number:  The lower absolute limit of values included in the test. \n";
  printf stderr "              Values in the output less than limit are ignored. \n";
  printf stderr " 2nd number:  The maximum relative difference allowed between local test\n";
  printf stderr "              results and the standard results (in percentage).\n";
  printf stderr "\n";
  printf stderr "If this is still unclear, check the source in test/test.pl.in.\n";
  printf stderr "\n";
  printf stderr "\n";
  printf stderr "Usage:\n";
  printf stderr "  --verbose                     : More verbose output ($verbose)\n";
  printf stderr "  --answers                     : Create expected output files ($answers)\n";
  printf stderr "  --valgrind                    : Use valgrind to call tests ($valgrind)\n";
  printf stderr "  --mc_photons_multiplier=<>    : Multiply used photons for MonteCarlo tests ($mc_photons_multiplier)\n";
  printf stderr "  --test=<>                     : Provide regular expression for test case names ($test)\n";
  printf stderr "  --global_atol=<>              : Use specified absolute tolerance for all tests ($global_atol)\n";
  printf stderr "  --global_rtol=<>              : Use specified relative tolerance for all tests ($global_rtol)\n";
  printf stderr "\n";
}
exit if $help;

my $inp = "TMP.INP";
my $out = "TMP.OUT";

#
# basic form of test entries is:
# [ TestName,  (uvspec inp file , outfile if different), Description, atol, rtol, testlvl (as in var $all), constraints for the test, i.e. if it should be possible to run it ]
#
my @standardTests= (
 ["uvspec_aerosol_moments"               ,["UVSPEC_AEROSOL_MOMENTS"                     ,                         ],"disort aerosol moments"                       ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_aerosol_opac"                  ,["UVSPEC_AEROSOL_OPAC"                        ,                         ],"disort aerosol opac"                          ,0.00001,0.001  ,1,$found_opac                                     ],
 ["uvspec_aerosol_refrac"                ,["UVSPEC_AEROSOL_REFRAC"                      ,                         ],"disort aerosol refractive index"              ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_aerosol"                       ,["UVSPEC_AEROSOL"                             ,                         ],"disort aerosol"                               ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_avhrr_solar_ch1_fdisort2"      ,["tests/UVSPEC_AVHRR_SOLAR_CH1_FDISORT2"      ,                         ],"AVHRR (Kratz 1995 )channel 1 solar fdisort2"  ,0.00001,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_solar_ch1"               ,["UVSPEC_AVHRR_SOLAR_CH1"                     ,                         ],"AVHRR (Kratz 1995 )channel 1 solar"           ,0.00001,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_avhrr_solar_ch2_fdisort2"      ,["tests/UVSPEC_AVHRR_SOLAR_CH2_FDISORT2"      ,                         ],"AVHRR (Kratz 1995 )channel 2 solar fdisort2"  ,0.001  ,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_solar_ch2"               ,["UVSPEC_AVHRR_SOLAR_CH2"                     ,                         ],"AVHRR (Kratz 1995 )channel 2 solar"           ,0.001  ,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_avhrr_solar_ch3_fdisort2"      ,["tests/UVSPEC_AVHRR_SOLAR_CH3_FDISORT2"      ,                         ],"AVHRR (Kratz 1995 )channel 3 solar fdisort2"  ,0.00001,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_solar_ch3"               ,["UVSPEC_AVHRR_SOLAR_CH3"                     ,                         ],"AVHRR (Kratz 1995 )channel 3 solar"           ,0.00001,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_avhrr_thermal_ch3_fdisort2"    ,["tests/UVSPEC_AVHRR_THERMAL_CH3_FDISORT2"    ,                         ],"AVHRR (Kratz 1995 )channel 3 thermal fdisort2",0.00001,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_thermal_ch3"             ,["UVSPEC_AVHRR_THERMAL_CH3"                   ,                         ],"AVHRR (Kratz 1995 )channel 3 thermal"         ,0.00001,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_avhrr_thermal_ch4_fdisort2"    ,["tests/UVSPEC_AVHRR_THERMAL_CH4_FDISORT2"    ,                         ],"AVHRR (Kratz 1995 )channel 4 thermal fdisort2",0.00001,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_thermal_ch4"             ,["UVSPEC_AVHRR_THERMAL_CH4"                   ,                         ],"AVHRR (Kratz 1995 )channel 4 thermal"         ,0.00001,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_avhrr_thermal_ch5_fdisort2"    ,["tests/UVSPEC_AVHRR_THERMAL_CH5_FDISORT2"    ,                         ],"AVHRR (Kratz 1995 )channel 5 thermal fdisort2",0.00001,0.001  ,2,$found_avhrr && $found_tests                    ],
 ["uvspec_avhrr_thermal_ch5"             ,["UVSPEC_AVHRR_THERMAL_CH5"                   ,                         ],"AVHRR (Kratz 1995 )channel 5 thermal"         ,0.00001,0.001  ,1,$found_avhrr                                    ],
 ["uvspec_clear"                         ,["UVSPEC_CLEAR"                               ,                         ],"disort clear sky"                             ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_cloudcover_ck_redistribution"  ,["UVSPEC_CLOUDCOVER_REDISTRIBUTION"           ,"UVSPEC_CLOUDCOVER"      ],"cloudcover correlated_k and redistribution"   ,0.001  ,0.018  ,1,"no_constraint"                                 ],
 ["uvspec_cloudcover_ck"                 ,["UVSPEC_CLOUDCOVER"                          ,                         ],"cloudcover and correlated_k"                  ,0.001  ,0.018  ,1,"no_constraint"                                 ],
 ["uvspec_c_twostr"                      ,["UVSPEC_TWOSTR"                              ,                         ],"c_twostr aerosol and water cloud"             ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_disort"                        ,["UVSPEC_DISORT"                              ,                         ],"disort aerosol"                               ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_fdisort2"                      ,["UVSPEC_FDISORT2"                            ,                         ],"fdisort2 aerosol"                             ,0.00001,0.001  ,2,"no_constraint"                                 ],
 ["uvspec_fluorescence"                  ,["UVSPEC_FLUORESCENCE"                        ,                         ],"fluorescence"                                 ,0.1    ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_f_twostr"                      ,["UVSPEC_FTWOSTR"                             ,                         ],"twostr aerosol and water cloud"               ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_fuliou_ic"                     ,["UVSPEC_FULIOU_IC"                           ,                         ],"correlated-k (Fu and Liou 1992/93 )ice clouds",0.01   ,0.015  ,1,"no_constraint"                                 ],
 ["uvspec_fuliou_thermal"                ,["UVSPEC_FULIOU_THERMAL"                      ,                         ],"Fu and Liou thermal irradiance disort"        ,0.00001,0.009  ,1,"no_constraint"                                 ],
 ["uvspec_fuliou"                        ,["UVSPEC_FULIOU"                              ,                         ],"correlated-k (Fu and Liou 1992/93 )disort"    ,0.01   ,0.015  ,1,"no_constraint"                                 ],
 ["uvspec_himalaja_fdisort2"             ,["tests/UVSPEC_HIMALAJA_FDISORT2"             ,                         ],"fdisort2 albedo and altitude map"             ,0.00001,0.007  ,2,$found_libnetcdf && $found_tests                ],
 ["uvspec_himalaja"                      ,["UVSPEC_HIMALAJA"                            ,                         ],"disort albedo and altitude map"               ,0.00001,0.007  ,1,"no_constraint"                                 ],
 ["uvspec_ic_baum_fdisort2"              ,["tests/UVSPEC_IC_BAUM_FDISORT2"              ,                         ],"fdisort2 ice cloud Baum et al. (2005)"        ,0.001  ,0.001  ,1,$found_baum && $found_libnetcdf && $found_tests ],
 ["uvspec_ic_baum"                       ,["UVSPEC_IC_BAUM"                             ,                         ],"disort ice cloud Baum et al. (2005)"          ,0.00001,0.001  ,1,$found_baum                                     ],
 ["uvspec_ic_wc_fdisort2"                ,["tests/UVSPEC_IC_WC_FDISORT2"                ,                         ],"fdisort2 water and ice clouds"                ,0.00001,0.002  ,1,$found_tests                                    ],
 ["uvspec_ic_wc"                         ,["UVSPEC_IC_WC"                               ,                         ],"disort water and ice clouds"                  ,0.00001,0.002  ,1,"no_constraint"                                 ],
 ["uvspec_ipa_ck"                        ,["UVSPEC_IPA_CK"                              ,                         ],"IPA and correlated_k"                         ,0.001  ,0.002  ,1,"no_constraint"                                 ],
 ["uvspec_kato2"                         ,["UVSPEC_KATO2"                               ,                         ],"correlated-k new Kato tables twostr"          ,0.001  ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_kato"                          ,["UVSPEC_KATO"                                ,                         ],"correlated-k (Kato et al.1999 )twostr"        ,0.001  ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_mystic_1dabs_thermal"          ,["UVSPEC_MC_ABS1D_THERMAL"                    ,                         ],"MYSTIC 3D abso thermal 1D benchmark (disort)" ,0.00001,0.001  ,1,$found_mystic3d                                 ],
 ["uvspec_mystic_1dabs"                  ,["UVSPEC_MC_ABS1D"                            ,                         ],"MYSTIC 3D abso 1D benchmark (disort)"         ,0.00001,0.001  ,1,$found_mystic3d                                 ],
 ["uvspec_mystic_3dabs_thermal"          ,["UVSPEC_MC_ABS3D_THERMAL"                    ,"UVSPEC_MC_ABS1D_THERMAL"],"MYSTIC 3D absorption thermal"                 ,0.01   ,0.1    ,1,$found_mystic3d                                 ],#
 ["uvspec_mystic_3dabs"                  ,["UVSPEC_MC_ABS3D"                            ,"UVSPEC_MC_ABS1D"        ],"MYSTIC 3D absorption"                         ,0.01   ,0.1    ,1,$found_mystic3d                                 ],#
 ["uvspec_mystic_albedo2d_spectral_file" ,["UVSPEC_MC_ALBEDO2D_SPECTRAL"                ,                         ],"MYSTIC 2D spectral albedo"                    ,0.0005 ,0.01   ,1,$found_mystic3d                                 ],
 ["uvspec_mystic_raytracing"             ,["UVSPEC_MC_RAYTRACING"                       ,                         ],"MYSTIC ice particle raytracing"               ,0.00001,0.02   ,1,$found_mystic3d && $found_opac                  ],
 ["uvspec_mystic"                        ,["UVSPEC_MC"                                  ,                         ],"MYSTIC water cloud"                           ,0.00001,0.03   ,1,$found_mystic3d                                 ],
 ["uvspec_overlap_max_fdisort2"          ,["tests/UVSPEC_OVERLAP_MAX_FDISORT2"          ,                         ],"fdisort2 cloud overlap maximum"               ,0.00001,0.002  ,2,$found_tests                                    ],
 ["uvspec_overlap_maxrandom_fdisort2"    ,["tests/UVSPEC_OVERLAP_MAXRAND_FDISORT2"      ,                         ],"fdisort2 cloud overlap maximum-random"        ,0.00001,0.002  ,1,$found_tests                                    ],
 ["uvspec_overlap_maxrandom"             ,["UVSPEC_OVERLAP_MAXRAND"                     ,                         ],"disort cloud overlap maximum-random"          ,0.00001,0.002  ,1,"no_constraint"                                 ],
 ["uvspec_overlap_max"                   ,["UVSPEC_OVERLAP_MAX"                         ,                         ],"disort cloud overlap maximum"                 ,0.00001,0.002  ,1,"no_constraint"                                 ],
 ["uvspec_overlap_random_fdisort2"       ,["tests/UVSPEC_OVERLAP_RAND_FDISORT2"         ,                         ],"fdisort2 cloud overlap random"                ,0.00001,0.002  ,1,$found_tests                                    ],
 ["uvspec_overlap_random"                ,["UVSPEC_OVERLAP_RAND"                        ,                         ],"disort cloud overlap random"                  ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_profiles1_fdisort2"            ,["tests/UVSPEC_PROFILES1_FDISORT2"            ,                         ],"profiles 1 fdisort2"                          ,0.00001,0.002  ,2,$found_tests                                    ],
 ["uvspec_profiles1"                     ,["UVSPEC_PROFILES1"                           ,                         ],"profiles 1"                                   ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_profiles2_fdisort2"            ,["tests/UVSPEC_PROFILES2_FDISORT2"            ,                         ],"profiles 2 fdisort2"                          ,0.00001,0.001  ,2,$found_tests                                    ],
 ["uvspec_profiles2"                     ,["UVSPEC_PROFILES2"                           ,                         ],"profiles 2"                                   ,0.00001,0.002  ,1,"no_constraint"                                 ],
 ["uvspec_profiles3_fdisort2"            ,["tests/UVSPEC_PROFILES3_FDISORT2"            ,                         ],"profiles 3 fdisort2"                          ,0.00001,0.001  ,2,$found_tests                                    ],
 ["uvspec_profiles3"                     ,["UVSPEC_PROFILES3"                           ,                         ],"profiles 3"                                   ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_profiles4_fdisort2"            ,["tests/UVSPEC_PROFILES4_FDISORT2"            ,                         ],"profiles 4 fdisort2"                          ,0.00001,0.001  ,2,$found_tests                                    ],
 ["uvspec_profiles4"                     ,["UVSPEC_PROFILES4"                           ,                         ],"profiles 4"                                   ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_radiances_fdisort2"            ,["tests/UVSPEC_RADIANCES_FDISORT2"            ,                         ],"fdisort2 radiances"                           ,0.00001,0.003  ,2,$found_tests                                    ],
 ["uvspec_radiances"                     ,["UVSPEC_RADIANCES"                           ,                         ],"disort radiances"                             ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_radiosonde"                    ,["UVSPEC_RADIOSONDE"                          ,                         ],"radiosonde"                                   ,0.00001,0.004  ,1,$found_avhrr                                    ],
 ["uvspec_raman"                         ,["UVSPEC_RAMAN"                               ,                         ],"raman"                                        ,0.1    ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_reptran_channel_thermal"       ,["UVSPEC_REPTRAN_CHANNEL_THERMAL"             ,                         ],"reptran_channel_thermal"                      ,0.01   ,0.001  ,1,$found_reptran_ext                              ],
 ["uvspec_reptran_solar"                 ,["UVSPEC_REPTRAN_SOLAR"                       ,                         ],"reptran_solar"                                ,0.01   ,0.002  ,1,$found_reptran_coarse                           ],
    ["uvspec_reptran_thermal"               ,["UVSPEC_REPTRAN_THERMAL"                     ,                      ],"reptran_thermal"                              ,0.01   ,0.001  ,1,$found_reptran_coarse                           ],
#  ["uvspec_repint_solar"                 ,["UVSPEC_REPINT_SOLAR"                       ,                          ],"repint_solar"                                ,0.01   ,0.002  ,1, "no_constraint"                         ],
 ["uvspec_repint_thermal"               ,["UVSPEC_REPINT_THERMAL"                     ,                           ],"repint_thermal"                              ,0.01   ,0.001  ,1, "no_constraint"                           ],
 ["uvspec_rodents"                       ,["UVSPEC_RODENTS"                             ,                         ],"rodents aerosol and water cloud"              ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_rgb_disort"                    ,["tests/UVSPEC_RGB"                           ,                         ],"rgb values disort"                            ,0.00001,0.001  ,1,$found_tests                                    ],
 ["uvspec_rodents_zout_thermal"          ,["UVSPEC_RODENTS_ZOUT_THERMAL"                ,                         ],"rodents aerosol and water cloud zout thermal" ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_rodents_zout"                  ,["UVSPEC_RODENTS_ZOUT"                        ,                         ],"rodents aerosol and water cloud zout"         ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_sbdart_ftwostr"                ,["tests/UVSPEC_LOWTRAN_SOLAR_FTWOSTR"         ,                         ],"LOWTRAN absorption parameterization"          ,0.00001,0.022  ,2,$found_tests && $found_lowtran                  ],
 ["uvspec_sbdart"                        ,["UVSPEC_LOWTRAN_SOLAR"                       ,                         ],"LOWTRAN absorption parameterization"          ,0.00001,0.022  ,1,$found_lowtran                                  ],
 ["uvspec_scia_approx_fdisort2"          ,["tests/UVSPEC_SCIA_APPROX_FDISORT2"          ,                         ],"fdisort2 SCIAMACHY HG approximation"          ,0.00001,0.003  ,2,$found_tests                                    ],
 ["uvspec_scia_approx"                   ,["UVSPEC_SCIA_APPROX"                         ,                         ],"disort SCIAMACHY HG approximation"            ,0.00001,0.003  ,1,"no_constraint"                                 ],
 ["uvspec_scia_fdisort2"                 ,["tests/UVSPEC_SCIA_FDISORT2"                 ,                         ],"fdisort2 SCIAMACHY simulation"                ,0.00001,0.003  ,2,$found_wcmie && $found_tests                    ],
 ["uvspec_scia"                          ,["UVSPEC_SCIA"                                ,                         ],"disort SCIAMACHY simulation"                  ,0.00001,0.003  ,1,$found_wcmie                                    ],
 ["uvspec_sdisort"                       ,["UVSPEC_SDISORT"                             ,                         ],"sdisort aerosol"                              ,0.00001,0.001  ,2,"no_constraint"                                 ],
 ["uvspec_simple"                        ,["UVSPEC_SIMPLE"                              ,                         ],"uvspec simple"                                ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_so2"                           ,["UVSPEC_SO2"                                 ,                         ],"disort SO2"                                   ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_spherical_albedo"              ,["UVSPEC_SPHER_ALB"                           ,                         ],"disort_spherical_albedo"                      ,0.01   ,0.001  ,1,$found_reptran_coarse                           ],
 ["uvspec_sslidar_pol"                   ,["UVSPEC_SSLIDAR_POL"                         ,                         ],"single scattering polarized lidar water cloud",0.00001,0.001  ,1,$found_wcmie                                    ],
 ["uvspec_sslidar"                       ,["UVSPEC_SSLIDAR"                             ,                         ],"single scattering lidar water cloud"          ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_test_brdf_ambrals_file"        ,["UVSPEC_BRDF_AMBRALS_FILE"                   ,                         ],"disort BRDF AMBRALS file"                     ,0.0005 ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_test_brdf_ambrals"             ,["UVSPEC_BRDF_AMBRALS"                        ,                         ],"disort BRDF AMBRALS"                          ,0.0005 ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_test_brdf_fdisort2"            ,["tests/UVSPEC_BRDF_FDISORT2"                 ,                         ],"fdisort2 BRDF"                                ,0.0005 ,0.001  ,1,$found_tests                                    ],
 ["uvspec_test_brdf_hapke"               ,["UVSPEC_BRDF_HAPKE"                          ,                         ],"disort BRDF Hapke"                            ,0.0005 ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_test_brdf_rossli"              ,["UVSPEC_BRDF_ROSSLI"                         ,                         ],"disort BRDF Ross-Li"                          ,0.0005 ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_test_brdf"                     ,["UVSPEC_BRDF"                                ,                         ],"disort BRDF"                                  ,0.0005 ,0.001  ,1,"no_constraint"                                 ],
#["uvspec_test_wc_ipa"                   ,["UVSPEC_WC_IPA"                              ,                         ],"wc_ipa"                                       ,0.00001,0.001  ,1,"no_constraint"                                 ],#broken
 ["uvspec_transmittance_wl_file_fdisort2",["tests/UVSPEC_TRANSMITTANCE_WL_FILE_FDISORT2",                         ],"transmittance_wl_file fdisort2"               ,0.001  ,0.001  ,2,$found_libnetcdf && $found_tests                ],
 ["uvspec_transmittance_wl_file"         ,["UVSPEC_TRANSMITTANCE_WL_FILE"               ,                         ],"transmittance_wl_file"                        ,0.001  ,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twomaxrnd3c"                   ,["UVSPEC_TWOMAXRND3C"                         ,                         ],"twomaxrand triple cloud"                      ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twomaxrnd3cwcic"               ,["UVSPEC_TWOMAXRND3CWCIC"                     ,                         ],"twomaxrand triple wc/ic cloud"                ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twomaxrnd"                     ,["UVSPEC_TWOMAXRND"                           ,                         ],"twomaxrand water cloud"                       ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twomaxrndwcic"                 ,["UVSPEC_TWOMAXRNDWCIC"                       ,                         ],"twomaxrand wc/ic cloud"                       ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twostrebe"                     ,["UVSPEC_TWOSTREBE"                           ,                         ],"twostrebe aerosol and water cloud"            ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_twostrpp"                      ,["UVSPEC_TWOSTRPP"                            ,                         ],"twostrpp aerosol and water cloud"             ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_tzs_channel"                   ,["UVSPEC_TZS_CHANNEL"                         ,                         ],"tzs_channel"                                  ,0.01   ,0.001  ,1,$found_reptran_ext                              ],
 ["uvspec_tzs"                           ,["UVSPEC_TZS"                                 ,                         ],"uvspec_tzs"                                   ,0.01   ,0.001  ,1,$found_reptran_coarse                           ],
 ["uvspec_wc_ic_ipa_files"               ,["UVSPEC_WC_IC_IPA_FILES"                     ,                         ],"wc_ipa_files ic_ipa_files and redistribution" ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_wc_ipa_cooling"                ,["UVSPEC_COOLING_IPA"                         ,                         ],"cooling rates and wc_ipa_files"               ,0.001  ,0.019  ,1,"no_constraint"                                 ],
 ["uvspec_wc_ipa_files_ck_redistribution",["UVSPEC_WC_IPA_FILES_REDISTRIBUTION"         ,"UVSPEC_WC_IPA_FILES"    ],"wc_ipa_files correlated_k and redistribution" ,0.001  ,0.00135,1,"no_constraint"                                 ],
 ["uvspec_wc_ipa_files_ck"               ,["UVSPEC_WC_IPA_FILES"                        ,                         ],"wc_ipa_files and correlated_k"                ,0.001  ,0.00135,1,"no_constraint"                                 ],
 ["uvspec_wc_ipa_heating"                ,["UVSPEC_HEATING_IPA"                         ,                         ],"heating rates and wc_ipa_files"               ,0.001  ,0.0017 ,1,"no_constraint"                                 ],
 ["uvspec_wc_mie_fdisort2"               ,["tests/UVSPEC_WC_MIE_FDISORT2"               ,                         ],"fdisort2 water cloud (Mie)"                   ,0.001  ,0.01   ,2,$found_wcmie && $found_libnetcdf && $found_tests],
 ["uvspec_wc_mie"                        ,["UVSPEC_WC_MIE"                              ,                         ],"disort water cloud (Mie)"                     ,0.00001,0.001  ,1,$found_wcmie                                    ],
 ["uvspec_wc_moments_fdisort2"           ,["tests/UVSPEC_WC_MOMENTS_FDISORT2"           ,                         ],"fdisort2 wc Legendre moments"                 ,0.00001,0.001  ,2,$found_tests                                    ],
 ["uvspec_wc_moments"                    ,["UVSPEC_WC_MOMENTS"                          ,                         ],"disort wc Legendre moments"                   ,0.00001,0.001  ,1,"no_constraint"                                 ],
 ["uvspec_wc"                            ,["UVSPEC_WC"                                  ,                         ],"disort water cloud"                           ,0.00001,0.001  ,1,"no_constraint"                                 ],
);
foreach my $test (@standardTests) { run_test_entry(\&uvspec_run_test, $test); }

my @mc_thermal_forward_tests = (
 ["mc_thermal_forward_disort_ckfu"        ,["mc_thermal_forward/DISORT_THERMAL_CKFU"       ,]                                               ,"mc_thermal_forward_disort_ckfu"        ,0.001  ,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_ckfu"            ,["mc_thermal_forward/MC_THERMAL_CKFU"           ,"mc_thermal_forward/DISORT_THERMAL_CKFU"]       ,"mc_thermal_forward_mc_ckfu"            ,0.1    ,0.01 ,1,"no_constraint"],
 ["mc_thermal_forward_disort_ckfu_hr"     ,["mc_thermal_forward/DISORT_THERMAL_CKFU_HR"    ,]                                               ,"mc_thermal_forward_disort_ckfu_hr"     ,0.001  ,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_ckfu_hr"         ,["mc_thermal_forward/MC_THERMAL_CKFU_HR"        ,"mc_thermal_forward/DISORT_THERMAL_CKFU_HR"]    ,"mc_thermal_forward_mc_ckfu_hr"         ,0.1    ,0.01 ,1,$found_mystic3d],
 ["mc_thermal_forward_disort_Ag0"         ,["mc_thermal_forward/DISORT_THERMAL_AG.0"       ,]                                               ,"mc_thermal_forward_disort_Ag0"         ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_Ag0"             ,["mc_thermal_forward/MC_THERMAL_AG.0"           ,"mc_thermal_forward/DISORT_THERMAL_AG.0"]       ,"mc_thermal_forward_mc_Ag0"             ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_disort_Ag0.5"       ,["mc_thermal_forward/DISORT_THERMAL_AG.5"       ,]                                               ,"mc_thermal_forward_disort_Ag0.5"       ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_Ag0.5"           ,["mc_thermal_forward/MC_THERMAL_AG.5"           ,"mc_thermal_forward/DISORT_THERMAL_AG.5"]       ,"mc_thermal_forward_mc_Ag0.5"           ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_disort_ssa1"        ,["mc_thermal_forward/DISORT_THERMAL_SSA1"       ,]                                               ,"mc_thermal_forward_disort_ssa1"        ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_ssa1"            ,["mc_thermal_forward/MC_THERMAL_SSA1"           ,"mc_thermal_forward/DISORT_THERMAL_SSA1"]       ,"mc_thermal_forward_mc_ssa1"            ,0.00001,0.01 ,1,"no_constraint"],
 ["mc_thermal_forward_disort_Ag0_noscat"  ,["mc_thermal_forward/DISORT_THERMAL_AG.0_NOSCAT",]                                               ,"mc_thermal_forward_disort_Ag0_noscat"  ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_Ag0_noscat"      ,["mc_thermal_forward/MC_THERMAL_AG.0_NOSCAT"    ,"mc_thermal_forward/DISORT_THERMAL_AG.0_NOSCAT"],"mc_thermal_forward_mc_Ag0_noscat"      ,0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_disort_Ag0.5_noscat",["mc_thermal_forward/DISORT_THERMAL_AG.5_NOSCAT",]                                               ,"mc_thermal_forward_disort_Ag0.5_noscat",0.00001,0.001,1,"no_constraint"],
 ["mc_thermal_forward_mc_Ag0.5_noscat"    ,["mc_thermal_forward/MC_THERMAL_AG.5_NOSCAT"    ,"mc_thermal_forward/DISORT_THERMAL_AG.5_NOSCAT"],"mc_thermal_forward_mc_Ag0.5_noscat"    ,0.00001,0.001,1,"no_constraint"],
);
foreach my $test (@mc_thermal_forward_tests) { run_test_entry(\&uvspec_run_flx_test, $test); }

my @schwarzschild_tests = (
 ["schwarzschild_thermal_Ag0"  ,["mc_thermal_forward/SCHWARZSCHILD_THERMAL_AG.0","mc_thermal_forward/DISORT_THERMAL_AG.0_NOSCAT"],"schwarzschild_thermal_Ag0"  ,0.00001,0.001,1,"no_constraint"],
 ["schwarzschild_thermal_Ag0.5",["mc_thermal_forward/SCHWARZSCHILD_THERMAL_AG.5","mc_thermal_forward/DISORT_THERMAL_AG.5_NOSCAT"],"schwarzschild_thermal_Ag0.5",0.00001,0.001,1,"no_constraint"],
);
foreach my $test (@schwarzschild_tests) { run_test_entry(\&uvspec_run_flx_test, $test); }

my @mc_flx_abs_tests = (
 ["mc_thermal_forward_mcflxabs_ckfu_hr",["mc_thermal_forward/MC_THERMAL_CKFU_HR",] ,"mc_thermal_forward_mcflxabs_ckfu_hr", 0.1, 0.05, 1, $found_mystic3d],
 ["mc_solar_forward_ckfu_hr"           ,["mc_solar_forward/MC_SOLAR_CKFU_HR"    ,] ,"mc_solar_forward_ckfu_hr"           , 1.0, 0.2 , 1, $found_mystic3d],
);
foreach my $test (@mc_flx_abs_tests) { run_test_entry(\&mystic_run_flx_abs_test, $test); }


my @sbdart_thermal_tests = (
 ["uvspec_sbdart_thermal"        ,["UVSPEC_LOWTRAN_THERMAL"              ,],"uvspec_sbdart_thermal"        ,0.00001,0.01,1,$found_lowtran],
 ["uvspec_sbdart_thermal_ftwostr",["tests/UVSPEC_LOWTRAN_THERMAL_FTWOSTR",],"uvspec_sbdart_thermal_ftwostr",0.00001,0.01,2,$found_lowtran],
);
foreach my $test (@sbdart_thermal_tests) { run_test_entry(\&uvspec_test_sbdart_thermal, $test); }

my @mystic_run_tests = (
 ["uvspec_mystic_backward_bpdf"        ,["UVSPEC_MC_BACKWARD_BPDF","UVSPEC_MC_BPDF"],"UVSPEC_MC_BPDF"                                 ,0.00001,0.01 ,1,$found_mystic                                  ],
 ["uvspec_mystic_backward_polarisation",["UVSPEC_MC_BACKWARD_POL" ,"UVSPEC_MC_POL" ],"UVSPEC_MC_POL"                                  ,0.00001,0.03 ,1,$found_mystic && $found_tests                  ],
 ["uvspec_mystic_bpdf"                 ,["UVSPEC_MC_BPDF"         ,                ],"MYSTIC polarized surface reflection (BPDF)"     ,0.00001,0.01 ,1,$found_mystic                                  ],
 ["uvspec_mystic_cam_polarized"        ,["UVSPEC_MC_CAM_POLARIZED",                ],"MYSTIC polarized surface reflection (BPDF) CAM ",0.00001,0.001,1,$found_mystic                                  ],
 ["uvspec_mystic_delta_scaling"        ,["UVSPEC_MC_DELTASCALING" ,                ],"MYSTIC delta scaling"                           ,0.00001,0.01 ,1,$found_vroom && $found_mystic3d&& $found_wcmie ],
 ["uvspec_mystic_no_polarisation"      ,["UVSPEC_MC_RAD"          ,                ],"MYSTIC without polarisation"                    ,0.00001,0.01 ,1,$found_mystic && $found_tests                  ],
 ["uvspec_mystic_polarisation"         ,["UVSPEC_MC_POL"          ,                ],"MYSTIC polarisation"                            ,0.00001,0.03 ,1,$found_mystic                                  ],

 ["uvspec_mystic_vroom"                ,["UVSPEC_MC_VROOM"        ,                ],"MYSTIC vroom"                                   ,0.00001,0.015,1,$found_vroom && $found_mystic3d && $found_wcmie],
);
foreach my $test (@mystic_run_tests) { run_test_entry(\&mystic_run_test_mc_rad, $test); }


my @mystic_run_rad_tests = (
    ["uvspec_mystic_core"                 ,["UVSPEC_MC_CORE_INPUT"           ,"UVSPEC_MC_CORE"],"MYSTIC hidden core delta scaling"  ,0.00001,0.02 ,1,$found_mystic3d && $found_wcmie                         ],
    ["uvspec_mystic_core_corefinder"      ,["UVSPEC_MC_CORE"                 ,                ],"MYSTIC hidden corefinder from file",0.00001,0.02 ,1,$found_mystic3d && $found_wcmie                         ],
    ["uvspec_mystic_sunphotometer"        ,["UVSPEC_MC_SUNPHOTOMETER"        ,                ],"MYSTIC sunphotometer"              ,0.00001,0.015,1,$found_vroom && $found_mystic3d                         ],
    ["uvspec_mystic_sunshape"             ,["UVSPEC_MC_SUNSHAPE"             ,                ],"MYSTIC sunshape"                   ,0.00001,0.03 ,1,$found_vroom && $found_mystic3d && $found_reptran_coarse],
    ["uvspec_mystic_spherical"            ,["UVSPEC_MC_SPHERICAL"            ,                ],"MYSTIC spherical"                  ,0.00001,0.02 ,1,$found_reptran_coarse                                   ],
    ["uvspec_mystic_ambrals_spectral_file",["UVSPEC_MC_AMBRALS_SPECTRAL_FILE",                ],"MYSTIC AMBRALS spectral file"      ,0.0005 ,0.001,1,$found_mystic3d                                         ],
    #broken ["uvspec_mystic_blitz" ,["UVSPEC_MC_BLITZ" ,] ,"MYSTIC blitz" ,0.00001,0.1 ,1,$found_vroom],
    ["uvspec_mystic_rgb_panorama_file"    ,["UVSPEC_MC_RGB_PANORAMA_FILE"    ,                ],"MYSTIC output_process rgb + pano"  ,0.0005 ,0.05 ,1,$found_mystic3d && $found_vroom                         ],
    ["uvspec_mystic_backward_file"        ,["UVSPEC_MC_BACKWARD_FILE"        ,                ],"MYSTIC backward rays from file"    ,0.0005 ,0.05 ,1,$found_mystic3d && $found_vroom                         ],
    ["uvspec_mystic_cloudglint"           ,["UVSPEC_MC_CLOUDGLINT"           ,                ],"MYSTIC cloudglint parameterization",0.5    ,0.05, 1,$found_mystic3d && $found_vroom && $found_baum_v36      ],
    ["uvspec_mystic_spherical3d_hamster"  ,["UVSPEC_MC_SPHERICAL3D_HAMSTER"  ,                ],"MYSTIC spherical3d HAMSTER albedo" ,0.5    ,0.05 ,1,$found_mystic3d && $found_vroom && $found_hamster       ],
    );
foreach my $test (@mystic_run_rad_tests) { run_test_entry(\&mystic_run_rad_test, $test); }


my @mystic_run_backward_tests = (
 ["uvspec_mystic_backward_edir",["UVSPEC_MC_BACKWARD_EDIR",],"MYSTIC backward direct irradiance"       ,0.00001,0.02,1,$found_mystic3d && $found_tests,3],
 ["uvspec_mystic_backward_edn" ,["UVSPEC_MC_BACKWARD_EDN" ,],"MYSTIC backward downwelling irradiance"  ,0.00001,0.02,1,$found_mystic3d && $found_tests,4],
 ["uvspec_mystic_backward_eup" ,["UVSPEC_MC_BACKWARD_EUP" ,],"MYSTIC backward upwelling irradiance"    ,0.00001,0.02,1,$found_mystic3d && $found_tests,5],
 ["uvspec_mystic_backward_fdir",["UVSPEC_MC_BACKWARD_FDIR",],"MYSTIC backward direct actinic flux"     ,0.00001,0.02,1,$found_mystic3d && $found_tests,6],
 ["uvspec_mystic_backward_fdn" ,["UVSPEC_MC_BACKWARD_FDN" ,],"MYSTIC backward downwelling actinic flux",0.00001,0.02,1,$found_mystic3d && $found_tests,7],
 ["uvspec_mystic_backward_fup" ,["UVSPEC_MC_BACKWARD_FUP" ,],"MYSTIC backward upwelling actinic flux"  ,0.00001,0.02,1,$found_mystic3d && $found_tests,8],
);
foreach my $test (@mystic_run_backward_tests) { run_test_entry(\&mystic_run_bw_test, $test); }


my @mystic_run_t_tests = (
 ["uvspec_mystic_backward_heat_hybrid" ,["UVSPEC_MC_BACKWARD_HEAT_HYBRID" ,],"MYSTIC backward heating rate hybrid" ,0.00001,0.01 ,1,$found_mystic3d && $found_tests,5],
 ["uvspec_mystic_backward_thermal_abs" ,["UVSPEC_MC_BACKWARD_THERMAL_ABS" ,],"MYSTIC backward thermal absorption"  ,0.00001,0.01 ,1,$found_mystic3d && $found_tests,5],
 ["uvspec_mystic_backward_thermal_act" ,["UVSPEC_MC_BACKWARD_THERMAL_ACT" ,],"MYSTIC backward thermal actinic flux",0.00001,0.01 ,1,$found_mystic3d && $found_tests,5],
 ["uvspec_mystic_backward_thermal_emis",["UVSPEC_MC_BACKWARD_THERMAL_EMIS",],"MYSTIC backward thermal emission"    ,0.00001,0.0  ,1,$found_mystic3d && $found_tests,5],
 ["uvspec_nca2016"                     ,["UVSPEC_HEAT_NCA2016"            ,],"NCA 1.0"                             ,0.00001,0.001,1,$found_mystic                  ,5],
 ["uvspec_nca2019_cuboid"              ,["UVSPEC_HEAT_NCA2019_CUBOID"     ,],"NCA 2.0 cuboid"                      ,0.00001,0.001,1,$found_mystic                  ,5],
 ["uvspec_nca2019_triangle"            ,["UVSPEC_HEAT_NCA2019_TRIANGLE"   ,],"NCA 2.0 triangle"                    ,0.00001,0.001,1,$found_mystic                  ,5],
);
foreach my $test (@mystic_run_t_tests) { run_test_entry(\&mystic_run_t_test, $test); }

my @uvspec_run_3d_abs_tests=(
    ["uvspec_dts_solar",["UVSPEC_DTS_SOLAR",],"dynamic tenstream solar"  ,0.00001,0.001,1,$found_dts,5],
    ["uvspec_dts_therm",["UVSPEC_DTS_THERM",],"dynamic tenstream thermal",0.00001,0.001,1,$found_dts,5],
    );
foreach my $test (@uvspec_run_3d_abs_tests) { run_test_entry(\&uvspec_run_3d_abs_test, $test); }
 
my @mystic_run_is_IQ_tests = (
 ["uvspec_mystic_spectral"         ,["UVSPEC_MC_SPECTRAL"         ,                    ],"MYSTIC spectral importance sampling (line by line and BPDF)",0.00001,0.05,1,$found_mystic && $found_opac ],
 ["uvspec_mystic_backward_spectral",["UVSPEC_MC_BACKWARD_SPECTRAL","UVSPEC_MC_SPECTRAL"],"MYSTIC backward spectral importance sampling"               ,0.00001,0.01,1,$found_mystic && $found_opac ],
 ["uvspec_mystic_spectral_reptran" ,["UVSPEC_MC_SPECTRAL_REPTRAN" ,                    ],"MYSTIC spectral importance sampling (REPTRAN)"              ,0.00001,0.05,1,$found_mystic                ],
 ["uvspec_mystic_earthshine"       ,["UVSPEC_MC_EARTHSHINE"       ,                    ],"MYSTIC Earthshine"                                          ,0.00001,0.05,1,$found_mystic3d && $found_hey],
);
foreach my $test (@mystic_run_is_IQ_tests) { run_test_entry(\&mystic_run_is_IQ_test, $test); }

# CE: these run also in 1D but should not be in public version, therefor found_mystic3d included 
my @mystic_jacobian_tests = (
 ["uvspec_mystic_jacobian_novroom",["UVSPEC_MC_JACOBIAN_NOVROOM",],"MYSTIC jacobians VROOM disabled",0.005,0.3,1,$found_mystic && $found_tests && $found_mystic3d               ],
 ["uvspec_mystic_jacobian_vroom"  ,["UVSPEC_MC_JACOBIAN_VROOM"  ,],"MYSTIC jacobians VROOM enabled" ,0.005,0.3,1,$found_mystic && $found_vroom && $found_tests && $found_mystic3d],
 ["uvspec_mystic_jacobian_edn"    ,["UVSPEC_MC_JACOBIAN_EDN"    ,],"MYSTIC jacobians irradiance"    ,0.005,0.5,1,$found_mystic && $found_vroom && $found_tests && $found_mystic3d],
);
foreach my $test (@mystic_jacobian_tests) { run_test_entry(\&mystic_run_jacobian_test, $test); }

my @mystic_mc_flx_spc_tests = (
 ["mc_thermal_forward_elevation"     ,["mc_thermal_forward/MC_THERMAL_ELEVATION"    ,""]   ,"mc_thermal_forward_elevation"     ,0.0001,0.01 ,1,$found_mystic3d],
 ["mc_thermal_forward_elevation_par" ,["mc_thermal_forward/MC_THERMAL_ELEVATION_PAR",""]   ,"mc_thermal_forward_elevation_par" ,0.0001,0.01 ,1,$found_mystic3d],
);
foreach my $test (@mystic_mc_flx_spc_tests) { run_test_entry(\&mystic_run_test_mc_flx_spc, $test); }

# thermal backward tests, iterating over heights and methods
my @mc_thermal_backward_heights = ( "5.000" );
if ($all == 3){
my @mc_thermal_backward_heights = ( "0.000",  "1.000",  "2.000",  "3.000",  "4.000",  "5.000",  "6.000",  "7.000",  "8.000",  "9.000",
                                    "10.000", "11.000", "12.000", "13.000", "14.000", "15.000", "16.000", "17.000", "18.000", "19.000");
}


my @mc_thermal_backward_method = ("hybrid", "emabs", "emabsopt", ); # "denet"); # not running denet because too noisy for heating rates

my @mc_thermal_backward_tests = ();
foreach my $bw_method (@mc_thermal_backward_method) {
  foreach my $height (@mc_thermal_backward_heights) {
    push (@mc_thermal_backward_tests, ["mc_thermal_backward_${bw_method}_edn_$height", ["mc_thermal_backward/MC_${bw_method}_$height.EDN",], "mc_thermal_backward_${bw_method}_edn_$height", 1e-5,0.05 ,1,"no_constraint"] );}

  foreach my $height (@mc_thermal_backward_heights) {
    push (@mc_thermal_backward_tests, ["mc_thermal_backward_${bw_method}_eup_$height", ["mc_thermal_backward/MC_${bw_method}_$height.EUP",], "mc_thermal_backward_${bw_method}_eup_$height", 1e-5,0.05 ,1,$found_mystic3d] );}
}
foreach my $test (@mc_thermal_backward_tests) { run_test_entry(\&mystic_run_test_mc_flx_spc, $test); }

my @mc_thermal_backward_heat_tests = ();
foreach my $bw_method (@mc_thermal_backward_method) {
  foreach my $height (@mc_thermal_backward_heights) {
    push (@mc_thermal_backward_heat_tests, ["mc_thermal_backward_${bw_method}_heat_$height", ["mc_thermal_backward/MC_${bw_method}_$height.HEAT",], "mc_thermal_backward_${bw_method}_heat_$height", 1e-5,0.1 ,1,$found_mystic3d,5] );}
}
foreach my $test (@mc_thermal_backward_heat_tests) { run_test_entry(\&mystic_run_t_test, $test); }


run_test_entry(\&mystic_run_amf_test  ,["uvspec_mystic_boxairmass"          ,["UVSPEC_MC_BOXAIRMASS"          ,],"MYSTIC boxairmass factor"                             ,0.00001,0.01,1,$found_mystic                                ]);
run_test_entry(\&mystic_run_irr_test  ,["uvspec_mystic_irradiance"          ,["UVSPEC_MC_IRRADIANCE"          ,],"MYSTIC irradiance"                                    ,0.00001,0.05,1,$found_mystic3d                              ]);
run_test_entry(\&mystic_run_aeris_test,["uvspec_mystic_aerosol_is_edn"      ,["UVSPEC_MC_AERIS_EDN"           ,],"MYSTIC aerosol importance sampling (solar irradiance)",0.001  ,0.05,1,$found_mystic3d                              ]);
run_test_entry(\&mystic_run_cb_test   ,["uvspec_mystic_coherent_backscatter",["UVSPEC_MC_COHERENT_BACKSCATTER",],"MYSTIC coherent backscatter"                          ,0.01   ,0.05,1,$found_mystic3d && $found_wcmie                ]);
run_test_entry(\&radar_run_test       ,["uvspec_mystic_radar"               ,["UVSPEC_MC_RADAR"               ,],"MYSTIC radar with ris factor"                         ,1.7e-14,0.01,1,$found_lidar && $found_wcmie                 ]);
run_test_entry(\&ssradar_test         ,["ssradar"                           ,["SSRADAR"                       ,],"SSRadar Test"                                         ,0.01   ,0.0 ,1,"no_constraint"                              ]);

# expensive, jacobians have huge noise, we are happy with factor 5:
run_test_entry(\&lidar_run_test, ["uvspec_mystic_lidar", ["UVSPEC_MC_LIDAR", ], "MYSTIC lidar" ,1e-14, 5, 1, $found_lidar && $found_wcmie]);


my @mie_tests = (
 ["mie_v0_wc_iapws"         ,["test_mie/MIE_V0_WC_IAPWS"      ,],"mie water IAPWS table"         ,0.00001,0.001  ,1,"no_constraint"],
 ["mie_v0_wc_segelstein"    ,["test_mie/MIE_V0_WC_SEGELSTEIN" ,],"mie water SEGELSTEIN table"    ,0.00001,0.005  ,1,"no_constraint"],
 ["mie_v0_ic"               ,["test_mie/MIE_V0_IC"            ,],"mie ice table"                 ,0.00001,0.001  ,1,"no_constraint"],
);
foreach my $test (@mie_tests) { run_test_entry(\&mie_run_test, $test); }



#
# Some more custom test function that dont fit the run_test_entry template:
#
#

if ( $test eq "make_slitfunction" || $all >= 1 ) {
  if ($dryrun) {
    print("\nWould run testcase make_slitfunction\n");
  } else{
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    printf stderr "make_slitfunction test\n";
    $fail += make_slitfunction_test($answers,$rtol,$atol);
  }
}

if ( ($test eq "make_angresfunc" || $all >= 1) && $found_libgsl) {
  if ($dryrun) {
    print("\nWould run testcase make_angresfunc\n");
  } else {
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    printf stderr "make_angresfunc test\n";
    $fail += make_angresfunc_test($answers,$rtol,$atol);
  }
}

if ( ($test eq "angres" || $all >= 1) && $found_libgsl) {
  if ($dryrun) {
    print("\nWould run testcase angres\n");
  } else {
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    printf stderr "angres test\n";
    $fail += angres_test($answers,$rtol,$atol);
  }
}

if ( $test eq "uvspec_redistribution" || $all >= 1 )  {
  if ($dryrun) {
    print("\nWould run testcase uvspec_redistribution\n");
  } else {
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    $fail += uvspec_test_redistribution($answers,$rtol,$atol);
  }
}

if ( ($test eq "atm_reflectivity" || $all >= 1) && $found_reptran_coarse ) {
  if ($dryrun) {
    print("\nWould run testcase atm_reflectivity\n");
  } else {
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    $fail += atm_reflectivity_test($answers,$rtol,$atol);
  }
 }

if ( $test eq "uvspec_saturate_clouds" ||
     $test eq "saturate clouds column A" ||
     $test eq "saturate clouds column B" ||
     $test eq "saturate clouds column C" ||
     $test eq "saturate clouds ipa" || $all >= 1 )  {
  if ($dryrun) {
    print("\nWould run testcase uvspec_saturate_clouds\n");
  } else {
    my ($atol, $rtol) = get_tolerances(0.00001, 0.001);
    $fail += uvspec_test_saturate_clouds($answers,$rtol,$atol,$test);
  }
}


# TENSTREAM Tests
my $tenstream_info = "";

# Check that we can call the tenstream lib with its 2stream solver.
# The n-stream solvers need Lookup tables which would need to be copied or created.
# However, this takes too long and we are therefore content at this point
# that we end up in the solver routine via the delta eddington twostream solver.
# Things we check here are, that the rank0toAll scatter and gathers work and we have a sensible domain decomposition.
my $tenstream_options="";

if ( ("tenstream_1D_box70_sw" =~ $test || $all >= 1) && $found_tenstream && $found_mystic )    {
  if ($dryrun) { print("\nWould run testcase tenstream_1D_box70_sw\n"); } else {

    my $spec = 'solar';
    my @szas = (0,60);        # run test for various zenith angles
    my @nrCPUs = (1,2);       # run test for various number of mpi nodes

    foreach my $sza (@szas) {
      for my $ncpu (@nrCPUs) {
        my ($atol, $rtol) = get_tolerances(0.00001, 0.08);

        my $extemplate  = "../examples/tenstream_test/BOX70.INP.TEMPLATE";
        my $exinp       = "../examples/tenstream_test/BOX70_TENSTR_1D.$spec.sza$sza.inp";
        my $exout       = "../examples/tenstream_test/BOX70_TENSTR_1D.$spec.sza$sza.out";

        `cat $extemplate           >  $exinp`;
        `echo source $spec         >> $exinp`;
        `echo sza $sza             >> $exinp`;
        if ( $answers == 1 ) {
          $ncpu = 1;
          `echo mc_ipa               >> $exinp`;
        } else {
          `echo mc_tenstream solver_2str >> $exinp`;
        }

        my $titlestr = sprintf("tenstream_1D_box70_sw %s with NCPU=%4d for sza=%3d",$spec,$ncpu,$sza) . $tenstream_info;
        $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
      } #ncpu
    } #sza
  }
}

if ( ("tenstream_1D_box70_lw" =~ $test || $all >= 1) && $found_tenstream && $found_mystic )    {
  if ($dryrun) { print("\nWould run testcase tenstream_1D_box70_lw\n"); } else {

    my $spec = 'thermal';
    my @nrCPUs = (1,2);       # run test for various number of mpi nodes

    for my $ncpu (@nrCPUs) {
      my ($atol, $rtol) = get_tolerances(0.00001, 0.05);

      my $extemplate  = "../examples/tenstream_test/BOX70.INP.TEMPLATE";
      my $exinp       = "../examples/tenstream_test/BOX70_TENSTR_1D.$spec.inp";
      my $exout       = "../examples/tenstream_test/BOX70_TENSTR_1D.$spec.out";

      `cat $extemplate           >  $exinp`;
      `echo source $spec         >> $exinp`;
      if ( $answers == 1 ) {
        $ncpu = 1;
        `echo mc_ipa               >> $exinp`;
      } else {
        $tenstream_options="-schwarzschild -schwarzschild_Nmu 2";
        `echo mc_tenstream solver_2str >> $exinp`;
      }

      my $titlestr = sprintf("tenstream_1D_box70_lw %s with NCPU=%4d",$spec,$ncpu) . $tenstream_info;
      $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
    } #ncpu
  }
}

if ( ("tenstream_3D_box70_sw" =~ $test || $all >= 4) && $found_tenstream && $found_mystic )    {
  if ($dryrun) { print("\nWould run testcase tenstream_3D_box70_sw\n"); } else {

    my $spec = 'solar';
    my @szas = (0,60);        # run test for various zenith angles
    my @nrCPUs = (1,2);       # run test for various number of mpi nodes
    # 3D
    foreach my $sza (@szas) {
      for my $ncpu (@nrCPUs) {
        my ($atol, $rtol) = get_tolerances(0.00001, 0.7);

        my $extemplate  = "../examples/tenstream_test/BOX70.INP.TEMPLATE";
        my $exinp       = "../examples/tenstream_test/BOX70_TENSTR_3D.$spec.sza$sza.inp";
        my $exout       = "../examples/tenstream_test/BOX70_TENSTR_3D.$spec.sza$sza.out";

        `cat $extemplate           >  $exinp`;
        `echo source $spec         >> $exinp`;
        `echo sza $sza             >> $exinp`;
        if ( $answers == 1 ) {
          $ncpu = 1;
        } else {
          `echo mc_tenstream solver_3_10 >> $exinp`;
        }

        my $titlestr = sprintf("tenstream_3D_box70_sw NCPU=%4d for sza=%3d",$ncpu,$sza) . $tenstream_info;
        $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
      } #ncpu
    } #sza
  }
}

if ( ("tenstream_3D_box70_lw" =~ $test || $all >= 4) && $found_tenstream && $found_mystic )    {
  if ($dryrun) { print("\nWould run testcase tenstream_3D_box70_lw\n"); } else {

  my $spec = 'thermal';
  my @nrCPUs = (1,2);       # run test for various number of mpi nodes
  # 3D
  for my $ncpu (@nrCPUs) {
    my ($atol, $rtol) = get_tolerances(0.000001, 0.15);

    my $extemplate  = "../examples/tenstream_test/BOX70.INP.TEMPLATE";
    my $exinp       = "../examples/tenstream_test/BOX70_TENSTR_3D.$spec.inp";
    my $exout       = "../examples/tenstream_test/BOX70_TENSTR_3D.$spec.out";

    `cat $extemplate           >  $exinp`;
    `echo source $spec         >> $exinp`;
    if ( $answers == 1 ) {
      $ncpu = 1;
    } else {
      `echo mc_tenstream solver_3_10 >> $exinp`;
    }

    my $titlestr = sprintf("tenstream_3D_box70_lw NCPU=%4d",$ncpu) . $tenstream_info;
    $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
  } #ncpu
}
}

# if ( ("tenstream_vs_mc_phillip1_sw" =~ $test || $all >= 1) && $found_tenstream && $found_mystic )    {
#   # Test Options:
#
#   my $spec = 'solar';
#   my @szas   = (60);                 # run test for various zenith angles
#   my @nrCPUs = (8);                  # run test for various number of mpi nodes
#   # 3D
#   foreach my $sza (@szas) {
#     for my $ncpu (@nrCPUs) {
#       my ($atol, $rtol) = get_tolerances(0.00001, 0.0001);
#
#       my $extemplate  = "../examples/tenstream_test/PHILLIP1.INP.TEMPLATE";
#       my $exinp       = "../examples/tenstream_test/PHILLIP1.$spec.sza$sza.inp";
#       my $exout       = "../examples/tenstream_test/PHILLIP1.$spec.sza$sza.out";
#
#       `cat $extemplate           >  $exinp`;
#       `echo source $spec         >> $exinp`;
#       `echo sza $sza             >> $exinp`;
#
#       if ( $answers == 1 ) {
#         $ncpu = 1
#       } else {
#         `echo mc_tenstream solver_3_10 >> $exinp`;
#       }
#
#       my $titlestr = sprintf("tenstream vs mystic ($spec) with NCPU=%4d for sza=%3d",$ncpu,$sza) . $tenstream_info;
#       $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
#     } #ncpu
#   } #sza
# }
#
# if ( ("tenstream_vs_mc_phillip1_lw" =~ $test || $all >= 1) && $found_tenstream && $found_mystic )    {
#   # Test Options:
#
#   my $spec = 'thermal';
#   my @nrCPUs = (8);                  # run test for various number of mpi nodes
#   # 3D
#   for my $ncpu (@nrCPUs) {
#     my ($atol, $rtol) = get_tolerances(0.00001, 0.0001);
#
#     my $extemplate  = "../examples/tenstream_test/PHILLIP1.INP.TEMPLATE";
#     my $exinp       = "../examples/tenstream_test/PHILLIP1.$spec.inp";
#     my $exout       = "../examples/tenstream_test/PHILLIP1.$spec.out";
#
#     `cat $extemplate           >  $exinp`;
#     `echo source $spec         >> $exinp`;
#
#     if ( $answers == 1 ) {
#       $ncpu = 1
#     } else {
#       `echo mc_tenstream solver_3_10 >> $exinp`;
#     }
#
#     my $titlestr = sprintf("tenstream vs mystic ($spec) with NCPU=%4d",$ncpu) . $tenstream_info;
#     $fail += uvspec_mpi_run_spc_test($inp, $out, $exinp, $exout, $titlestr, $answers,$rtol,$atol,$ncpu,$tenstream_options);
#   } #ncpu
# }


# Triangle Surface Tests
if ( ("triangle_surface_ncflx_forward_flat_surface_solar" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_ncflx_forward_flat_surface_solar\n"); } else {
    my ($atol, $rtol) = (30, 2e-2);
    my $exinp = "../examples/triangular_srfc/forward_flat_surface/uvspec.inp";
    my $exout = "../examples/triangular_srfc/forward_flat_surface/99_flat_surface_solar.testuite_output.nc";
    my ($f, $fdiff) = mystic_run_triangle_flx_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_ncflx_forward_flat_surface_solar", $answers, $atol, $rtol);
    $fail += $f; $fail_diff += $fdiff;
  }
}
if ( ("triangle_surface_ncflx_forward_flat_surface_thermal_noatm" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_ncflx_forward_flat_surface_thermal_noatm\n"); } else {
    my ($atol, $rtol) = (1e-5, 1e-7);
    my $exinp = "../examples/triangular_srfc/forward_flat_surface_thermal/no_atmosphere.inp";
    my $exout = "../examples/triangular_srfc/forward_flat_surface_thermal/99_flat_surface_thermal_no_atmosphere.testuite_output.nc";
    my ($f, $fdiff) = mystic_run_triangle_flx_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_ncflx_forward_flat_surface_thermal_noatm", $answers, $atol, $rtol);
    $fail += $f; $fail_diff += $fdiff;
  }
}
if ( ("triangle_surface_ncflx_forward_flat_surface_thermal_atm" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_ncflx_forward_flat_surface_thermal_atm\n"); } else {
    my ($atol, $rtol) = (1e-3, 1e-2);
    my $exinp = "../examples/triangular_srfc/forward_flat_surface_thermal/atmosphere.inp";
    my $exout = "../examples/triangular_srfc/forward_flat_surface_thermal/99_flat_surface_thermal_atmosphere.testuite_output.nc";
    my ($f, $fdiff) = mystic_run_triangle_flx_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_ncflx_forward_flat_surface_thermal_atm", $answers, $atol, $rtol);
    $fail += $f; $fail_diff += $fdiff;
  }}

if ( ("triangle_surface_camera_flat_surface_solar" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_camera_flat_surface_solar\n"); } else {
    my ($atol, $rtol) = get_tolerances(1e-2, 1e-2);
    my $exinp = "../examples/triangular_srfc/pano_flat_surface_solar/uvspec.inp";
    my $exout = "../examples/triangular_srfc/pano_flat_surface_solar/99_flat_surface_solar.testuite_output";
    my ($f, $fdiff) = mystic_run_rad_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_camera_flat_surface_solar", $answers, $rtol, $atol);
    $fail += $f; $fail_diff += $fdiff;
  }
}
if ( ("triangle_surface_camera_flat_surface_thermal" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_camera_flat_surface_thermal\n"); } else {
    my ($atol, $rtol) = get_tolerances(1e-2, 1e-2);
    my $exinp = "../examples/triangular_srfc/pano_flat_surface_thermal/uvspec.inp";
    my $exout = "../examples/triangular_srfc/pano_flat_surface_thermal/99_flat_surface_thermal.testuite_output";
    my ($f, $fdiff) = mystic_run_rad_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_camera_flat_surface_thermal", $answers, $rtol, $atol);
    $fail += $f; $fail_diff += $fdiff;
  }
}
if ( ("triangle_surface_camera_vertical_wall_solar" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_camera_vertical_wall_solar\n"); } else {
    my ($atol, $rtol) = get_tolerances(1e-0, 5e-2);
    my $exinp = "../examples/triangular_srfc/pano_vertical_wall_solar/uvspec.inp";
    my $exout = "../examples/triangular_srfc/pano_vertical_wall_solar/99_vertical_wall_solar.testuite_output";
    my ($f, $fdiff) = mystic_run_rad_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_camera_vertical_wall_solar", $answers, $rtol, $atol);
    $fail += $f; $fail_diff += $fdiff;
  }
}
if ( ("triangle_surface_camera_vertical_wall_thermal" =~ $test || $all >= 1) && $found_mystic && $found_mystic3d) {
  if ($dryrun) { print("\nWould run testcase triangle_surface_camera_vertical_wall_thermal\n"); } else {
    my ($atol, $rtol) = get_tolerances(1e-3, 2e-2);
    my $exinp = "../examples/triangular_srfc/pano_vertical_wall_thermal/uvspec.inp";
    my $exout = "../examples/triangular_srfc/pano_vertical_wall_thermal/99_vertical_wall_thermal.testuite_output";
    my ($f, $fdiff) = mystic_run_rad_test($inp, $out, $exinp, $exout, "MYSTIC triangle_surface_camera_vertical_wall_thermal", $answers, $rtol, $atol);
    $fail += $f; $fail_diff += $fdiff;
  }
}

{
  my ($atol, $rtol) = get_tolerances(1e-5, 1e-3);
  if ( $test eq "airmass_test_bro"  || $test eq "airmass_test_bro" || $all > 3 )   { printf stderr "airmass.pl test\n"; }
  if ( $test eq "airmass_test_bro"  || $all > 3 ) {
    if ($dryrun) { print("\nWould run testcase airmass_test_bro\n"); } else {
      $fail += airmass_test_bro ($answers,$rtol,$atol);
    }
  }
  if ( $test eq "airmass_test_hcho" || $all > 3 ) {
    if ($dryrun) { print("\nWould run testcase airmass_test_hcho\n"); } else {
      $fail += airmass_test_hcho($answers,$rtol,$atol);
    }
  }

  if ( $test eq "calc_j_test" || $all > 3 ) {
    if ($dryrun) { print("\nWould run testcase calc_j_test\n"); } else {
      $fail += calc_j_test($answers,$rtol,$atol);
    }
  }
}

# run optical properties only if full distribution is available, not in public distribution (therefore included mystic3d)
if ( ("optical_properties_test" =~ $test && $all == 3)
  && $found_mystic3d
  && ( $found_nco || $found_nccmp )
  && $found_hey
  && $found_baum
  && $found_wcmie
) {
  if ($dryrun) { print("\nWould run testcase optical_properties_test\n"); } else {
    printf stderr "Optical properties test\n";
    if(!$valgrind) {
      $fail += system("bash ../test/run_optical_properties_test.sh");
    } else {
      $fail += system("bash ../test/run_optical_properties_test.sh --valgrind");
    }
  }
}

if ($dryrun) { exit 0; }

my $cmp = $test cmp "uvspec";
if ( $cmp >= 0 || $all >= 1 ) { report_results($answers,$fail,"uvspec calls"); }
if ( $cmp >= 0 || $all >= 1 ) { report_results($answers,$fail_diff,"difference"); }

my $global_time = time - $global_t_start;
print("\nTest suite took $global_time seconds to run\n");

exit $fail + $fail_diff;
#
#
# End of test runner
#
#

sub get_tolerances {
  my ($test_atol,  $test_rtol) = @_;
  if ($global_atol >= 0) { $test_atol = $global_atol; }
  if ($global_rtol >= 0) { $test_rtol = $global_rtol; }
  return ($test_atol,  $test_rtol);
}

sub run_test_entry { # runs a testdefinition with a given test function
  my ($testmethod, $test_definition) = @_;
  my ($testname, $expfiles, $descr, $test_atol, $test_rtol, $testlvl, $constraint, $col) = @{$test_definition};
  my $num_args = scalar @{$test_definition};
  if ($num_args < 7) { die "run_test_entry has to have at least 6 arguments, $testname has only $num_args" };
  my ($inpfile, $outfile) = @{$expfiles};
  if (!$outfile) { $outfile = $inpfile; }

  if ( ($testname =~ $test || $all >= $testlvl) && $constraint ) {
    my ($atol, $rtol) = get_tolerances($test_atol, $test_rtol);
    my $exinp  = "../examples/".$inpfile.".INP";
    my $exout  = "../examples/".$outfile.".OUT";

    if ($dryrun) {
      print("\nWould run testcase $testname : input file $exinp and compare to $exout\n");
      return;
    }
    if ($verbose) {
      print("\nRunning testcase $testname : input file $exinp and compare to $exout\n");
    }
    my ($f, $fdiff) = $testmethod->($inp, $out, $exinp, $exout, $descr, $answers, $rtol, $atol, $valgrind, $col);
    $fail += $f; $fail_diff += $fdiff;
  }
}


sub atm_reflectivity_test {
    my $answers = shift @_;
    my $failure=0;
    my $pi = 4*atan2(1,1);
    my ($atol, $rtol) = get_tolerances(1e-5, 1e-3);

    # Create solar flux file with unit flux
    open(FP,">TMP.SUN");
    for(my $i=300;$i<=340;$i++) {printf FP "%f 1.0\n", $i;}
    close(FP);

    my $test = "atmospheric reflectivity";
    my $inp = "TMP.INP";
    my $out = "TMP.OUT1";
    my $exinp  = "../examples/UVSPEC_REF1.INP";
    printtest($test,$rtol,$atol);
    system("cat  $exinp > $inp");
    my $test_t_start = time;
    my $ret = system("$UVTools::uvspec < $inp > $out");
    if ($ret != 0) {carp "uvspec failed"; return 1;}

    $exinp  = "../examples/UVSPEC_REF2.INP";
    system("cat  $exinp > $inp");
    $out = "TMP.OUT2";
    $ret = system("$UVTools::uvspec < $inp > $out");
    if ($ret != 0) {carp "uvspec failed"; return 1;}

    $exinp  = "../examples/UVSPEC_REF3.INP";
    system("cat  $exinp > $inp");
    $out = "TMP.OUT3";
    $ret = system("$UVTools::uvspec < $inp > $out");
    if ($ret != 0) {carp "uvspec failed"; return 1;}
    my $test_time = time - $test_t_start;

    open(FP1,"TMP.OUT1");
    open(FP2,"TMP.OUT2");
    open(FP3,"TMP.OUT3");
    while (my $line1 = <FP1>) {
      my $line2 = <FP2>;
      my $line3 = <FP3>;
      $line1 =~ s/^[ ]+//; $line2 =~ s/^[ ]+//; $line3 =~ s/^[ ]+//;
      my ($wvn1,$dir1,$dn1,$up1) = split(/ +/,$line1);
      my ($wvn2,$dir2,$dn2,$up2) = split(/ +/,$line2);
      my ($wvn3,$dir3,$dn3,$up3) = split(/ +/,$line3);
      $up1 = $up1/$pi; # Divide by pi to get reflectivity
      if (abs(1-$up1/(1-($dir2+$dn2)/($dir3+$dn3)))>0.0001) {
        $failure++
      }
    }
    close(FP1);
    close(FP2);
    close(FP3);
    printf stderr " (%6.1f s) ", $test_time;
    if ($failure==0) {  printf stderr " ok.\n"; }
    else { printf stderr " %d serious differences.\n", $failure;}

    `rm -f TMP.OUT1 TMP.OUT2 TMP.OUT3 TMP.INP TMP.SUN`;
    $failure;
}

sub uvspec_test_redistribution {
    my $answers = shift @_;
    my $rtol = shift @_;
    my $atol   = shift @_;
    my $uvspec_failure=0;

    my $test = "redistribution";
    my $inp  = "TMP.INP";
    my $out1 = "TMP.OUT1";
    my $exinp  = "../examples/UVSPEC_REDISTRIBUTION.INP";
    printtest($test,$rtol,$atol);
    system("cat  $exinp > $inp");
    my $test_t_start = time;
    my $ret = system("$UVTools::uvspec < $inp > $out1");
    if ($ret != 0) {carp "uvspec failed"; return 1;}

    $exinp  = "../examples/UVSPEC_REDISTRIBUTION.INP";
    open(FPI,$exinp);
    open(FPT,">".$inp);
    while (my $lsi = <FPI>) {
      $_ = $lsi;
      s/altitude 0.0/altitude 0.0 0.5/;
      printf FPT;
    }
    close (FPI);
    close (FPT);
    my $out2 = "TMP.OUT2";
    $ret = system("$UVTools::uvspec < $inp > $out2");
    if ($ret != 0) {carp "uvspec failed"; return 1;}
    my $test_time = time - $test_t_start;

    my $nox = "";
    $uvspec_failure+=report_differences($answers,$out1,$out2,$inp,$test,$nox,$rtol,$atol,$test_time);
    `rm -f TMP.OUT1 TMP.OUT2 TMP.INP`;

    $uvspec_failure;
}

sub uvspec_test_saturate_clouds {
    my $answers = shift @_;
    my $rtol = shift @_;
    my $atol   = shift @_;
    my $test    = shift @_;
    my $uvspec_failure=0;
    my $lambda = 0;
    my $dir = 0;
    my $ddn = 0;
    my $dup = 0;
    my $uavg = 0;
    my @dir_sum;
    my @ddn_sum;
    my @dup_sum;
    my @uavg_sum;

    # run column A
    if ( $test eq "uvspec_test_saturate_clouds" || $test eq "saturate clouds column A" )  {
      my $exinp  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_A.INP";
      my $exout  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_A.OUT";
      my $test   = "saturate clouds column A";
      printtest($test,$rtol,$atol);
      system("cat  $exinp > $inp");
      ($ret, my $test_time) = run_uvspec($inp, $out, $valgrind);
      if ($ret != 0) {carp "uvspec failed"; return 1;}
      open (COLUMN_A, "$out") || die "couldn't open the file!";
      while (<COLUMN_A>) {
       ($lambda, $dir, $ddn, $dup, $uavg) = split;
        $dir_sum[$lambda]  = 0.5 * $dir;
        $ddn_sum[$lambda]  = 0.5 * $ddn;
        $dup_sum[$lambda]  = 0.5 * $dup;
        $uavg_sum[$lambda] = 0.5 * $uavg;
      }
      close(COLUMN_A);
      my $nox = "";
      $uvspec_failure+=report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    }

    # run column B
    if ( $test eq "uvspec_test_saturate_clouds" || $test eq "saturate clouds column B" )  {
      my $exinp  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_B.INP";
      my $exout  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_B.OUT";
      my $test   = "saturate clouds column B";
      printtest($test,$rtol,$atol);
      system("cat  $exinp > $inp");
      ($ret, my $test_time) = run_uvspec($inp, $out, $valgrind);
      if ($ret != 0) {carp "uvspec failed"; return 1;}
      open (COLUMN_B, "$out") || die "couldn't open the file!";
      while (<COLUMN_B>) {
       ($lambda, $dir, $ddn, $dup, $uavg) = split;
        $dir_sum[$lambda]  += 0.3 * $dir;
        $ddn_sum[$lambda]  += 0.3 * $ddn;
        $dup_sum[$lambda]  += 0.3 * $dup;
        $uavg_sum[$lambda] += 0.3 * $uavg;
      }
      close(COLUMN_B);
      my $nox = "";
      $uvspec_failure+=report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    }

    # run column C
    if ( $test eq "uvspec_test_saturate_clouds" || $test eq "saturate clouds column C" )  {
      my $exinp  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_C.INP";
      my $exout  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_COLUMN_C.OUT";
      my $test   = "saturate clouds column C";
      printtest($test,$rtol,$atol);
      system("cat  $exinp > $inp");
      ($ret, my $test_time) = run_uvspec($inp, $out, $valgrind);
      if ($ret != 0) {carp "uvspec failed"; return 1;}
      my $nox = "";
      open (COLUMN_C, "$out") || die "couldn't open the file!";
      open (SUM, "> TMP_IPA_SUM.OUT") || die "couldn't open the file for writing!";
      while (<COLUMN_C>) {
        ($lambda, $dir, $ddn, $dup, $uavg) = split;
         $dir_sum[$lambda]  +=  0.2 * $dir;
         $ddn_sum[$lambda]  +=  0.2 * $ddn;
         $dup_sum[$lambda]  +=  0.2 * $dup;
         $uavg_sum[$lambda] +=  0.2 * $uavg;
         printf SUM ("  %5.3f %13.6e %13.6e %13.6e %13.6e\n",$lambda,$dir_sum[$lambda],$ddn_sum[$lambda],$dup_sum[$lambda],$uavg_sum[$lambda]);
      }
      close(COLUMN_C);
      close(SUM);
      $uvspec_failure+=report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    }

    # run ipa calculation
    if ( $test eq "uvspec_test_saturate_clouds" || $test eq "saturate clouds ipa" )  {
      my $exinp  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_FILES.INP";
      my $exout  = "../examples/UVSPEC_SATURATE_CLOUDS_WC_IC_FILES.OUT";
      my $test   = "saturate clouds ipa";
      printtest($test,$rtol,$atol);
      system("cat  $exinp > $inp");
      ($ret, my $test_time) = run_uvspec($inp, $out, $valgrind);
      if ($ret != 0) {carp "uvspec failed"; return 1;}
      my $ret = system("cp $out TMP_IPA_ALL.OUT");
      my $nox = "";
      $uvspec_failure+=report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    }

    if ( $test eq "uvspec_test_saturate_clouds" )  {
      my $nox = "";
      my $test   = "saturate clouds (intern ipa comparison)";
      printtest($test,$rtol,$atol);
      $uvspec_failure+=report_differences($answers,"TMP_IPA_SUM.OUT","TMP_IPA_ALL.OUT",$inp,$test,$nox,$rtol,$atol);
    }

    $uvspec_failure;
}

sub angres_test {
    my $answers = shift @_;
    my $rtol = shift @_;
    my $atol   = shift @_;
    my $diffout= "TMP_DIFF.OUT";
    my $test = "angres";
    my $nox = "";
    printtest($test,$rtol,$atol);
    my $inp = "";
    my $out = "TMP.OUT";
    my $exout  = "../examples/ANGRES.OUT";
    my $angfile = "../examples/ANGRES_1_ANG.DAT";
    my $raddisfile = "../examples/ANGRES_RADDIS_1.DAT";
    my $cmd = "../bin/angres -a -t 0 -r 0 -s 32 -z 0 $angfile $raddisfile";
    open(FP,">".$out) || carp "Could not open file $out";
    printf FP "# Angular response function made by executing:\n";
    printf FP "#    $cmd\n";
    close(FP);
    my $test_t_start = time;
    my $ret = system("$cmd >> $out");
    my $test_time = time - $test_t_start;
    if ($ret != 0) {carp "angres failed"; return 1;}
    my $nox = "--nox";
    return report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
}

sub make_angresfunc_test {
    my $answers = shift @_;
    my $rtol = shift @_;
    my $atol   = shift @_;
    my $diffout= "TMP_DIFF.OUT";
    my $test = "make_angresfunc";
    my $nox = "";
    printtest($test,$rtol,$atol);
    my $inp = "";
    my $out = "TMP.OUT";
    my $exout  = "../examples/ANGRES_1_ANG.DAT";
    my $cmd = "../bin/make_angresfunc -t 1 -r 1";
    open(FP,">".$out) || carp "Could not open file $out";
    printf FP "# Angular response function made by executing:\n";
    printf FP "#    $cmd\n";
    close(FP);
    my $test_t_start = time;
    my $ret = system("$cmd >> $out");
    my $test_time = time - $test_t_start;
    if ($ret != 0) {carp "make_angresfunc failed"; return 1;}
    return report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
}
sub make_slitfunction_test {
    my $answers = shift @_;
    my $rtol = shift @_;
    my $atol   = shift @_;
    my $diffout= "TMP_DIFF.OUT";
    my $test = "make_slitfunction";
    my $nox = "";
    printtest($test,$rtol,$atol);
    my $inp = "";
    my $out = "TMP.OUT";
    my $exout  = "../examples/TRI_SLIT.DAT";
    my $cmd = "../bin/make_slitfunction -f 0.75 -r 0.01 -t 1";
    open(FP,">".$out) || carp "Could not open file $out";
    printf FP "# Slit function made by executing:\n";
    printf FP "#    $cmd\n";
    close(FP);
    my $test_t_start = time;
    my $ret = system("$cmd >> $out");
    my $test_time = time - $test_t_start;
    if ($ret != 0) {carp "make_slitfunction failed"; return 1;}
    return report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
}

sub uvspec_test_sbdart_thermal {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );
    printtest($test,$rtol,$atol);
    system("cat  $exinp > $inp");
    my $test_t_start = time;
    my $current_dir = getcwd; chdir($inpdir);
    my $ret = system("$uvspec_abs_path -f $inp_abs_path > $current_dir/TMPA.OUT");
    chdir($current_dir); # change back to dir where we started
    my $test_time = time - $test_t_start;
    print_randomseed($ret, $inpdir);
    if ($ret != 0) {carp "uvspec failed"; return 1;}

    open(FP1,"TMPA.OUT");
    open(FP2,">$out");
    while (my $line1 = <FP1>) {
      $line1 =~ s/^[ ]+//;
      my ($wvn1,$dir1,$dn1,$up1, $avg) = split(/ +/,$line1);
      printf FP2 "%f %f\n", $wvn1, $dn1*1.e+3;
    }
    close(FP1);
    close(FP2);

    #    system("gawk '{print \$1,\$3*1e3}' TMPA.OUT > $out");  # awk treats $3 as integer on iMac, do not known why
    # Replaced by above perl stuff by Arve

    `rm -f TMPA.OUT`;

    my $nox = "";
    my $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    print_randomseed($diff, $inpdir);
    return $ret, $diff;
}


sub run_uvspec {
  my ($inpfile, $outfile, $valgrind, $binary_path) = @_;
  if (not $binary_path) { $binary_path = $UVTools::uvspec; }
  if ( ! -e $binary_path ) { carp("\ncould not find binary: $binary_path\n"); return 1; }
  if ( ! -e $inpfile     ) { carp("\ncould not find input file: $inpfile\n"); return 1; }
  if($verbose) {
    my $current_dir = getcwd;
    print("Running in folder $current_dir");
  }
  my $test_t_start = time;
  my $cmd;
  if ($valgrind == 0) {
    $cmd = "$binary_path -f $inpfile > $outfile";
  } else {
    $cmd = "valgrind -q --leak-check=no --error-exitcode=9 $binary_path -f $inpfile > $outfile";
  }
  if($verbose) { print("\ncalling $cmd\n"); }
  my $ret = system($cmd);
  my $test_time = time - $test_t_start;
  if($verbose) { print("Returned $ret (time=$test_time)\n"); }
  return ($ret, $test_time);
}

sub run_uvspec_mpi {
  my ($inpfile, $outfile, $valgrind, $binary_path, $ncpu, $options) = @_;
  if (not $binary_path) { $binary_path = $UVTools::uvspec_mpi; }
  if ( ! -e $binary_path ) { carp("\ncould not find binary: $binary_path\n"); return 1; }
  if ( ! -e $inpfile     ) { carp("\ncould not find input file: $inpfile\n"); return 1; }
  if($verbose) {
    my $current_dir = getcwd;
    print("Running in folder $current_dir");
  }
  my $test_t_start = time;
  my $cmd;
  if ($valgrind == 0) {
    $cmd = "mpirun -np $ncpu $binary_path $inpfile $options > $outfile";
  } else {
    $cmd = "mpirun -np $ncpu valgrind -q --leak-check=no --error-exitcode=9 $binary_path $inpfile $options > $outfile";
  }
  if($verbose) { print("\ncalling $cmd\n"); }
  my $ret = system($cmd);
  my $test_time = time - $test_t_start;
  if($verbose) { print("Returned $ret (time=$test_time)\n"); }
  return ($ret, $test_time);
}


# Running uvspec and comparing the default output
sub uvspec_run_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier);
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, $out_abs_path, $valgrind, $uvspec_abs_path);
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed :: $ret"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
}

# Running uvspec and comparing the default output, but only the first 4 columns
sub uvspec_run_flx_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier);
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy.out", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1 , \$2 , \$3 , \$4}' dummy.out > $out_abs_path; rm dummy.out");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed :: $ret"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
}


sub apply_mystic_photon_multiplier {
  my ($inpfile, $multiplier, $force_present) = @_;
  my $current_photons = -1; #100000; # mystic default value

  # search mc_photons line and parse value
  open my $fh, '<', $inpfile or die "\nCould not open input file :: $!\n";
  while (my $line = <$fh>) { if ($line =~ /^\s*mc_photons\s+([\+]?[\d.]+(?:[eE]?[-\+]?\d+)?)/) { $current_photons = $1; } }
  close($fh) or die "\nCouldn't close file properly\n";

  if ($current_photons < 0) {
    if ($force_present) {
      die "\nCould not find mc_photons in test file $inpfile -- however, every montecarlo test has to set it.\n";
    } else {
      return;
    }
  }

  my $photons = $current_photons * $multiplier;

  if ($verbose) {
    print("\nIncreasing mc_photons for inpfile $inpfile with factor $multiplier\n");
    print("mc_photons was $current_photons is now $photons\n");
  }

  # append mc_photons line
  open my $fh, '>>', $inpfile or die "could not open input file for appending:: $!\n";
  print $fh "\nmc_photons $photons\n";
  close($fh) or die "Couldn't close file properly";
  return $photons;
}

sub mystic_run_test_mc_flx_spc {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' mc.flx.spc > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }


sub mystic_run_test_mc_rad {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$3,\$4,\$8}' mc.rad > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_bw_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind, $col) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1, \$$col}' mc.flx > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_t_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind, $col) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1, \$$col}' mc.abs.spc > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub uvspec_run_3d_abs_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind, $col) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    system("cat  $exinp > $inp");
    printtest($test,$rtol,$atol);
    # change to directory where input file resides,
    # running there with local paths allows for deeper dir structures in tests
    my $current_dir = getcwd; chdir($inpdir);

    ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
    system("gawk '{print \$1, \$$col}' mc.abs.spc > $out_abs_path; rm dummy");
    chdir($current_dir); # change back to dir where we started

    if ($ret != 0) {carp "uvspec failed"; return 1;}

    my $nox = "";
    $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    print_randomseed($diff, $inpdir);
    return $ret, $diff;
  }

sub mystic_run_irr_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk 'FNR==1 {print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' mc.flx > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_cb_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      printtest($test,$rtol,$atol);
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk 'FNR==2 {print \$1,\$4}' mc.mish.cb > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      my $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_rad_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$2,\$3,\$5}' mc.rad.spc > $out_abs_path; rm -f dummy mc.flx.spc mc.jac.spc mc.rad mc.rad.spc");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_flx_abs_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("paste mc.flx.spc mc.abs.spc | gawk '{print \$4,\$5,\$6,\$7,\$15}' > $out_abs_path"); # print iz, edir, edn, eup, absorption
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }




sub mystic_run_triangle_flx_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $atol, $rtol) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      # change to directory where input file resides,
      # running there with local paths allows for deeper dir structures in tests
      my $current_dir = getcwd; chdir($inpdir);

      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("rm -f dummy mc.flx.spc mc.jac.spc mc.rad mc.rad.spc mc1.flx mc.flx");
      chdir($current_dir); # change back to dir where we started

      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $mc_basename = "mc";
      open my $fh, '<', $inp or die "\nCould not open input file :: $!\n";
      while (my $line = <$fh>) { if ($line =~ /^\s*mc_basename\s+([^\s]+)/) { $mc_basename = $1; }}
      close($fh) or die "\nCouldn't close file properly\n";
      if ($verbose && ($mc_basename ne "mc")) { print("found mc_basename in inp file: $mc_basename\n"); }

      my $test_output = File::Spec->catfile($inpdir, $mc_basename.'.triangle.nc');
      my @variable_names = ('edir', 'edn', 'eup');
      $diff = nc_check_output(\@variable_names, $atol, $rtol, $test_output, $exout, $answers, $test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff
  }


sub uvspec_mpi_run_spc_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $ncpu, $options) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec_mpi );

    my $ret = -1;
    system("cat  $exinp > $inp");
    if ($verbose_uvspec) { `echo verbose >> $inp`; }
    apply_mystic_photon_multiplier($inp, $mc_photons_multiplier, "force_present");
    printtest($test,$rtol,$atol);

    my $current_dir = getcwd; chdir($inpdir);
    ($ret, my $test_time) = run_uvspec_mpi($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path, $ncpu, $options);
    if ($ret != 0) {carp "call to uvspec_mpi failed"; open my $fh, "dummy"; print <$fh>;}
    system("paste mc.flx.spc mc.abs.spc | gawk '{print \$4,\$5+\$6,\$7,\$15}' > $out_abs_path && rm -f dummy"); # print iz, edir+edn, eup, absorption
    chdir($current_dir); # change back to dir where we started

    print_randomseed($ret, $inpdir);
    if ($ret != 0) {carp "call to uvspec_mpi failed"; return 1;}

    my $nox = "";
    my $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    print_randomseed($diff, $inpdir);
    return $ret, $diff;
  }

sub mystic_run_is_IQ_test {
  # test MYSTIC spectral importance sampling (I and Q components)
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1,\$4,\$5}' mc.is.spc > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_aeris_test {
    # test MYSTIC concentration importance sampling
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$2, \$6}' mc.flx.is.spc > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
}

sub mystic_run_amf_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$3,\$4}' mc.amf > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub mystic_run_jacobian_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("gawk '{print \$1,\$2,\$3,\$4,\$5,\$6}' mc0.jac > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub lidar_run_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      my $photons = apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("cat mc.loc mc.loc_t mc.jac_t > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret != 0) {carp "uvspec failed"; return 1;}

      # We are now replacing the number of photons which is in the default output
      # so that we can run the test case with any number of photons
      # Read output file handle and unlink (we will write to the same filename in a second)
      open my $FH_old, "<", $out_abs_path or die "Unable to open $out_abs_path for read exited $? $!";
      unlink $out_abs_path;
      # Parse and replace text in same file
      open my $FH, ">", $out_abs_path or die "Unable to open $out_abs_path for write exited $? $!";
      while(my $line = <$FH_old>) {
        $line =~ s/$photons/Nphotons_replaced_for_tests/g;
        print $FH $line;
      }

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub radar_run_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my $inp_abs_path    = File::Spec->rel2abs( $inp );
    my $out_abs_path    = File::Spec->rel2abs( $out );
    my $uvspec_abs_path = File::Spec->rel2abs( $UVTools::uvspec );

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      apply_mystic_photon_multiplier($inp, $iter * $iter * $mc_photons_multiplier, "force_present");
      printtest($test,$rtol,$atol);

      my $current_dir = getcwd; chdir($inpdir);
      ($ret, my $test_time) = run_uvspec($inp_abs_path, "dummy", $valgrind, $uvspec_abs_path);
      system("cat mc.loc_t > $out_abs_path; rm dummy");
      chdir($current_dir); # change back to dir where we started
      print_randomseed($ret, $inpdir);
      if ($ret) { carp "uvspec failed"; return 1; }

      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff, $inpdir);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

sub ssradar_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);

    my ($ret, $diff) = -1, -1;
    for my $iter (1..$ITER_MC_TESTS) {
      system("cat  $exinp > $inp");
      printtest($test,$rtol,$atol);

      my $test_t_start = time;
      if ($valgrind == 0) { $ret = system("../bin/ssradar $inp > dummy"); }
      else { $ret = system("valgrind -q  --leak-check=no ../bin/ssradar $inp > dummy"); }
      my $test_time = time - $test_t_start;
      system("gawk 'FNR==9,FNR==14 {print \$3,\$4}' dummy > $out; rm dummy");
      print_randomseed($ret);
      if ($ret != 0) {carp "ssradar failed"; return 1;}
      my $nox = "";
      $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
      print_randomseed($diff);
      last if ($diff == 0 or $ret);
    }
    return $ret, $diff;
  }

# Running mie tool and comparing the default output
sub mie_run_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol, $valgrind) = @_;
    my ($inpvol, $inpdir, $inpfile) = File::Spec->splitpath($exinp);
    my $mie_abs_path = File::Spec->rel2abs( $UVTools::mie );
    system("cat  $exinp > $inp");
    printtest($test,$rtol,$atol);
    my ($ret, $test_time) = run_uvspec($inp, $out, $valgrind, $mie_abs_path);
    print_randomseed($ret);
    if ($ret != 0) {carp "mie failed :: $ret"; return 1;}
    my $nox = "";
    my $diff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
    print_randomseed($diff);
    return $ret, $diff;
}

sub calc_j_test {
    my $answers = shift @_;
    my $rtol    = shift @_;
    my $atol    = shift @_;
    my $test    = "j values";
    printtest($test,$rtol,$atol);
    my $inp = "TMP.INP";
    my $out = "TMP.OUT";
    my $exinp  = "../examples/J_UVSPEC.INP";
    my $exout  = "../examples/J.OUT";
    my $diff = calc_j_run_test($inp,$out,$exinp,$exout,$test,$answers,$rtol,$atol);
    `rm -f TMP_J_UVSPEC.OUT randomseed`;
    return $diff;
}

sub calc_j_run_test {
    my ($inp, $out, $exinp, $exout, $test, $answers, $rtol, $atol) = @_;
    my $diffout= "TMP_DIFF.OUT";
    system("cat  $exinp > $inp");
    my $uvspec_out = "./TMP_J_UVSPEC.OUT";

    ($ret, my $test_time) = run_uvspec($inp, $uvspec_out, $valgrind);
    if ($ret != 0) {carp "uvspec failed";return 1;}
    my $temp = 300.0;
    my $string = "'self_surface_temperature(K) $temp'";
    system("echo $string >> $uvspec_out");
    $ret = system("perl ../src/Calc_J.pl $uvspec_out > $out");
    if ($ret != 0) {carp "Calc_J.pl failed"; return  1;}
    my $nox = "--nox";
    return $ret + report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol,$test_time);
}

sub get_disort_mxd_value {
  my $f   = shift @_;
  my $key = shift @_;
  open(FPINC,$f) || carp "Could not open file $f";
  my @ls = <FPINC>;
  close(FPINC);

  my @keys    = grep (/$key/,@ls);
  (@keys)    = split(/\s+/,join(' ',@keys));
  my @keys   = grep (/[0-9]/,@keys);
  my $keyval = $keys[0];
  $keyval    =~ s/[,)]//g;
  $keyval;
}

sub airmass_test_bro {
  # This test uses the input density matrix option of sdisort. It requires
  # larger values of MXCLY and MXSZA in ../src_f/DISORT.MXD than the default
  # values. Hence we check and run the test if enough space is available.
  # Otherwise it is skipped.
    my $answers = shift @_;
    my $rtol    = shift @_;
    my $atol    = shift @_;
    my $uvspec_failure =0;
    my $test = "BrO airmass and slant column";
    printtest($test,$rtol,$atol);
    # Check if test can be run
    my $inc = "../src_f/DISORT.MXD";
    my $mxsza = get_disort_mxd_value($inc," MXSZA =");
    my $mxcly = get_disort_mxd_value($inc," MXCLY =");

    my $modify_disort_mxd=0;
    if ($mxsza < 48 || $mxcly < 480 ) {
      $modify_disort_mxd=1;
      my $ret = system("mv ../src_f/DISORT.MXD ./DISORT.MXD");
      open(FP,">../src_f/DISORT.MXD") || carp "Could not open file ../src_f/DISORT.MXD";
      printf FP "      PARAMETER ( MXCLY = 480,\n";
      printf FP "     &            MXULV = 2*MXCLY, \n";
      printf FP "     &            MXCMU = 16, \n";
      printf FP "     &            MXUMU = 2,\n";
      printf FP "     &            MXPHI = 2,\n";
      printf FP "     &            MXSZA = 48 )\n";
      close(FP);
      chdir("..");
#      $ret = system("make -s uvspec &> tmp_comp");
      `export MAKE='make -s'; make -s VERBOSE=0 uvspec &> tmp_comp`;
      chdir("test");
    }
    my $inp = "";
    my $out = "TMP.OUT";
    chdir("../src") || carp "Could not change directory to ../src";
    $ret = system("perl ../src/airmass.pl --amfX=bro --wvl=352 --sza_file=../examples/AMF_sza_am_1a.dat > $out");
   if ($ret != 0) {carp "airmass.pl failed"; return  1;}
    my $exout  = "../examples/AMF_BRO.OUT";
    my $nox = "--nox";
    my $ndiff = report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol);
    `rm -f $out tmp_*`;
    chdir("../test") || carp "Could not change directory to ../test";

    if (  $modify_disort_mxd ) {
      my $ret = system("mv ./DISORT.MXD ../src_f/DISORT.MXD ");
      my $ret = system("touch ../src_f/DISORT.MXD ");
      chdir("..");
      $ret = system("make -s VERBOSE=0 uvspec &> tmp_comp");
      chdir("test");
      $modify_disort_mxd = 0;
    }

    $uvspec_failure = $ndiff;
    $uvspec_failure;
}

sub airmass_test_hcho {
  # If you want to run this example for larger solar zenith angles the
  # number of vertical layers should be increased. A plausible vertical
  # resolution is 0.25km. This is easily achived with the altitude option.
  # (But also requires recompiling if MXSZA and MXCLY must be increased)
  my $answers = shift @_;
  my $rtol    = shift @_;
  my $atol    = shift @_;
  my $uvspec_failure =0;
  my $test = "HCHO airmass and slant column";
  printtest($test,$rtol,$atol);

  my $inc = "../src_f/DISORT.MXD";
  my $mxsza = get_disort_mxd_value($inc," MXSZA =");
  my $mxcly = get_disort_mxd_value($inc," MXCLY =");

  if ($mxsza < 48 || $mxcly < 126 ) {
    printf stderr " skipped.\n";
    printf stderr " MXSZA = %d, MXCLY = %d.\n", $mxsza, $mxcly;
    printf stderr " Increase MXSZA to 48 and MXCLY to 126 in $inc to run this test.\n";
    $uvspec_failure = 0;
  }
  else {

    my $inp = "";
    my $out = "TMP.OUT";

    # Just to make sure so we now where we are
    chdir("../test") || carp "Could not change directory to ../test";

    # cp standard airmass.pl to tmp_airmass.pl
    my $tmpamfpl = 'tmp_airmass.pl';
    `cp ../src/airmass.pl $tmpamfpl`;

    # modify tmp_airmass.pl to run test case
    open(FP,$tmpamfpl) || carp "Could not open file $tmpamfpl";
    my @ls = <FP>;
    close(FP);

    # Remove all the stuff we later want to modify
    @ls = grep(!/^dens_file/,@ls);
    @ls = grep(!/^crs_file/,@ls);
    #@ls = grep(!/^dens_column/,@ls);
    @ls = grep(!/^atmosphere_file/,@ls);

    # Add lines for this test:

    #  $ls[$#ls+1] = "dens_file O3     ../examples/AMF_o3_pro_1a_uvspec.dat\n";
    #  $ls[$#ls+1] = "dens_file NO2    ../examples/AMF_no2_pro_1a_uvspec.dat\n";
    $ls[$#ls+1] = "dens_file HCHO   ../examples/AMF_hcho_pro_uvspec.dat\n";
    $ls[$#ls+1] = "crs_file O3      ../examples/AMF_o3_025.xs\n";
    $ls[$#ls+1] = "crs_file NO2     ../examples/AMF_no2_025.xs\n";
    $ls[$#ls+1] = "crs_file HCHO    ../examples/AMF_hcho.xs\n";
    $ls[$#ls+1] = "dens_column HCHO 2.0e+17 CM\_2\n";

    $ls[$#ls+1] = "atmosphere_file  ../examples/AMF_atmmod.dat\n";
    $ls[$#ls+1] = "altitude 1.14\n";

    open(FP,">$tmpamfpl") || carp "Could not open file $tmpamfpl";
    foreach my $l (@ls) {  printf FP "%s", $l;}
    close(FP);

    if (-e "../src/$tmpamfpl") {`rm ../src/$tmpamfpl`; }
    `cp $tmpamfpl ../src/$tmpamfpl`;
    chdir("../src") || carp "Could not change directory to ../src";
    $ret = system("perl ../src/$tmpamfpl --amfX=hcho --wvl=356 --sza0=65 --sza1=95 --dsza=1 > $out");
    if ($ret != 0) {carp "airmass.pl failed"; return  1;}
    my $exout  = "../examples/AMF_HCHO.OUT";
    my $nox = "--nox";

    return report_differences($answers,$out,$exout,$inp,$test,$nox,$rtol,$atol);
    `rm -f $out tmp_*`;
    chdir("../test") || carp "Could not change directory to ../test";
    $uvspec_failure;
  }
  $uvspec_failure;

}

sub report_differences {
  my $answers = shift @_;
  my $out     = shift @_;
  my $exout   = shift @_;
  my $inp     = shift @_;
  my $test    = shift @_;
  my $nox     = shift @_;
  my $rtol    = shift @_;
  my $atol    = shift @_;
  my $time    = shift @_;
  my $diffout = "TMP_DIFF.OUT";

  if ( $answers == 1 ) {
    my $ans;
    if ($force_answers) {
      $ans = "y\n";
    } else {
      printf stderr "Replace $exout? [y/n]\n"; $ans = <STDIN>;
    }
    if ( $ans eq "y\n" ) {
      system("mv $out $exout");
      if ( -e $inp ) {system("rm -f $inp");}
      printf stderr " $exout replaced\n";
    }
    else {
      printf stderr " $exout NOT replaced\n";
    }
    return 0;
  }
  else {
    if (!-e $exout ) {croak "Example output file $exout not found\n";}
    if (!-e $out ) {croak "Test  output file $out not found\n";}
    if ($verbose) { print("Comparing $exout to $out\n"); }
    my $ret = system("perl ../src/ndiff $nox --quiet --limit $atol --maxdiff $rtol $exout $out > $diffout");
    if ($ret != 0) { carp "ndiff failed"; return 1; }
    open(FILE,$diffout);
    my $line = <FILE>;
    close(FILE);
    chop $line;
    my ($ndiffs,$largest, $absval0, $absval1) = split(/ +/,$line);
    if($time) { printf stderr " (%6.1f s) ", $time; }
    if ($ndiffs==0) {
      if ($verbose) {
        printf stderr " %4d serious differences. Maximum difference: %12.6e (%10.6f\%). Absolute values (ans,test): (%12.6e, %12.6e).\n",
        $ndiffs, $absval0-$absval1, $largest*100, $absval0, $absval1;
      } else {
        printf stderr " ok.\n";
      }
    } else {
      printf stderr " %4d serious differences. Maximum difference: %12.6e (%10.6f\%). Absolute values (ans,test): (%12.6e, %12.6e).\n",
      $ndiffs, $absval0-$absval1, $largest*100, $absval0, $absval1;
    }
    if ($ndiffs==0) {
      system("rm -f $diffout $inp $out");
    }
    return $ndiffs;
  }
}

sub report_results{
  my $answers = shift @_;
  my $failure = shift @_;;
  my $test    = shift @_;
  if ( $answers != 1 ) {
      printf stderr "$failure of the $test checks failed.\n";
  }
  return $failure;
}

sub printtest {
  my $test = shift @_;
  my $rtol = shift @_;
  my $atol = shift @_;
  printf stderr "$test (%7.5f, %5.1f%%)", $atol, $rtol*100; # Multiply rtol with 100 to show percentage.
  my $maxlength = 82;
  my $pnts = $maxlength-length($test);
  for (my $i=0;$i<$pnts;$i++) {
    printf stderr ".";
  }
}

# run nco commands to generate a new netcdf file with error metrics and give the results back as string (from ncdump)
sub nc_compare_variable_nco {
  my ($varname, $atol, $rtol, $inp, $target) = @_;
  my $cmd = <<"END";
    rm -f diff_tmp.nc;
     -v $varname --no_tmp_fl $inp $target -O -o diff_tmp.nc;
     -v --no_tmp_fl -s '${varname}_difference=${varname}' diff_tmp.nc -O -o diff.nc;
     --no_tmp_fl -A $inp diff.nc;
     -v --no_tmp_fl -s '
${varname}_max_err=${varname}_difference.fabs().max();
${varname}_max_rel_err=(${varname}_difference/(1e-300+${varname})).fabs().max();
${varname}_mean_err=${varname}_difference.avg();
${varname}_rel_err=(${varname}_difference/(1e-300+${varname})).fabs().max();
${varname}_mean=${varname}.avg();
${varname}_count_abs_err=(${varname}_difference.fabs() > ${atol}).total();
${varname}_count_rel_err=((${varname}_difference/(1e-300+${varname})).fabs() > ${rtol}).total();
' diff.nc -O -o diff_tmp.nc;
END
  if ($verbose) { printf stderr "Running netcdf comparison with nco command: \n%s\n", $cmd; }
  my $ret = system($cmd);
  if ($ret != 0) { return $ret; }
  my $diff_output = `ncdump diff_tmp.nc`;
  unlink "diff.nc";
  unlink "diff_tmp.nc";
  return $diff_output;
}

# check if the error metrics for $varname exceed absolute or relative errors and report the differences
my $nc_check_variable_limits_break_line = 0;
sub nc_check_variable_limits_nco {
   my ($varname, $atol, $rtol, $inp, $target) = @_;

   my $diff_output = nc_compare_variable_nco($varname, $atol, $rtol, $inp, $target);

   my ($var_count_abs  ) = ($diff_output =~ m/${varname}_count_abs_err\s+=\s+(\d+)\s+;/);
   my ($var_count_rel  ) = ($diff_output =~ m/${varname}_count_rel_err\s+=\s+(\d+)\s+;/);
   my ($var_mean_err   ) = ($diff_output =~ m/${varname}_mean_err\s+=\s+(\d+\.\d+)\s+;/);
   my ($var_rel_err    ) = ($diff_output =~ m/${varname}_rel_err\s+=\s+(\d+\.\d+)\s+;/);
   my ($var_max_err    ) = ($diff_output =~ m/${varname}_max_err\s+=\s+(\d+\.\d+)\s+;/);
   my ($var_max_rel_err) = ($diff_output =~ m/${varname}_max_rel_err\s+=\s+(\d+\.\d+)\s+;/);
   my ($var_mean       ) = ($diff_output =~ m/${varname}_mean\s+=\s+(\d+\.\d+)\s+;/);

   if ($var_count_abs) {
     if ($nc_check_variable_limits_break_line==0) { $nc_check_variable_limits_break_line = 1; printf "\n"; }
     printf "Variable %-8s exceeded the absolute tolerance (%8g):".
            "\tcount=%8d".
            "\tmean_err=%8g".
            "\trel_err=%8g%%".
            "\tmax_err=%8g".
            "\tmax_rel_err=%8g%%".
            "\t(variable_mean=%8g)".
            "\n",
            $varname,
            $atol,
            $var_count_abs,
            $var_mean_err,
            $var_rel_err*100,
            $var_max_err,
            $var_max_rel_err*100,
            $var_mean;
   } elsif ($var_count_rel) {
     if ($nc_check_variable_limits_break_line==0) { $nc_check_variable_limits_break_line = 1; printf "\n"; }
     printf "Variable %-8s exceeded the relative tolerance (%8g):".
            "\tcount=%8d".
            "\tmean_err=%8g".
            "\trel_err=%8g%%".
            "\tmax_err=%8g".
            "\tmax_rel_err=%8g%%".
            "\t(variable_mean=%8g)".
            "\n",
            $varname,
            $rtol,
            $var_count_rel,
            $var_mean_err,
            $var_rel_err*100,
            $var_max_err,
            $var_max_rel_err*100,
            $var_mean;
   }
   return $var_count_abs + $var_count_rel;
}

sub nc_check_variable_limits_nccmp {
   my ($varname, $atol, $rtol, $inp, $target) = @_;

   my $rtol_percent = $rtol * 100.;
   my $cmd = " -dmgf --tolerance=$atol --Tolerance=$rtol_percent --variable=$varname $inp $target";
   if ($verbose) { printf stderr "Running netcdf comparison with nccmp command: \n%s\n", $cmd; }
   my $ret = system($cmd);
   return $ret;
}

# check a netcdf output file for a number of variables
sub nc_check_output() {
    my ($variable_names_ref, $atol, $rtol, $test_output, $target, $answers, $time) = @_;
    my @variable_names = @$variable_names_ref;

    if ( ! -e $test_output ) { carp("expected test output does not exist: $test_output\n"); return 1; }
    if ( ! -e $target )      { carp("expected target output does not exist: $target\n"); return 2; }

    if ($answers) {
      my $ans;
        if ($force_answers) {
          $ans = "y\n";
        } else {
          printf stderr "Replace $target? [y/n]\n"; $ans = <STDIN>;
        }
        if ( $ans eq "y\n" ) {
            copy($test_output, $target);
            printf stderr " $target replaced\n";
        } else {
            printf stderr " $target NOT replaced\n";
        }
    }

    my $ret = 0;
    if($found_nco) {
        foreach my $varname (@variable_names) {
            $ret += nc_check_variable_limits_nco($varname, $atol, $rtol, $test_output, $target);
        }
    }
    if($found_nccmp) {
        foreach my $varname (@variable_names) {
          $ret += nc_check_variable_limits_nccmp($varname, $atol, $rtol, $test_output, $target);
        }
    }

    if ( -e $test_output ) { unlink $test_output; }
    if($time) { printf stderr " (%6.1f s) ", $time; }
    if(! ($found_nco || $found_nccmp)) {
      print(" skipping netcdf comparison\n");
      return $ret;
    }
    if ( $ret == 0 ) { print(" ok."); }
    print("\n");

    return $ret;
}

# in case that the $ret value != 0, try to find the randomseed file and report its contents
sub print_randomseed {
  my ($ret, $inpdir) = @_;
  if ($ret != 0) {
    if(!$inpdir) { $inpdir = getcwd; }
    my $randomseedfile = File::Spec->catfile($inpdir, 'randomseed');
    if (open (my $fh, '<', $randomseedfile)) {
      my $r = <$fh>;
      close $fh;
      chomp($r);
      print("RandomSeed = $r\n");
      unlink $randomseedfile;
    } else {
      #print("Could not open randomseed file: $randomseedfile\n");
    }
  }
}
# Local Variables:
# mode: Perl
# End:
