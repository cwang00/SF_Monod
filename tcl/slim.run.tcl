#
# executable & I/O
#

set slim_executable $env(SLIM_DIR)/bin/SLIM.exe

#
#
# Write to the Parameter file:
#
set fileId [open $parameter_file w 0600]
set nfileId [open $npar_file w 0600]

puts $fileId "\"$run_name\"						! Title of Run"
puts $fileId "\"$log_file\"						! Log File Name"

if {$saturated == "yes"} {
   puts $fileId "1                                            !saturated"
 } elseif { $saturated == "const" } {
   puts $fileId "2                                          !const V. G. pars"
 } elseif { $saturated == "const_sat" } {
   puts $fileId  "3                            ! const saturation and V.G. pars"
 } else {
   puts $fileId "0                                             !saturated"
 }


if {$v_type == "vbin"} {
puts $fileId "1							! reading V's from vbin type file"
puts $fileId "\"$vel_file\"					! Vbin File Name" 
if {$saturated != "yes"} {
  if { $saturated == "const" } {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sres					 ! VG Sres"
  } elseif { $saturated == "const_sat"} {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sat				                 ! const sat"
  } else {
     puts $fileId "\"$vga_file\"								! VG alpha File Name" 
     puts $fileId "\"$vgn_file\"								! VG n File Name" 
     puts $fileId "\"$sres_file\"								! VG Sres File Name"
  } 
 }
 puts $fileId "$phi									! phi value \[-\]"
}


if {$v_type == "calc"} {
puts $fileId "2							! calcing V's internally, steady state"
puts $fileId "\"$kx_file\"								! Kx File Name" 
puts $fileId "\"$ky_file\"								! Ky File Name" 
puts $fileId "\"$kz_file\"								! Kz File Name" 
## @RMM 8-6-08 added to make consistent w/ varsat
if {$saturated != "yes"} {
  if { $saturated == "const" } {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sres					 ! VG Sres"
  } elseif { $saturated == "const_sat"} {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sat				                 ! const sat"
  } else {
     puts $fileId "\"$vga_file\"								! VG alpha File Name" 
     puts $fileId "\"$vgn_file\"								! VG n File Name" 
     puts $fileId "\"$sres_file\"								! VG Sres File Name"
  } 
}
puts $fileId "$press								! head (0) or pressure (1) flag" 
puts $fileId "\"$head_file\"								! Head File Name" 
if {$phi_type == "constant"} {
puts $fileId "1										! const phi"
puts $fileId "$phi									! phi value \[-\]"
}
if {$phi_type == "PFBFile"} {
puts $fileId "2										! const phi"
puts $fileId "\"$phi_file\"									! phi value \[-\]"
}
}
if {$v_type == "calc_trans"} {
puts $fileId "3							! calcing v's internally, transient"
puts $fileId "\"$kx_file\"								! Kx File Name" 
puts $fileId "\"$ky_file\"								! Ky File Name" 
puts $fileId "\"$kz_file\"								! Kz File Name" 
if {$saturated != "yes"} {
  if {$saturated == "const"} {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sres					 ! VG Sres"
  } elseif { $saturated == "const_sat"} {
     puts $fileId "$vga				                 ! VG alpha" 
     puts $fileId "$vgn					         ! VG n " 
     puts $fileId "$sat				                 ! const sat"
  } else {
    puts $fileId "\"$vga_file\"								! VG alpha File Name" 
    puts $fileId "\"$vgn_file\"								! VG n File Name" 
    puts $fileId "\"$sres_file\"								! VG Sres File Name" 
  }
}
puts $fileId "$press								! head (0) or pressure (1) flag" 
puts $fileId "\"$time_file\"								! File w/ timesteps " 
puts $fileId "\"$head_list_file\"								! File w/  pressure or head" 
if {$phi_type == "constant"} {
puts $fileId "1										! const phi"
puts $fileId "$phi									! phi value \[-\]"
}
if {$phi_type == "PFBFile"} {
puts $fileId "2										! const phi"
puts $fileId "\"$phi_file\"									! phi value \[-\]"
}
}

if {[info exists write_vel_field_gnuplot ]} {
	if { $write_vel_field_gnuplot == "yes" } {
           puts $fileId "1                        !if write velocity field to a gnuplot file, 1 - yes, 0 - no"

           puts $fileId "\"$vel_field_file\"                        !velocity field gnuplot file name"
	} else {
           puts $fileId "0                        !if write velocity field to a gnuplot file, 1 - yes, 0 - no"
           puts $fileId "\"*********\"                        !velocity field gnuplot file name"
	}
} else {
           puts $fileId "0                        !if write velocity field to a gnuplot file, 1 - yes, 0 - no"
           puts $fileId "\"*********\"                        !velocity field gnuplot file name"
}

puts $fileId "$nx 						! number of x-nodes in the domain \[-\]"
puts $fileId "$ny 						! number of y-nodes in the domain \[-\]"
puts $fileId "$nz 						! number of z-nodes in the domain \[-\]"
puts $fileId "$dx 						! delta-x \[m\]"
puts $fileId "$dy 						! delta-y \[m\]"
puts $fileId "$dz 						! delta-z \[m\]"
puts $fileId "$num_constituents 				! number of constituents"
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "$half_life($jj)				        ! radioactive half life,  \[d\]" 
}
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $nfileId "$jj   $constituent_names($jj)				        ! number to  constituent names,  \[d\]" 
}
puts $fileId "$alpha_l						! longitudinal disersivity, alpha_l, \[m\]"
puts $fileId "$alpha_t						! transverse disersivity, alpha_t, \[m\]"

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
if {$sorption_type($jj) == "constant"} {
puts $fileId "$sorption_value($jj)					! Linear Retardation Coeff, R, \[-\]" 
}
if {$sorption_type($jj) == "min_file"} {
puts $fileId "-9999					! min_file" 
puts $fileId "\"$mineralization_file($jj)\""
}
if {$sorption_type($jj) == "kd_file"} {
puts $fileId "-1					! kd_file" 
puts $fileId "\"$kd_file($jj)\""
puts $fileId "$bulk_density($jj) 				! bulk density"

}
if {$sorption_type($jj) == "R_file"} {
puts $fileId "-2					! R_file" 
puts $fileId "\"$R_file($jj)\""

}

}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
if {$attach_type($jj) == "constant"} {
puts $fileId "$attachment_value($jj)					! Attachment Coeff, Katt, \[1/d\]" 
puts $fileId "$detachment_value($jj)					! Detachment Coeff, Kdet, \[1/d\]" 
}
if {$attach_type($jj) == "file"} {
puts $fileId "-9999					! file" 
puts $fileId "-9999					! file" 
puts $fileId "\"$att_det_file($jj)\""
puts $fileId "\"$perm_catagory_file($jj)\""
}

#if {$detach_type($jj) == "constant"} {
#puts $fileId "$detachment_value($jj)					! Detachment Coeff, Kdet, \[1/d\]" 
#}
#if {$detach_type($jj) == "file"} {
#puts $fileId "-9999					! file" 
#}

}


puts $fileId "$time_increment					! time increment until reporting particle vales, Tnext \[d\]"
puts $fileId "$num_increments					! number of time increments, nt, (total runtime=nt*tnext in days)"

if {$write_well_breakthrough == "yes"} {if {$well_overwrite == "no"} {puts $fileId "1				! appending well BTC to file:"
} else { puts $fileId "2				! writing out well BTC to file:"}
} else { puts $fileId "0				! not writing out well BTC"}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$well_out_file($jj)\""
}

puts $fileId "$well_breakthrough_dt				! time discretization for wells \[days\]"
puts $fileId "$number_well_steps				! number of time discretizations for wells"
if {$write_concentrations == "yes"} {puts $fileId "1				! writing out concentrations using header:"
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$concentration_header($jj)\""
}
} elseif {$write_concentrations == "vtk"} {
puts $fileId "2				! writing out concentrations in vtk using header:"
puts $fileId "\"$concentration_filename\""
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$concentration_header($jj)\""
}
} elseif {$write_concentrations == "no"} { puts $fileId "0				! not writing out concentrations"}

if {$write_out_particles == "yes"} {puts $fileId "1				! writing out particle locations to file:"
} else { puts $fileId "0				! not writing out particle locations"}
puts $fileId "\"$particle_out_file\""

if {$write_out_moments == "yes"} {puts $fileId "1				! writing out moments to file:"
} else { puts $fileId "0				! not writing out moments"}
puts $fileId "\"$moment_out_file\""

if {$simulation_mode == "forward"} {puts $fileId "0				! forward simulation"
} else { puts $fileId "1				! backward simulation"}

puts $fileId "$diffusivity						! molucular diffusivity,  \[m**2/day\]"

puts $fileId "$npmax						! maximum number of particles"
puts $fileId "$give_up						! maximum number of particle steps per time loop"
puts $fileId "$epsilon						! epsilon/zero value drop tolerance"


puts $fileId "$number_wells						! number of wells \[x,y,zt,zb,q\]"
for {set ii 1} {$ii <= $number_wells} {incr ii 1} {
puts $fileId "$well_x_location($ii) , $well_y_location($ii) , $well_screen_top($ii) , $well_screen_bottom($ii) , $well_pumping_rate($ii) "
}

if {[info exists bdy_cond ]} {
   if { $bdy_cond == "const_conc" } {
       puts $fileId "const_conc                      !const conc boundary"
   } else {
       puts $fileId  "const_flux                    !const flux boundary" 
   }
} else {
puts $fileId  "const_flux                          !const flux boundary" 
}

puts $fileId "$modelname                           !select which reaction model"

if {[info exists parflow_mask ]} {
	puts $fileId "1                               ! has parFlow mask file"
	puts $fileId "\"$parflow_mask\"                 ! ParFlow mask file"
} else {
puts $fileId  "0                          ! does not have parFlow mask file"
puts $fileId "\"*********\"                ! ParFlow mask file place holder"
}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$plane_out_file($jj)\""
}


puts $fileId "$number_planes						! number of planes \[x/y/z,coord\]"

for {set ii 1} {$ii <= $number_planes} {incr ii 1} {
if {$plane_xyz($ii) == "x"} {puts $fileId "1				! plane in x"}
if {$plane_xyz($ii) == "y"} {puts $fileId "2				! plane in y"}
if {$plane_xyz($ii) == "z"} {puts $fileId "3				! plane in z"}
puts $fileId "$plane_location($ii)					! location of plane \[m\]"
}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {

if {$slim_initial_condition_type($jj) == "well"} {
puts $fileId "1							! Well-Type IC"
puts $fileId "$well_ic_x($jj)				! Well IC: X \[m\] "
puts $fileId "$well_ic_y($jj)				! 	Y \[m\]"
puts $fileId "$well_ic_bot($jj)				! Zlower of well screen \[m\]"
puts $fileId "$well_ic_top($jj)				! Zupper of well screen \[m\]"
puts $fileId "$well_ic_pump_rate($jj)				! pumping rate \[vol/day\] "
puts $fileId "$well_ic_num_part_per($jj)				! number of particles per node "
}

if {$slim_initial_condition_type($jj) == "pulse"} {
puts $fileId "2							! Pulse-Type IC"
puts $fileId "$slim_pulseIC_Xlower($jj)				! Pulse IC: Xlower \[m\] "
puts $fileId "$slim_pulseIC_Xupper($jj)                              ! Xupper \[m\] "
puts $fileId "$slim_pulseIC_Ylower($jj)				! Ylower \[m\]"
puts $fileId "$slim_pulseIC_Yupper($jj)                             ! Yupper \[m\]"
puts $fileId "$slim_pulseIC_Zlower($jj)				! Zlower \[m\]"
puts $fileId "$slim_pulseIC_Zupper($jj)				! Zupper \[m\]"
puts $fileId "$slim_pulseIC_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_pulseIC_number_particles($jj)				! num particles total "
puts $fileId "$slim_pulseIC_decay_rate($jj)				! pulse decay rate "
puts $fileId "$slim_pulseIC_decay_time($jj)				! pulse decay time \[days\] "
puts $fileId "$slim_pulseIC_decay_timesteps($jj)				! decay timesteps \[time/timesteps\] "
}
if {$slim_initial_condition_type($jj) == "ind_file"} {
puts $fileId "3							! Indicator file-Type IC"
puts $fileId "\"$slim_ic_ind_file($jj)\""
puts $fileId "$slim_ic_ind_cat($jj)						! indicator cat"
puts $fileId "$slim_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_number_particles($jj)				! num particles total "
}

if {$slim_initial_condition_type($jj) == "ind_mult"} {
puts $fileId "4							! Mulitple Indicator; reading in file-Type IC"
puts $fileId "\"$slim_ic_ind_file($jj)\""
}

if {$slim_initial_condition_type($jj) == "ind_mult_background"} {
puts $fileId "5							! Mulitple Indicator; reading in second-Type IC, use background conc to initialize"
puts $fileId "\"$slim_ic_ind_file($jj)\""
puts $fileId "\"$slim_ic_ind_file_bg($jj)\""
puts $fileId "$slim_ic_ind_cat($jj)						! indicator cat"
puts $fileId "$slim_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_number_particles($jj)				! num particles total "
}

if {$slim_initial_condition_type($jj) == "none"} {
puts $fileId "0							! No IC"
}

if {$slim_initial_condition_type($jj) == "conc_at_bnd"} {
puts $fileId "6							! given conc at boundary reading in file-Type IC"
puts $fileId "\"$slim_ic_ind_file($jj)\""
puts $fileId "\"$slim_ic_ind_file_bg($jj)\""
puts $fileId "$slim_ic_ind_cat($jj)"
puts $fileId "$slim_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_number_particles($jj)				! num particles total "
}
}

# particle splitting?
if {$part_split == "yes"} {puts $fileId "1							! Particle Splitting"
puts $fileId "$min_conc						! min conc"
} else { puts $fileId "0							! Fixed Number of Particles" }
# temporal averaging?
if {$temp_averaging == "yes"} {puts $fileId "1							! Temporal Averaging"
} else { puts $fileId "0							! Not temporally avging for concentration" }
puts $fileId "$vel_nskip						! number of timesteps to skip"
puts $fileId "$Xup_Ref                               ! X upper reflection?"
puts $fileId "$bnd_Xup                               !X boundary location"
puts $fileId "$Xdown_Ref                               ! X down reflection?"
puts $fileId "$bnd_Xdown                               !X boundary location"
puts $fileId "$Yup_Ref                               ! Y upper reflection?"
puts $fileId "$bnd_Yup                               !Y boundary location"
puts $fileId "$Ydown_Ref                               ! Y down reflection?"
puts $fileId "$bnd_Ydown                               !Y boundary location"
puts $fileId "$Zup_Ref                               ! Z upper reflection?"
puts $fileId "$bnd_Zup                               !Z boundary location"
puts $fileId "$Zdown_Ref                               ! Z donw reflection?"
puts $fileId "$bnd_Zdown                               !Z boundary location"
puts $fileId " "
puts $fileId " "

puts $fileId " "
puts $fileId " "

if {$modelname == "MacQ1990"} {
  puts $nfileId "$kmax      ! maximum rate of substrate utilization"
  puts $nfileId "$Ks        ! substrate half-saturation constant"
  puts $nfileId "$Ka        ! electron acceptor half-saturation constant"
  puts $nfileId "$Rs        ! substrate retardation"
  puts $nfileId "$Rm        ! electron acceptor retardation"
  puts $nfileId "$Y         ! biomass yield coefficient"
  puts $nfileId "$r         ! Utilizaton ratio, X"
  puts $nfileId "$decay     ! biomass decay coefficient"
} elseif { $modelname == "MacQ1990unsat" } {
  puts $nfileId "$kmax      ! maximum rate of substrate utilization"
  puts $nfileId "$Ks        ! substrate half-saturation constant"
  puts $nfileId "$Ka        ! electron acceptor half-saturation constant"
  puts $nfileId "$Rs        ! substrate retardation"
  puts $nfileId "$Rm        ! electron acceptor retardation"
  puts $nfileId "$Y         ! biomass yield coefficient"
  puts $nfileId "$r1         ! Utilizaton ratio, X"
  puts $nfileId "$r2         ! Utilizaton ratio, X"
  puts $nfileId "$decay     ! biomass decay coefficient"
  puts $nfileId "$Kb        ! biomass inhibition coefficient"
  puts $nfileId "$SH        ! substrate Henry's constant"
  puts $nfileId "$EH        ! electron acceptor Henry's constant"
  puts $nfileId "$SD0        ! substract free air diffusion coefficient"
  puts $nfileId "$ED0        ! electorn acceptor free air diffusion coefficient"
} elseif { $modelname == "MacQ"} {
puts $nfileId "$biomass1_concentration_mg_l ! biomass1 concentration mg/l of porous media"
puts $nfileId "$biomass2_concentration_mg_l ! biomass2 concentration mg/l of porous media"
puts $nfileId "$Y1_M_biomass_M_substrate !microbial yield coefficient for biomass 1 (M biomass/ M substrate)"
puts $nfileId "$Y2_M_biomass_M_substrate !microbial yield coefficient for biomass 2 (M biomass/ M substrate)"
puts $nfileId "$Kb1_mg_l !emperical biomass 1 inhibition constant (M biomass/L^3 porous medium) mg/l"
puts $nfileId "$Kb2_mg_l !emperical biomass 2 inhibition constant (M biomass/L^3 porous medium) mg/l"
puts $nfileId "$K_CH2O_mg_l !CH2O half-saturation constant (M species/L^3 water) mg/l"
puts $nfileId "$K_O2_mg_l !O2 half-saturation constant (M species/L^3 water) mg/l"
puts $nfileId "$K_NH4_mg_l !NH4 half-saturation constant (M species/L^3 water) mg/l"
puts $nfileId "$K_NO3_mg_l !NO3 half-saturation constant (M species/L^3 water) mg/l"
puts $nfileId "$K_O2I_mg_l !O2 inhibition coefficient(M species/L^3 water) mg/l"
puts $nfileId "$Kd_CH2O_cm3_g !CH2O distribution coefficient cm^3/g"
puts $nfileId "$Kd_NH4_cm3_g !NH4 distribution coefficient cm^3/g"
puts $nfileId "$k1d_1_day ! specific biomass decay or maintenance constant (1/T) 1/day for biomass 1"
puts $nfileId "$k2d_1_day ! specific biomass decay or maintenance constant (1/T) 1/day for biomass 2"
puts $nfileId "$H_O2 !O2 Henry's constant"
puts $nfileId "$H_N2 !N2 Henry's constant"
puts $nfileId "$H_CO2 !CO2 Henry's constant"
puts $nfileId "$Dg1_m2_day !substrate free air diffusion coefficient m^2/day"
puts $nfileId "$Dg2_m2_day !Electron acceptor free air diffusion coefficient m^2/day"
puts $nfileId "$k_max_ox_1_day ! organic carbon oxidation maximum primary substrate utilization rate 1/day"
puts $nfileId "$k_max_denit_1_day ! denitrification maximum primary substrate utilization rate 1/day"
puts $nfileId "$k_max_nit_1_day ! nitrification maximum primary substrate utilization rate 1/day"
puts $nfileId "$k_f_co2_1_day ! CO2 forward rate 1/day" 
puts $nfileId "$k_b_co2_1_mday ! CO2 backward rate 1/(M*day)"
puts $nfileId "$k_f_H_1_mday ! CO2+H=HCO3 forward rate 1/(M*day)"
puts $nfileId "$k_b_H_1_day ! CO2+H=HCO3 backward rate 1/day"
puts $nfileId "$k_f_HCO3_1_day ! CO3 +H2O = HCO3 +H forward rate 1/day"
puts $nfileId "$k_b_HCO3_1_mday ! CO3 +H2O = HCO3 +H backward rate 1/(M*day)"
puts $nfileId "$k_f_H20_1_day ! H2O=H+OH forward rate 1/day"
puts $nfileId "$k_b_H2O_1_mday ! H2O=H+OH backward rate 1/(M*day)"
puts $nfileId "$k_f_1_cm_day !CaCO3(s) + H = HCO3+Ca forward rate mol/(cm^2*day)"
puts $nfileId "$k_b_1_cm_mday !CaCO3(s) + H = HCO3+Ca forward rate"
puts $nfileId "$k_f_2_cm_day !CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca forward rate mol(cm^2*day)"
puts $nfileId "$k_b_2_cm_m2day !CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca backward rate mol(cm^2*day)"
puts $nfileId "$k_f_3_mol_cm2day !CaCO3(s) = CO3 + Ca forward rate mol(cm^2*day)"
puts $nfileId "$k_b_3_cm_mday !CaCO3(s) = CO3 + Ca backward rate mol(cm^2*day)"
puts $nfileId "$Acc_1_cm ! reactive mineral surface area per unit vol. of porous media estimated  1/cm"
puts $nfileId "$soil_bulk_den_g_cm3 ! soil bulk density cm^3/g"

} elseif {$modelname == "VG"} {

puts $nfileId "$firstorderdecay_1_day ! first-order decay"
puts $nfileId "$zeroorderdecay_1_day ! zero-order decay"

} elseif {$modelname == "Chen1992"} {
 puts $nfileId "$kmax_tol ! maximum rate of substrate utilization"
 puts $nfileId "$kmax_ben ! maximum rate of substrate utilization"
 puts $nfileId   "$Ktol   ! toluene half-saturation constant"
 puts $nfileId    "$Kben   ! benzene half-saturation constant"
 puts $nfileId    "$Kdo   ! electron acceptor half-saturation constant"
 puts $nfileId    "$Rtol   ! substrate retardation"
 puts $nfileId    "$Rben   ! substrate retardation"
 puts $nfileId    "$Rm   ! electron acceptor retardation"
 puts $nfileId    "$Y    ! biomass yield coefficient"
 puts $nfileId    "$r_tol_do    ! Utilizaton ratio, X" 
 puts $nfileId    "$r_ben_do    ! Utilizaton ratio, X"
 puts $nfileId    "$decay ! biomass decay coefficient"
}
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
   puts $nfileId "$slim_background_conc_type($jj)       !backgorund type, const or PFBfile"
   if { $slim_background_conc_type($jj) == "const" } {
     puts $nfileId "$slim_background_conc($jj)            ! background conc"
   } else {
     puts $nfileId "\"$slim_background_conc($jj)\"          ! background conc"
   }
}


close $fileId
close $nfileId
#
#  Run SLIM 
#
set fileId [open "slim_dummy_file.par" w 0600] 
puts $fileId "$parameter_file"
puts $fileId "$npar_file"
close $fileId

#exec $slim_executable  < slim_dummy_file.par > slim.out.txt

