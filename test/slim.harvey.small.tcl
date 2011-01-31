
#
#  Slim Fast Input Script
# Cape Cod Bacterial Injection Experiment from
# Harvey and Garabedian (1991); file detailed in example section of 
# SLIMFAST manual.
#

file mkdir "slim.harvey.small"
cd "slim.harvey.small"

set n_slim_runs 1

#
#  Loop through runs
#
for {set k 1} {$k <= $n_slim_runs} {incr k 1} {

#
#	run id and other info
#
set parameter_file "slim.par"
set run_name  "slim_harvey_small EXAMPLE $k"
set log_file  "slog.harvey_small.EX.$k.txt"
#
#	input files for head, perm and porosity
#
set v_type calc
#set kx_file "../harvey_flow.$k.out.perm_x.pfb"
#set ky_file "../harvey_flow.$k.out.perm_x.pfb"
#set kz_file "../harvey_flow.$k.out.perm_x.pfb"
#set saturated "yes"
set kx_file "/home/cui1/weltyc_common/zcui/parflow_tcl/harvey_flow_small/harvey_flow.$k.out.perm_x.pfb"
set ky_file "/home/cui1/weltyc_common/zcui/parflow_tcl/harvey_flow_small/harvey_flow.$k.out.perm_x.pfb"
set kz_file "/home/cui1/weltyc_common/zcui/parflow_tcl/harvey_flow_small/harvey_flow.$k.out.perm_x.pfb"

## @RMM added to make consistent w/ var sat
set vga_file "../harvey_flow.vga.pfb"
set vgn_file "../harvey_flow.vgn.pfb"
set sres_file "../harvey_flow.sres.pfb"

set press 0
set head_file "/home/cui1/weltyc_common/zcui/parflow_tcl/harvey_flow_small/harvey_flow.$k.head.pfb"
set phi_type "constant"
set phi 0.39

set saturated "yes"

#
#	grid and domain info
#

set nx 50
set ny 8
set nz 100 
set dx .34 
set dy .34 
set dz .038


# particle number and numerics
set npmax 7500000
set give_up 100000
set epsilon 1E-10

#
# number of constituents
#
set num_constituents 11

set constituent_names(1) NH4
set constituent_names(2) NO3
set constituent_names(3) CH2O
set constituent_names(4) O2
set constituent_names(5) CO2
set constituent_names(6) HCO3
set constituent_names(7) H
set constituent_names(8) CO3
set constituent_names(9) Ca2 
set constituent_names(10) N2 
set constituent_names(11) OH 


#
# radioactive decay for constituents
#


set half_life(1) 0.
set half_life(2) 0.
set half_life(3) 0.
set half_life(4) 0.
set half_life(5) 0.
set half_life(6) 0.
set half_life(7) 0.
set half_life(8) 0.
set half_life(9) 0.
set half_life(10) 0.
set half_life(11) 0.
#
#	Dispersion and Diffusion Parameters
#

set alpha_l 0.001
set alpha_t 0.0001
set diffusivity 0.0
#
#	set linear sorption parameters
#
set sorption_type(1) "constant"
set sorption_value(1) 1.0

set sorption_type(2) "constant"
set sorption_value(2) 1.0

set sorption_type(3) "constant"
set sorption_value(3) 1.0

set sorption_type(4) "constant"
set sorption_value(4) 1.0

set sorption_type(5) "constant"
set sorption_value(5) 1.0

set sorption_type(6) "constant"
set sorption_value(6) 1.0

set sorption_type(7) "constant"
set sorption_value(7) 1.0

set sorption_type(8) "constant"
set sorption_value(8) 1.0

set sorption_type(9) "constant"
set sorption_value(9) 1.0

set sorption_type(10) "constant"
set sorption_value(10) 1.0

set sorption_type(11) "constant"
set sorption_value(11) 1.0


#
#	set matrix diffusion or Att/det parameters
#
set attach_type(1) "constant"
set attachment_value(1)  .0462
# kinetic k_forward rate [1/d]
set detach_type(1) "constant"
set detachment_value(1) 0.0869
# kinetic k_reverse rate [1/d]

set attach_type(2) "constant"
set attachment_value(2) 0.
set detach_type(2) "constant"
set detachment_value(2) 0.00

set attach_type(3) "constant"
set attachment_value(3) 0.02
set detach_type(3) "constant"
set detachment_value(3) 0.01

set attach_type(4) "constant"
set attachment_value(4) 0.
set detach_type(4) "constant"
set detachment_value(4) 0.

set attach_type(5) "constant"
set attachment_value(5) 0.
set detach_type(5) "constant"
set detachment_value(5) 0.

set attach_type(6) "constant"
set attachment_value(6) 0.
set detach_type(6) "constant"
set detachment_value(6) 0.

set attach_type(7) "constant"
set attachment_value(7) 0.
set detach_type(7) "constant"
set detachment_value(7) 0.

set attach_type(8) "constant"
set attachment_value(8) 0.
set detach_type(8) "constant"
set detachment_value(8) 0.

set attach_type(9) "constant"
set attachment_value(9) 0.
set detach_type(9) "constant"
set detachment_value(9) 0.

set attach_type(10) "constant"
set attachment_value(10) 0.
set detach_type(10) "constant"
set detachment_value(10) 0.

set attach_type(11) "constant"
set attachment_value(11) 0.
set detach_type(11) "constant"
set detachment_value(11) 0.


#
#	Set Run Timing Info
#
set time_increment 1
set num_increments 10


set vel_nskip  0
#
#	SlimNG Output
#
set write_well_breakthrough "yes"
set well_overwrite "yes"
set well_out_file(1) "well.EX.$k.txt"
set well_out_file(2) "well.EX2.$k.txt"
set well_out_file(3) "well.EX3.$k.txt"
set well_out_file(4) "well.EX4.$k.txt"
set well_out_file(5) "well.EX5.$k.txt"
set well_out_file(6) "well.EX6.$k.txt"
set well_out_file(7) "well.EX7.$k.txt"
set well_out_file(8) "well.EX8.$k.txt"
set well_out_file(9) "well.EX9.$k.txt"
set well_out_file(10) "well.EX10.$k.txt"
set well_out_file(11) "well.EX11.$k.txt"


set well_breakthrough_dt 1.
set number_well_steps 40

set write_concentrations "vtk"
set concentration_filename "conc.$k"
set concentration_header(1) "NH4"
set concentration_header(2) "NO3"
set concentration_header(3) "CH2O"
set concentration_header(4) "O2"
set concentration_header(5) "CO2"
set concentration_header(6) "HCO3"
set concentration_header(7) "H"
set concentration_header(8) "CO3"
set concentration_header(9) "Ca"
set concentration_header(10) "N2"
set concentration_header(11) "OH"


set write_out_particles  "no"
set particle_out_file  "parttrace.txt"
set write_out_moments "no"
set moment_out_file  "moments.txt"
set temp_averaging "no"
#
#	mode of operation
#   expects "forward" or "backward"
set simulation_mode "forward"
#set part_split "yes"
set part_split "no"
set min_conc 1E-20
set temp_averaging "no"
#
#	well information
#
# note, Z=0, model coords, corresponds to 5.96 m elev above msl and 12.1 m bls
#
set number_wells  0
# monitoring port on M01 at 8.510 m bls
set well_x_location(1) 12.07
set well_y_location(1) 5.27
set well_screen_top(1) 1.9
set well_screen_bottom(1) 1.862
set well_pumping_rate(1) 0.000

# monitoring port on M01 at 9.089 m bls
set well_x_location(2) 12.07
set well_y_location(2) 5.27
set well_screen_top(2) 1.330
set well_screen_bottom(2) 1.292
set well_pumping_rate(2) 0.000

# monitoring port on M07 at 8.504 m bls
set well_x_location(3) 9.83
set well_y_location(3) 4.96
set well_screen_top(3) 1.9
set well_screen_bottom(3) 1.862
set well_pumping_rate(3) 0.000

# monitoring port on M07 at 9.083 m bls

set well_x_location(4) 9.83
set well_y_location(4) 4.96
set well_screen_top(4) 1.330
set well_screen_bottom(4) 1.292
set well_pumping_rate(4) 0.000

#
#	Plane Information
#
set number_planes  0
set plane_out_file(1) "plane.EX.$k.txt"
set plane_xyz(1)  "x"
set plane_location(1) 11.9

set plane_out_file(2) "plane.EX.$k.txt"
set plane_out_file(3) "plane.EX.$k.txt"
set plane_out_file(4) "plane.EX.$k.txt"
set plane_out_file(5) "plane.EX.$k.txt"
set plane_out_file(6) "plane.EX.$k.txt"
set plane_out_file(7) "plane.EX.$k.txt"
set plane_out_file(8) "plane.EX.$k.txt"
set plane_out_file(9) "plane.EX.$k.txt"
set plane_out_file(10) "plane.EX.$k.txt"
set plane_out_file(11) "plane.EX.$k.txt"

#
#	Slim Initial Condition
#
set slim_initial_condition_type(1) "pulse"
set slim_initial_condition_type(2) "pulse"
set slim_initial_condition_type(3) "pulse"
set slim_initial_condition_type(4) "pulse"
set slim_initial_condition_type(5) "pulse"
set slim_initial_condition_type(6) "pulse"
set slim_initial_condition_type(7) "pulse"
set slim_initial_condition_type(8) "pulse"
set slim_initial_condition_type(9) "pulse"
set slim_initial_condition_type(10) "pulse"
set slim_initial_condition_type(11) "pulse"

#NH4
set slim_pulseIC_Xlower(1) 0.68
set slim_pulseIC_Xupper(1) 1.02 
set slim_pulseIC_Ylower(1) 0.68 
set slim_pulseIC_Yupper(1) 2.04
set slim_pulseIC_Zlower(1) 0.076 
set slim_pulseIC_Zupper(1) 3.724


set slim_pulseIC_Initial_Concentration(1) 40.
#set slim_pulseIC_number_particles(1) 10
set slim_pulseIC_number_particles(1) 1000 
set slim_pulseIC_decay_rate(1) 0.00001
set slim_pulseIC_decay_time(1) .000
set slim_pulseIC_decay_timesteps(1) 1

#NO3
set slim_pulseIC_Xlower(2) 0.68
set slim_pulseIC_Xupper(2) 1.02 
set slim_pulseIC_Ylower(2) 0.68
set slim_pulseIC_Yupper(2) 2.04 
set slim_pulseIC_Zlower(2) 0.076
set slim_pulseIC_Zupper(2) 3.724
set slim_pulseIC_Initial_Concentration(2) 5.76 
#set slim_pulseIC_number_particles(2) 10
set slim_pulseIC_number_particles(2) 1000
set slim_pulseIC_decay_rate(2) 0.00001
set slim_pulseIC_decay_time(2) .000
set slim_pulseIC_decay_timesteps(2) 1

# CH2O

set slim_pulseIC_Xlower(3) 0.68
set slim_pulseIC_Xupper(3) 1.02 
set slim_pulseIC_Ylower(3) 0.68 
set slim_pulseIC_Yupper(3) 2.04 
set slim_pulseIC_Zlower(3) 0.076 
set slim_pulseIC_Zupper(3) 3.724
set slim_pulseIC_Initial_Concentration(3) 82.0
#set slim_pulseIC_number_particles(3) 10
set slim_pulseIC_number_particles(3) 1000
set slim_pulseIC_decay_rate(3) 0.00001
set slim_pulseIC_decay_time(3) .000
set slim_pulseIC_decay_timesteps(3) 1

# O2

set slim_pulseIC_Xlower(4) 0.68
set slim_pulseIC_Xupper(4) 1.02
set slim_pulseIC_Ylower(4) 0.68 
set slim_pulseIC_Yupper(4) 2.04 
set slim_pulseIC_Zlower(4) 0.076 
set slim_pulseIC_Zupper(4) 3.724
set slim_pulseIC_Initial_Concentration(4) 6.0
#set slim_pulseIC_number_particles(4) 10
set slim_pulseIC_number_particles(4) 500 
set slim_pulseIC_decay_rate(4) 0.00001
set slim_pulseIC_decay_time(4) .000
set slim_pulseIC_decay_timesteps(4) 1

# CO2

set slim_pulseIC_Xlower(5) 0.68
set slim_pulseIC_Xupper(5) 1.02
set slim_pulseIC_Ylower(5) 0.68 
set slim_pulseIC_Yupper(5) 2.04 
set slim_pulseIC_Zlower(5) 0.076 
set slim_pulseIC_Zupper(5) 3.724
set slim_pulseIC_Initial_Concentration(5) 57.2
#set slim_pulseIC_number_particles(5) 10
set slim_pulseIC_number_particles(5) 1000 
set slim_pulseIC_decay_rate(5) 0.00001
set slim_pulseIC_decay_time(5) .000
set slim_pulseIC_decay_timesteps(5) 1

# HCO3

set slim_pulseIC_Xlower(6) 0.68
set slim_pulseIC_Xupper(6) 1.02
set slim_pulseIC_Ylower(6) 0.68 
set slim_pulseIC_Yupper(6) 2.04 
set slim_pulseIC_Zlower(6) 0.076 
set slim_pulseIC_Zupper(6) 3.724
set slim_pulseIC_Initial_Concentration(6) 445.3
#set slim_pulseIC_number_particles(6) 10
set slim_pulseIC_number_particles(6) 5000 
set slim_pulseIC_decay_rate(6) 0.00001
set slim_pulseIC_decay_time(6) .000
set slim_pulseIC_decay_timesteps(6) 1

# H

set slim_pulseIC_Xlower(7) 0.68
set slim_pulseIC_Xupper(7) 1.02 
set slim_pulseIC_Ylower(7) 0.68 
set slim_pulseIC_Yupper(7) 2.04 
set slim_pulseIC_Zlower(7) 0.076 
set slim_pulseIC_Zupper(7) 3.724
set slim_pulseIC_Initial_Concentration(7) 6.8e-5
#set slim_pulseIC_number_particles(7) 1
set slim_pulseIC_number_particles(7) 100
set slim_pulseIC_decay_rate(7) 0.00001
set slim_pulseIC_decay_time(7) .000
set slim_pulseIC_decay_timesteps(7) 1

# CO3

set slim_pulseIC_Xlower(8) 0.68
set slim_pulseIC_Xupper(8) 1.02 
set slim_pulseIC_Ylower(8) 0.68 
set slim_pulseIC_Yupper(8) 2.04 
set slim_pulseIC_Zlower(8) 0.076 
set slim_pulseIC_Zupper(8) 3.724
set slim_pulseIC_Initial_Concentration(8) 0.24
#set slim_pulseIC_number_particles(8) 1
set slim_pulseIC_number_particles(8) 800
set slim_pulseIC_decay_rate(8) 0.00001
set slim_pulseIC_decay_time(8) .000
set slim_pulseIC_decay_timesteps(8) 1

# Ca2 

set slim_pulseIC_Xlower(9) 0.68
set slim_pulseIC_Xupper(9) 1.02 
set slim_pulseIC_Ylower(9) 0.68 
set slim_pulseIC_Yupper(9) 2.04 
set slim_pulseIC_Zlower(9) 0.076 
set slim_pulseIC_Zupper(9) 3.724
set slim_pulseIC_Initial_Concentration(9) 37.2
#set slim_pulseIC_number_particles(9) 10
set slim_pulseIC_number_particles(9) 1000
set slim_pulseIC_decay_rate(9) 0.00001
set slim_pulseIC_decay_time(9) .000
set slim_pulseIC_decay_timesteps(9) 1

# N2 
set slim_pulseIC_Xlower(10) 0.68
set slim_pulseIC_Xupper(10) 1.02 
set slim_pulseIC_Ylower(10) 0.68 
set slim_pulseIC_Yupper(10) 2.04 
set slim_pulseIC_Zlower(10) 0.076 
set slim_pulseIC_Zupper(10) 3.724
set slim_pulseIC_Initial_Concentration(10) 0.0
set slim_pulseIC_number_particles(10) 0
set slim_pulseIC_decay_rate(10) 0.00001
set slim_pulseIC_decay_time(10) .000
set slim_pulseIC_decay_timesteps(10) 1

# OH 
set slim_pulseIC_Xlower(11) 0.68
set slim_pulseIC_Xupper(11) 1.02 
set slim_pulseIC_Ylower(11) 0.68 
set slim_pulseIC_Yupper(11) 2.04 
set slim_pulseIC_Zlower(11) 0.076 
set slim_pulseIC_Zupper(11) 3.724
set slim_pulseIC_Initial_Concentration(11) 1.122e-3
#set slim_pulseIC_number_particles(11) 1
set slim_pulseIC_number_particles(11) 100
set slim_pulseIC_decay_rate(11) 0.00001
set slim_pulseIC_decay_time(11) .000
set slim_pulseIC_decay_timesteps(11) 1
#Nitrogen transport and transformation

set npar_file "slim.npar"

#mg/l of porous media
set biomass1_concentration_mg_l 0.1
#mg/l of porous media
set biomass2_concentration_mg_l 0.1
#microbial yield coefficient for biomass 1 (M biomass/ M substrate)
set Y1_M_biomass_M_substrate 0.17 
#microbial yield coefficient for biomass 2 (M biomass/ M substrate)
set Y2_M_biomass_M_substrate 0.5
#emperical biomass 1 inhibition constant (M biomass/L^3 porous medium) mg/l
set Kb1_mg_l 1.0 
#emperical biomass 2 inhibition constant (M biomass/L^3 porous medium) mg/l
set Kb2_mg_l 0.5 
#CH2O half-saturation constant (M species/L^3 water) mg/l
set K_CH2O_mg_l 10 
#O2 half-saturation constant (M species/L^3 water) mg/l
set K_O2_mg_l 0.1
#NH4 half-saturation constant (M species/L^3 water) mg/l
set K_NH4_mg_l 0.1
 #NO3 half-saturation constant (M species/L^3 water) mg/l
set K_NO3_mg_l 0.5
#O2 inhibition coefficient(M species/L^3 water) mg/l
set K_O2I_mg_l 1.0 
#CH2O distribution coefficient cm^3/g
set Kd_CH2O_cm3_g  1.5 
#NH4 distribution coefficient cm^3/g
set Kd_NH4_cm3_g  0.34 
# specific biomass decay or maintenance constant (1/T) 1/day for biomass 1
set k1d_1_day 0.0 
# specific biomass decay or maintenance constant (1/T) 1/day for biomass 2
set k2d_1_day 0.0 
 #O2 Henry's constant dimensionless
set H_O2  28.2 
#N2 Henry's constant dimensionless
set H_N2  56.6  
 #CO2 Henry's constant dimensionless
set H_CO2 0.95  
#substrate free air diffusion coefficient m^2/day
set Dg1_m2_day  1.0  
 #Electron acceptor free air diffusion coefficient m^2/day
set Dg2_m2_day  5.0 
 # organic carbon oxidation maximum primary substrate utilization rate 1/day
set k_max_ox_1_day 10.0
 # denitrification maximum primary substrate utilization rate 1/day
set k_max_denit_1_day 10.0
# nitrification maximum primary substrate utilization rate 1/day
set k_max_nit_1_day 1.0 
 # CO2 forward rate 1/day
#Chen's number
set k_f_co2_1_day 2.5920
#MacQuarrie's number
#set k_f_co2_1_day 2592.0
# CO2 backward rate 1/(M*day)
#Chen's number
set k_b_co2_1_mday 6.818E6 
#MacQuarrie's number
#set k_b_co2_1_mday 6.818E9 
# CO2+OH=HCO3 forward rate 1/(M*day)
set k_f_H_1_mday 7.344E8 
# CO2+H=HCO3 backward rate 1/day
set k_b_H_1_day 8.63 
# CO3 +H2O = HCO3 +H forward rate 1/day
#Chen's number
set k_f_HCO3_1_day 1.08E4 
#MacQuarrie's number
#set k_f_HCO3_1_day 1.11E11 
# CO3 +H2O = HCO3 +H backward rate 1/(M*day)
#Chen's number
set k_b_HCO3_1_mday 9.216E7 
#MacQuarrie's number
#set k_b_HCO3_1_mday 9.22E14 
# H2O=H+OH forward rate 1/day
#Chen's number
set k_f_H20_1_day 1.1E-4 
#MacQuarrie's number
#set k_f_H20_1_day 110.0 
# H2O=H+OH backward rate 1/(M*day)
#Chen's number
set k_b_H2O_1_mday 2.46E10 
#MacQuarrie's number
#set k_b_H2O_1_mday 2.46E16 
#CaCO3(s) + H = HCO3+Ca2+ forward rate mol/(cm^2*day)
set k_f_1_cm_day 7.69    
 #CaCO3(s) + H = HCO3+Ca2+ backward rate
set k_b_1_cm_mday 7.69E-2   
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca forward rate mol(cm^2*day)
set k_f_2_cm_day 4.32E-3 
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca backward rate mol(cm^2*day)
#Chen's number
set k_b_2_cm_m2day 1.14E4 
#MacQuarrie's number
#set k_b_2_cm_m2day 1.14E2 
#CaCO3(s) = CO3 + Ca forward rate mol(cm^2*day)
#Chen's number
set k_f_3_mol_cm2day 6.0E-6 
#MacQuarrie's number
#set k_f_3_mol_cm2day 5.62E-6 
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca backward rate mol(cm^2*day)
#Chen's number
set k_b_3_cm_mday 1.512E3 
#MacQuarrie's number
#set k_b_3_cm_mday 1.512E-3 
# reactive mineral surface area per unit vol. of porous media estimated  1/cm
# 3.8 is estimated value
set Acc_1_cm 3.8 

#NH4
set slim_background_conc(1) 0.0
#NO3
set slim_background_conc(2) 0.0
#CH2O
set slim_background_conc(3) 0.0
#O2
set slim_background_conc(4) 6.0
# CO2
set slim_background_conc(5) 23.32
# HCO3
set slim_background_conc(6) 207.4
# H
set slim_background_conc(7) 5.9E-5
# CO3
set slim_background_conc(8) 0.132
# Ca2 
set slim_background_conc(9) 68.0
# N2 
set slim_background_conc(10) 0.0
# OH 
set slim_background_conc(11) 1.292e-3


source $env(SLIM_DIR)/tcl/slim.run.tcl

}

cd ".."
