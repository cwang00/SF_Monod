
#
#  Slim Fast Input Script
# Cape Cod Bacterial Injection Experiment from
# Harvey and Garabedian (1991); file detailed in example section of 
# SLIMFAST manual.
#

file mkdir "slim.harvey"
cd "slim.harvey"

set n_slim_runs 1

#
#  Loop through runs
#
for {set k 1} {$k <= $n_slim_runs} {incr k 1} {

#
#	run id and other info
#
set parameter_file "slim.par"
set run_name  "slim_harvey_bact EXAMPLE $k"
set log_file  "slog.harvey_bact.EX.$k.txt"
#
#	input files for head, perm and porosity
#
set v_type calc
set kx_file "../harvey_flow.$k.out.perm_x.pfb"
set ky_file "../harvey_flow.$k.out.perm_x.pfb"
set kz_file "../harvey_flow.$k.out.perm_x.pfb"
set saturated "yes"
## @RMM added to make consistent w/ var sat
set vga_file "../harvey_flow.vga.pfb"
set vgn_file "../harvey_flow.vgn.pfb"
set sres_file "../harvey_flow.sres.pfb"

set press 0
set head_file "../harvey_flow.$k.head.pfb"
set phi_type "constant"
set phi 0.39

set saturated "no"

#
#	grid and domain info
#

set nx 50
set ny 30
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
set num_constituents 2

set constituent_names(1) NH4
set constituent_names(2) NO3
set constituent_names(3) CH2O
set constituent_names(4) O2
set constituent_names(5) CO2
set constituent_names(6) HCO3
set constituent_names(7) Hplus
set constituent_names(8) CO3
set constituent_names(9) Ca2 
set constituent_names(10) N2 
set constituent_names(11) OH 


#
# radioactive decay for constituents
#


set half_life(1) 0.
set half_life(2) 0.
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


#
#	Set Run Timing Info
#
set time_increment 4.
set num_increments 10


set vel_nskip  0
#
#	SlimNG Output
#
set write_well_breakthrough "yes"
set well_overwrite "yes"
set well_out_file(1) "well.EX.$k.txt"
set well_out_file(2) "well.EX2.$k.txt"


set well_breakthrough_dt 1.
set number_well_steps 40

set write_concentrations "vtk"
set concentration_filename "conc.$k"
set concentration_header(1) "Metal"
set concentration_header(2) "Tracer"


set write_out_particles  "no"
set particle_out_file  "parttrace.txt"
set write_out_moments "no"
set moment_out_file  "moments.txt"
set temp_averaging "no"
#
#	mode of operation
#   expects "forward" or "backward"
set simulation_mode "forward"
set part_split "yes"
set min_conc 1E-20
set temp_averaging "no"
#
#	well information
#
# note, Z=0, model coords, corresponds to 5.96 m elev above msl and 12.1 m bls
#
set number_wells  2
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

#
#	Slim Initial Condition
#
set slim_initial_condition_type(1) "pulse"
set slim_initial_condition_type(2) "pulse"


set slim_pulseIC_Xlower(1) 4.76
set slim_pulseIC_Xupper(1) 5.25 
set slim_pulseIC_Ylower(1) 5.1 
set slim_pulseIC_Yupper(1) 5.44 
set slim_pulseIC_Zlower(1) 1.1 
set slim_pulseIC_Zupper(1) 2.1


set slim_pulseIC_Initial_Concentration(1) 1.
set slim_pulseIC_number_particles(1) 100000
set slim_pulseIC_decay_rate(1) 0.00001
set slim_pulseIC_decay_time(1) .000
set slim_pulseIC_decay_timesteps(1) 1


set slim_pulseIC_Xlower(2) 4.76
set slim_pulseIC_Xupper(2) 5.25 
set slim_pulseIC_Ylower(2) 5.1
set slim_pulseIC_Yupper(2) 5.44 
set slim_pulseIC_Zlower(2) 1.1
set slim_pulseIC_Zupper(2) 2.1
set slim_pulseIC_Initial_Concentration(2) 1.
set slim_pulseIC_number_particles(2) 100000
set slim_pulseIC_decay_rate(2) 0.00001
set slim_pulseIC_decay_time(2) .000
set slim_pulseIC_decay_timesteps(2) 1

#Nitrogen transport and transformation

set npar_file "slim.npar"

#mg/l of porous media
set biomass1_concentration 0.1 
#mg/l of porous media
set biomass2_concentration 0.1  
#microbial yield coefficient for biomass 1 (M biomass/ M substrate)
set Y1 0.17 
#microbial yield coefficient for biomass 2 (M biomass/ M substrate)
set Y2 0.5
#emperical biomass 1 inhibition constant (M biomass/L^3 porous medium) mg/l
set Kb1 1.0 
#emperical biomass 2 inhibition constant (M biomass/L^3 porous medium) mg/l
set Kb2 0.5 
#CH2O half-saturation constant (M species/L^3 water) mg/l
set K_CH2O 10 
#O2 half-saturation constant (M species/L^3 water) mg/l
set K_O2 0.1
#NH4 half-saturation constant (M species/L^3 water) mg/l
set K_NH4 0.1
 #NO3 half-saturation constant (M species/L^3 water) mg/l
set K_NO3 0.5
#O2 inhibition coefficient(M species/L^3 water) mg/l
set K_O2I 1.0 
#CH2O distribution coefficient cm^3/g
set Kd_CH2O  1.5 
#NH4 distribution coefficient cm^3/g
set Kd_NH4  0.34 
# specific biomass decay or maintenance constant (1/T) 1/day for biomass 1
set k1d 0.0 
# specific biomass decay or maintenance constant (1/T) 1/day for biomass 2
set k2d 0.0 
 #O2 Henry's constant
set H_O2  28.2 
#N2 Henry's constant
set H_N2  56.6  
 #CO2 Henry's constant
set H_CO2 0.95  
#substrate free air diffusion coefficient m^2/day
set Dg1  1.0  
 #Electron acceptor free air diffusion coefficient m^2/day
set Dg2  5.0 
 # organic carbon oxidation maximum primary substrate utilization rate 1/day
set k_max_ox 10.0
 # denitrification maximum primary substrate utilization rate 1/day
set k_max_denit 10.0
# nitrification maximum primary substrate utilization rate 1/day
set k_max_nit 1.0 
 # CO2 forward rate
set k_f_co2 2.5920
# CO2 backward rate 1/(M*day)
set k_b_co2 6.818E6 
# CO2+H=HCO3 forward rate 1/(M*day)
set k_f_H 7.344E8 
# CO2+H=HCO3 backward rate 1/day
set k_b_H 8.63 
# CO3 +H2O = HCO3 +H forward rate 1/day
set k_f_HCO3 1.08E4 
# CO3 +H2O = HCO3 +H backward rate 1/(M*day)
set k_b_HCO3 9.216E7 
# H2O=H+OH forward rate 1/day
set k_f_H20 1.1E-4 
# H2O=H+OH backward rate 1/(M*day)
set k_b_H2O 2.46E10 
#CaCO3(s) + H = HCO3+Ca forward rate mol/(cm^2*day)
set k_f_1 7.69    
 #CaCO3(s) + H = HCO3+Ca forward rate
set k_b_1 7.69E-2   
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca forward rate mol(cm^2*day)
set k_f_2 4.32E-3 
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca backward rate mol(cm^2*day)
set k_b_2 1.14E4 
#CaCO3(s) = CO3 + Ca forward rate mol(cm^2*day)
set k_f_3 6.0E-6 
#CaCO3(s) + H2O + CO2(aq) = 2HCO3 + Ca backward rate mol(cm^2*day)
set k_b_3 1.512E3 
# reactive mineral surface area per unit vol. of porous media estimated  1/cm
# 3.8 is estimated value
set Acc 3.8 



source $env(SLIM_DIR)/tcl/slim.run.tcl

}

cd ".."
