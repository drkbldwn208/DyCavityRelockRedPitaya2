# scripts/hls_csim.tcl — Run C simulation of dy_cavity_relocker_2
# Usage: cd build/sim_workspace && vitis_hls -f ../../scripts/hls_csim.tcl

open_project -reset dy_cavity_relocker_2_csim

# Add source files
add_files ../../src/dy_cavity_relocker_2.cpp
add_files ../../src/dy_cavity_relocker_2.h
add_files ../../src/hinf_filter.hpp
add_files ../../src/hinf_coeffs.h
add_files ../../src/fsm_receiver.cpp
add_files ../../src/fsm_receiver.h

# Add testbench
add_files -tb ../../tb/tb_freq_response.cpp

# Set top-level function
set_top dy_cavity_relocker_2

# Create solution
open_solution -reset solution1

# Set target device (same as your actual hardware)
set_part xc7z010clg400-1
create_clock -period 8 -name default

# Run C simulation
csim_design -clean

exit