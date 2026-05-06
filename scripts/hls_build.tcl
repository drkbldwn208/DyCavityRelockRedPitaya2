open_project dy_cavity_relocker_2
set_top dy_cavity_relocker_2

set source_files [glob ../../src/*.cpp ../../src/*.h ../../src/*.hpp]

foreach file $source_files {
    add_files $file
}

open_solution "solution1" -flow_target vivado

set_part {xc7z010clg400-1}

create_clock -period 8 -name default

csynth_design
export_design -flow syn -rtl verilog -format ip_catalog
export_design -format ip_catalog -output ../ip_repo/relocker_2_ip.zip

exit