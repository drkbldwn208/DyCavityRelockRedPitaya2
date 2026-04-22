create_project dy_cavity_relocker_2 . -part xc7z010clg400-1 -force

set_property ip_repo_paths {../hls_workspace/dy_cavity_relocker_2/solution1/impl/ip} [current_project]
update_ip_catalog

add_files -norecurse ../../extern/red-pitaya-notes/cores/axis_red_pitaya_adc.v
add_files -norecurse ../../src/axis_red_pitaya_dac.v

add_files -fileset constrs_1 -norecurse ../../src/constraints.xdc

update_compile_order -fileset sources_1

create_bd_design "system"

  # Create interface ports
  set DDR [ create_bd_intf_port -mode Master -vlnv xilinx.com:interface:ddrx_rtl:1.0 DDR ]

  set FIXED_IO [ create_bd_intf_port -mode Master -vlnv xilinx.com:display_processing_system7:fixedio_rtl:1.0 FIXED_IO ]


  # Create ports
  set adc_enc_p_o [ create_bd_port -dir O adc_enc_p_o ]
  set adc_clk_n_i [ create_bd_port -dir I adc_clk_n_i ]
  set adc_enc_n_o [ create_bd_port -dir O adc_enc_n_o ]
  set adc_csn_o [ create_bd_port -dir O adc_csn_o ]
  set adc_dat_a_i [ create_bd_port -dir I -from 15 -to 0 adc_dat_a_i ]
  set adc_dat_b_i [ create_bd_port -dir I -from 15 -to 0 adc_dat_b_i ]
  set adc_clk_p_i [ create_bd_port -dir I adc_clk_p_i ]
  set dac_clk_o [ create_bd_port -dir O dac_clk_o ]
  set dac_rst_o [ create_bd_port -dir O dac_rst_o ]
  set dac_sel_o [ create_bd_port -dir O dac_sel_o ]
  set dac_wrt_o [ create_bd_port -dir O dac_wrt_o ]
  set dac_dat_o [ create_bd_port -dir O -from 13 -to 0 dac_dat_o ]


  # Create instance: processing_system7_0, and set properties
  set processing_system7_0 [ create_bd_cell -type ip -vlnv xilinx.com:ip:processing_system7:5.5 processing_system7_0 ]
  set_property CONFIG.PCW_FPGA_FCLK0_ENABLE {1} $processing_system7_0


  # Create instance: clk_wiz_0, and set properties
  set clk_wiz_0 [ create_bd_cell -type ip -vlnv xilinx.com:ip:clk_wiz:6.0 clk_wiz_0 ]
  set_property -dict [list \
    CONFIG.CLKIN1_JITTER_PS {80.0} \
    CONFIG.CLKOUT1_JITTER {119.348} \
    CONFIG.CLKOUT1_PHASE_ERROR {96.948} \
    CONFIG.CLKOUT1_REQUESTED_OUT_FREQ {125} \
    CONFIG.CLKOUT2_JITTER {119.348} \
    CONFIG.CLKOUT2_PHASE_ERROR {96.948} \
    CONFIG.CLKOUT2_REQUESTED_OUT_FREQ {125} \
    CONFIG.CLKOUT2_USED {true} \
    CONFIG.CLKOUT3_JITTER {119.348} \
    CONFIG.CLKOUT3_PHASE_ERROR {96.948} \
    CONFIG.CLKOUT3_REQUESTED_OUT_FREQ {125} \
    CONFIG.CLKOUT3_USED {true} \
    CONFIG.CLKOUT4_JITTER {104.759} \
    CONFIG.CLKOUT4_PHASE_ERROR {96.948} \
    CONFIG.CLKOUT4_REQUESTED_OUT_FREQ {250} \
    CONFIG.CLKOUT4_USED {true} \
    CONFIG.MMCM_CLKFBOUT_MULT_F {8.000} \
    CONFIG.MMCM_CLKIN1_PERIOD {8.000} \
    CONFIG.MMCM_CLKOUT0_DIVIDE_F {8.000} \
    CONFIG.MMCM_CLKOUT1_DIVIDE {8} \
    CONFIG.MMCM_CLKOUT2_DIVIDE {8} \
    CONFIG.MMCM_CLKOUT3_DIVIDE {4} \
    CONFIG.NUM_OUT_CLKS {4} \
    CONFIG.PRIM_IN_FREQ {125} \
    CONFIG.PRIM_SOURCE {Differential_clock_capable_pin} \
    CONFIG.USE_RESET {false} \
  ] $clk_wiz_0


  # Create instance: axis_red_pitaya_adc_0, and set properties
  set block_name axis_red_pitaya_adc
  set block_cell_name axis_red_pitaya_adc_0
  if { [catch {set axis_red_pitaya_adc_0 [create_bd_cell -type module -reference $block_name $block_cell_name] } errmsg] } {
     catch {common::send_gid_msg -ssname BD::TCL -id 2095 -severity "ERROR" "Unable to add referenced block <$block_name>. Please add the files for ${block_name}'s definition into the project."}
     return 1
   } elseif { $axis_red_pitaya_adc_0 eq "" } {
     catch {common::send_gid_msg -ssname BD::TCL -id 2096 -severity "ERROR" "Unable to referenced block <$block_name>. Please add the files for ${block_name}'s definition into the project."}
     return 1
   }
  
  # Create instance: proc_sys_reset_0, and set properties
  set proc_sys_reset_0 [ create_bd_cell -type ip -vlnv xilinx.com:ip:proc_sys_reset:5.0 proc_sys_reset_0 ]

  # Create instance: dy_cavity_relocker_2_0, and set properties
  set dy_cavity_relocker_2_0 [ create_bd_cell -type ip -vlnv xilinx.com:hls:dy_cavity_relocker_2:1.0 dy_cavity_relocker_2_0 ]

  # Create instance: axis_red_pitaya_dac_0, and set properties
  set block_name axis_red_pitaya_dac
  set block_cell_name axis_red_pitaya_dac_0
  if { [catch {set axis_red_pitaya_dac_0 [create_bd_cell -type module -reference $block_name $block_cell_name] } errmsg] } {
     catch {common::send_gid_msg -ssname BD::TCL -id 2095 -severity "ERROR" "Unable to add referenced block <$block_name>. Please add the files for ${block_name}'s definition into the project."}
     return 1
   } elseif { $axis_red_pitaya_dac_0 eq "" } {
     catch {common::send_gid_msg -ssname BD::TCL -id 2096 -severity "ERROR" "Unable to referenced block <$block_name>. Please add the files for ${block_name}'s definition into the project."}
     return 1
   }
  
  # Create instance: axis_clock_converter_0, and set properties
  set axis_clock_converter_0 [ create_bd_cell -type ip -vlnv xilinx.com:ip:axis_clock_converter:1.1 axis_clock_converter_0 ]
  set_property -dict [list \
    CONFIG.IS_ACLK_ASYNC {0} \
    CONFIG.TDATA_NUM_BYTES {4} \
  ] $axis_clock_converter_0


  # Create instance: proc_sys_reset_1, and set properties
  set proc_sys_reset_1 [ create_bd_cell -type ip -vlnv xilinx.com:ip:proc_sys_reset:5.0 proc_sys_reset_1 ]

  # Create instance: ps7_0_axi_periph, and set properties
  set ps7_0_axi_periph [ create_bd_cell -type ip -vlnv xilinx.com:ip:axi_interconnect:2.1 ps7_0_axi_periph ]
  set_property CONFIG.NUM_MI {1} $ps7_0_axi_periph


  # Create interface connections
  connect_bd_intf_net -intf_net axis_clock_converter_0_M_AXIS [get_bd_intf_pins axis_clock_converter_0/M_AXIS] [get_bd_intf_pins axis_red_pitaya_dac_0/s_axis]
  connect_bd_intf_net -intf_net axis_red_pitaya_adc_0_m_axis [get_bd_intf_pins axis_red_pitaya_adc_0/m_axis] [get_bd_intf_pins dy_cavity_relocker_2_0/adc_in]
  connect_bd_intf_net -intf_net dy_cavity_relocker_2_0_dac_out [get_bd_intf_pins dy_cavity_relocker_2_0/dac_out] [get_bd_intf_pins axis_clock_converter_0/S_AXIS]
  connect_bd_intf_net -intf_net processing_system7_0_DDR [get_bd_intf_ports DDR] [get_bd_intf_pins processing_system7_0/DDR]
  connect_bd_intf_net -intf_net processing_system7_0_FIXED_IO [get_bd_intf_ports FIXED_IO] [get_bd_intf_pins processing_system7_0/FIXED_IO]
  connect_bd_intf_net -intf_net processing_system7_0_M_AXI_GP0 [get_bd_intf_pins processing_system7_0/M_AXI_GP0] [get_bd_intf_pins ps7_0_axi_periph/S00_AXI]
  connect_bd_intf_net -intf_net ps7_0_axi_periph_M00_AXI [get_bd_intf_pins ps7_0_axi_periph/M00_AXI] [get_bd_intf_pins dy_cavity_relocker_2_0/s_axi_control]

  # Create port connections
  connect_bd_net -net adc_clk_n_i_1 [get_bd_ports adc_clk_n_i] [get_bd_pins clk_wiz_0/clk_in1_n]
  connect_bd_net -net adc_clk_p_i_1 [get_bd_ports adc_clk_p_i] [get_bd_pins clk_wiz_0/clk_in1_p]
  connect_bd_net -net adc_dat_a_i_1 [get_bd_ports adc_dat_a_i] [get_bd_pins axis_red_pitaya_adc_0/adc_dat_a]
  connect_bd_net -net adc_dat_b_i_1 [get_bd_ports adc_dat_b_i] [get_bd_pins axis_red_pitaya_adc_0/adc_dat_b]
  connect_bd_net -net axis_red_pitaya_adc_0_adc_csn [get_bd_pins axis_red_pitaya_adc_0/adc_csn] [get_bd_ports adc_csn_o]
  connect_bd_net -net axis_red_pitaya_dac_0_dac_dat [get_bd_pins axis_red_pitaya_dac_0/dac_dat] [get_bd_ports dac_dat_o]
  connect_bd_net -net axis_red_pitaya_dac_0_dac_clk [get_bd_pins axis_red_pitaya_dac_0/dac_clk] [get_bd_ports dac_clk_o]
  connect_bd_net -net axis_red_pitaya_dac_0_dac_rst [get_bd_pins axis_red_pitaya_dac_0/dac_rst] [get_bd_ports dac_rst_o]
  connect_bd_net -net axis_red_pitaya_dac_0_dac_sel [get_bd_pins axis_red_pitaya_dac_0/dac_sel] [get_bd_ports dac_sel_o]
  connect_bd_net -net axis_red_pitaya_dac_0_dac_wrt [get_bd_pins axis_red_pitaya_dac_0/dac_wrt] [get_bd_ports dac_wrt_o]
  connect_bd_net -net clk_wiz_0_clk_out1 [get_bd_pins clk_wiz_0/clk_out1] [get_bd_pins processing_system7_0/M_AXI_GP0_ACLK] [get_bd_pins axis_red_pitaya_adc_0/aclk] [get_bd_pins proc_sys_reset_0/slowest_sync_clk] [get_bd_pins dy_cavity_relocker_2_0/ap_clk] [get_bd_pins axis_clock_converter_0/s_axis_aclk] [get_bd_pins ps7_0_axi_periph/ACLK] [get_bd_pins ps7_0_axi_periph/S00_ACLK] [get_bd_pins ps7_0_axi_periph/M00_ACLK]
  connect_bd_net -net clk_wiz_0_clk_out2 [get_bd_pins clk_wiz_0/clk_out2] [get_bd_pins axis_red_pitaya_dac_0/ddr_clk]
  connect_bd_net -net clk_wiz_0_clk_out3 [get_bd_pins clk_wiz_0/clk_out3] [get_bd_pins axis_red_pitaya_dac_0/wrt_clk]
  connect_bd_net -net clk_wiz_0_clk_out4 [get_bd_pins clk_wiz_0/clk_out4] [get_bd_pins axis_red_pitaya_dac_0/aclk] [get_bd_pins proc_sys_reset_1/slowest_sync_clk] [get_bd_pins axis_clock_converter_0/m_axis_aclk]
  connect_bd_net -net clk_wiz_0_locked [get_bd_pins clk_wiz_0/locked] [get_bd_pins axis_red_pitaya_dac_0/locked]
  connect_bd_net -net proc_sys_reset_0_peripheral_aresetn [get_bd_pins proc_sys_reset_0/peripheral_aresetn] [get_bd_pins axis_clock_converter_0/s_axis_aresetn] [get_bd_pins dy_cavity_relocker_2_0/ap_rst_n] [get_bd_pins ps7_0_axi_periph/S00_ARESETN] [get_bd_pins ps7_0_axi_periph/M00_ARESETN] [get_bd_pins ps7_0_axi_periph/ARESETN]
  connect_bd_net -net proc_sys_reset_1_peripheral_aresetn [get_bd_pins proc_sys_reset_1/peripheral_aresetn] [get_bd_pins axis_clock_converter_0/m_axis_aresetn]
  connect_bd_net -net processing_system7_0_FCLK_RESET0_N [get_bd_pins processing_system7_0/FCLK_RESET0_N] [get_bd_pins proc_sys_reset_0/ext_reset_in] [get_bd_pins proc_sys_reset_1/ext_reset_in]

  # Create address segments
  assign_bd_address -offset 0x40000000 -range 0x00010000 -target_address_space [get_bd_addr_spaces processing_system7_0/Data] [get_bd_addr_segs dy_cavity_relocker_2_0/s_axi_control/Reg] -force

validate_bd_design

make_wrapper -files [get_files ./dy_cavity_relocker_2.srcs/sources_1/bd/system/system.bd] -top
add_files -norecurse ./dy_cavity_relocker_2.srcs/sources_1/bd/system/hdl/system_wrapper.v

launch_runs impl_1 -to_step write_bitstream -jobs 10
wait_on_run impl_1

exit