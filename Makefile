SHELL := /bin/bash

SRC_DIR = src
SCRIPT_DIR = scripts
EXTERN_DIR = extern
BUILD_DIR = build

SOURCE_FILES = $(wildcard $(SRC_DIR)/*.cpp)
TESTBENCH = $(wildcard tb/*.cpp)

SIM_DIR = $(BUILD_DIR)/sim_workspace
HLS_DIR = $(BUILD_DIR)/hls_workspace
IP_REPO_DIR = $(BUILD_DIR)/ip_repo
VIVADO_DIR = $(BUILD_DIR)/vivado_workspace
MONITOR_SCRIPT = $(SCRIPT_DIR)/move_bit.sh
MONITOR_LOG = $(BUILD_DIR)/monitor.log

# Path to Vitis HLS csim output (where CSVs land)
CSIM_BUILD_DIR = $(SIM_DIR)/dy_cavity_relocker_2_csim/solution1/csim/build

.PHONY: all csim csim_plot hls vivado external_cores clean start_monitor stop_monitor

all: external_cores hls vivado

# Run C simulation through Vitis HLS (required for ap_fixed / HLS types)
csim:
	@echo "Running C simulation via Vitis HLS..."
	@mkdir -p $(SIM_DIR)
	cd $(SIM_DIR) && \
	source /home/levlabcukomen/tools/Vitis_HLS/2024.1/settings64.sh && \
	source /home/levlabcukomen/tools/Vivado/2024.1/settings64.sh && \
	vitis_hls -f ../../$(SCRIPT_DIR)/hls_csim.tcl
	@echo ""
	@echo "=== C simulation completed ==="
	@echo "CSV outputs in: $(CSIM_BUILD_DIR)/"
	@ls -la $(CSIM_BUILD_DIR)/*.csv 2>/dev/null || echo "(no CSV files found — check sim log)"

# Run csim then plot results
csim_plot: csim
	@echo "Plotting results..."
	python3 $(SCRIPT_DIR)/plot_response.py $(CSIM_BUILD_DIR)

external_cores:
	@echo "Cloning external cores repository..."
	git submodule update --init --recursive

hls:
	@echo "Running HLS synthesis..."
	@mkdir -p $(HLS_DIR)
	@mkdir -p $(IP_REPO_DIR)
	cd $(HLS_DIR) && \
	source /home/levlabcukomen/tools/Vitis_HLS/2024.1/settings64.sh && \
	source /home/levlabcukomen/tools/Vivado/2024.1/settings64.sh && \
	vitis_hls -f ../../$(SCRIPT_DIR)/hls_build.tcl

vivado: start_monitor
	@echo "Running Vivado synthesis and implementation..."
	@mkdir -p $(VIVADO_DIR)
	cd $(VIVADO_DIR) && \
	source /home/levlabcukomen/tools/Vivado/2024.1/settings64.sh && \
	vivado -mode batch -source ../../$(SCRIPT_DIR)/vivado_build.tcl

clean:
	@echo "Cleaning build directories..."
	rm -rf $(BUILD_DIR)/
	rm -f *.log *.jou *.str
	@echo "Clean completed."


MONITOR_PID = $(BUILD_DIR)/monitor.pid

start_monitor:
	@echo "Checking bitstream monitor..."
	@mkdir -p $(BUILD_DIR)
	@if [ -f $(MONITOR_PID) ] && kill -0 $$(cat $(MONITOR_PID)) 2>/dev/null; then \
		echo "Monitor is already running (PID $$(cat $(MONITOR_PID)))."; \
	else \
		echo "Starting daemon..."; \
		/bin/bash $(MONITOR_SCRIPT) > $(MONITOR_LOG) 2>&1 & echo $$! > $(MONITOR_PID); \
		echo "Monitor started with PID $$(cat $(MONITOR_PID))."; \
	fi

stop_monitor:
	@echo "Stopping bitstream monitor..."
	@if [ -f $(MONITOR_PID) ]; then \
		kill $$(cat $(MONITOR_PID)) && rm $(MONITOR_PID) && echo "Stopped."; \
	else \
		echo "No PID file found. Checking for ghosts..."; \
		pkill -f "inotifywait.*build" || echo "Nothing found."; \
	fi