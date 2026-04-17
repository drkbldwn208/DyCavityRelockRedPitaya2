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

.PHONY: all csim hls vivado external_cores clean

all: external_cores hls vivado

csim:
	@echo "Running C simulation..."
	@mkdir -p $(SIM_DIR)
	g++ $(SOURCE_FILES) $(TESTBENCH) -o $(SIM_DIR)/sim_extecutable
	cd $(SIM_DIR) && ./sim_extecutable
	@echo "C simulation completed."

external_cores:
	@echo "Cloning external cores repository..."
	git submodule update --init --recursive

hls:
	@echo "Running HLS synthesis..."
	@mkdir -p $(HLS_DIR)
	mkdir -p $(IP_REPO_DIR)
	cd $(HLS_DIR) && vitis_hls -f ../../$(SCRIPT_DIR)/hls_build.tcl

vivado:
	@echo "Running Vivado synthesis and implementation..."
	@mkdir -p $(VIVADO_DIR)
	cd $(VIVADO_DIR) && vivado -mode batch -source ../../$(SCRIPT_DIR)/vivado_build.tcl

clean:
	@echo "Cleaning build directories..."
	rm -rf $(BUILD_DIR)/
	rm -f *.log *.jou *.str
	@echo "Clean completed."
