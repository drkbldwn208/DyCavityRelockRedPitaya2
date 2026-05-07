
module pdm_wrapper #(
    parameter DWC=8,
    parameter CHN=4
    )(
    input wire clk,
    input wire rstn,
    input wire [DWC*CHN-1:0] cfg,
    input wire ena,
    input wire [DWC-1:0] rng,
    output wire [CHN-1:0] pdm
    );
    
red_pitaya_pdm #(
    .DWC(DWC),
    .CHN(CHN)
    ) pdm_inst (
    .clk(clk),
    .rstn(rstn),
    .cfg_flat(cfg),
    .ena(ena),
    .rng(rng),
    .pdm(pdm)   
);
endmodule