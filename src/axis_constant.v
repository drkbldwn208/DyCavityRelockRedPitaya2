module axis_constant(
    input wire clk,
    output wire [31:0] m_axis_tdata,
    output wire m_axis_tvalid
    );
    
    
    assign m_axis_tdata = 32'h1770_1770;
    assign m_axis_tvalid = 1'b1;
    
    
endmodule
