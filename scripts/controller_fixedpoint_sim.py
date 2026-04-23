import numpy as np

# Your exact coefficients from hinf_coeffs.h
HINF_SOS = [
    [ 178688417,    2118835,          0,   48158576,          0],
    [ 536870912, -816207689,  294793619, -860243371,  335489729],
    [ 536870912,-1006794266,  472026296,-1070291980,  533427175],
    [ 536870912,-1053405531,  516534750,-1064714812,  527851246],
    [ 536870912,-1072900259,  536736258,-1072882850,  536726797],
    [ 536870912,-1073257134,  536752120,-1073298003,  536798746],
    [ 536870912,-1073551492,  536825045,-1073533866,  536808406]
]
SCALE = 536870912.0 # 2^29

def run_fp_sim(x_in):
    x_hist = np.zeros((7, 2))
    y_hist = np.zeros((7, 2))
    y_out = np.zeros_like(x_in, dtype=float)
    
    for n in range(len(x_in)):
        pipe = float(x_in[n])
        for sec in range(7):
            b0, b1, b2, a1, a2 = HINF_SOS[sec]
            
            # Simulate the 64-bit Direct Form I accumulator
            acc = (b0 * pipe) + (b1 * x_hist[sec, 0]) + (b2 * x_hist[sec, 1]) \
                  - (a1 * y_hist[sec, 0]) - (a2 * y_hist[sec, 1])
            
            # Truncate back to pipe width
            pipe_out = int(acc / SCALE)
            
            # Shift history
            x_hist[sec, 1] = x_hist[sec, 0]
            x_hist[sec, 0] = pipe
            y_hist[sec, 1] = y_hist[sec, 0]
            y_hist[sec, 0] = pipe_out
            
            pipe = pipe_out
            
        y_out[n] = pipe
    return y_out

# Run the 200-sample Step Response
step_in = np.ones(200) * 50
step_out = run_fp_sim(step_in)

print("Python Fixed-Point Step Response (Last 10 samples):")
for i in range(190, 200):
    print(f" n={i}  out_ch1= {int(step_out[i])}")