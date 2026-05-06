#include "dy_cavity_relocker_2.h"
#include "hinf_filter.hpp"
#include <ap_int.h>
#include <cstdint>

void dy_cavity_relocker_2(hls::stream<axis_t> &adc_in,
                          hls::stream<axis_t> &dac_out,
                          volatile bool *gpio_in,
                          volatile int *servo_offset,
                          volatile int *servo_arm)
{
    #pragma HLS INTERFACE axis port=adc_in
    #pragma HLS INTERFACE axis port=dac_out
    #pragma HLS INTERFACE s_axilite port=servo_offset
    #pragma HLS INTERFACE s_axilite port=servo_arm
    #pragma HLS INTERFACE ap_none port=gpio_in
    #pragma HLS INTERFACE ap_ctrl_none port=return

    // Decimation accumulators
    static int32_t acc4   = 0; static ap_uint<2> cnt4  = 0;
    static int64_t acc32  = 0; static ap_uint<5> cnt32 = 0;
    static HinfFilter controller;
    static short dac1_held = 0;

    // Servo FSM state
    static State current_state = IDLE;
    static short held_voltage  = 0;

    while (true) {
        #pragma HLS PIPELINE II=1

        axis_t v = adc_in.read();
        short ch1 = (short)(v.data & 0xFFFF);
        short ch2 = (short)((v.data >> 16) & 0xFFFF);

        // CIC-ish decimate by 4
        acc4 += ch1;
        bool d4_valid = (cnt4 == 3);
        short d4 = (short)(acc4 >> 2);
        if (d4_valid) { acc4 = 0; cnt4 = 0; } else { cnt4++; }

        // Decimate by 32 on top of the by-4 stream → 1/128 overall
        if (d4_valid) {
            acc32 += d4;
            bool d32_valid = (cnt32 == 31);
            if (d32_valid) {
                short d32 = (short)(acc32 >> 5);
                sig_t u; u.range() = d32;
                dac1_held = (short)controller.process(u).range();
                acc32 = 0; cnt32 = 0;
            } else { cnt32++; }
        }

        // Servo path (every cycle)
        bool toggle = false;
        fsm_receiver(*gpio_in, toggle);
        short dac2_v;
        switch (current_state) {
          case IDLE:
            dac2_v = (short)(*servo_offset);
            if (*servo_arm && toggle) { current_state = SERVO; held_voltage = ch2; }
            break;
          case SERVO:
          default: {
            short err = ch2 - held_voltage;
            dac2_v = (short)(*servo_offset + (err >> GAIN_RIGHT_SHIFT));
            if (toggle) current_state = IDLE;
            break;
          }
        }

        // Pack & emit every cycle → tvalid is essentially constant-1
        axis_t o;
        o.data = (((uint32_t)dac2_v   & 0xFFFF) << 16)
               | (((uint32_t)dac1_held) & 0xFFFF);
        o.keep = 0xF; o.strb = 0xF; o.last = 0;
        dac_out.write(o);
    }
}