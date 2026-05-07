#include "dy_cavity_relocker_2.h"
#include "hinf_filter.hpp"
#include <ap_int.h>
#include <cstdint>

struct fast_t  { short ch2; };           // every cycle
struct slow_t  { short d128; };          // every 128 cycles  
struct ctrl_t  { short dac1; };          // every 128 cycles

// Process 1: ADC + decimation. Same rate as input on ch2_s, gated on slow_s.
static void front_end(hls::stream<axis_t> &adc_in,
                      hls::stream<fast_t> &ch2_s,
                      hls::stream<slow_t> &slow_s) {
    static int32_t acc4 = 0; static ap_uint<2> cnt4 = 0;
    static int64_t acc32 = 0; static ap_uint<5> cnt32 = 0;
    while (true) {
        #pragma HLS PIPELINE II=1
        axis_t v = adc_in.read();
        short ch1 = (short)(v.data & 0xFFFF);
        ch2_s.write({ (short)((v.data >> 16) & 0xFFFF) });

        acc4 += ch1;
        if (cnt4 == 3) {
            acc32 += (acc4 >> 2); acc4 = 0; cnt4 = 0;
            if (cnt32 == 31) {
                slow_s.write({ (short)(acc32 >> 5) });
                acc32 = 0; cnt32 = 0;
            } else cnt32++;
        } else cnt4++;
    }
}

// Process 2: controller. Blocking read → only runs when data arrives (1/128).
// II=8 gives the synthesizer 8 cycles to close the recurrence.
static void controller_stage(hls::stream<slow_t> &slow_s,
                             hls::stream<ctrl_t> &ctrl_s) {
    static HinfFilter controller;
    while (true) {
        #pragma HLS PIPELINE II=8
        slow_t in = slow_s.read();
        sig_t u; u.range() = in.d128;
        ctrl_s.write({ (short)controller.process(u).range() });
    }
}

// Process 3: hold-last + servo + pack. II=1, this is where read_nb is correct.
static void back_end(hls::stream<fast_t> &ch2_s,
                     hls::stream<ctrl_t> &ctrl_s,
                     hls::stream<axis_t> &dac_out,
                     volatile bool *gpio_in,
                     volatile int *servo_offset,
                     volatile int *servo_arm) {
    static short dac1_held = 0;
    static State current_state = IDLE;
    static short held_voltage = 0;
    while (true) {
        #pragma HLS PIPELINE II=1
        fast_t f = ch2_s.read();        // blocking — same rate as us
        ctrl_t c;
        if (ctrl_s.read_nb(c))          // nonblocking — hold-last
            dac1_held = c.dac1;

        bool toggle = false;
        fsm_receiver(*gpio_in, toggle);
        short dac2_v;
        switch (current_state) {
          case IDLE:
            dac2_v = (short)(*servo_offset);
            if (*servo_arm && toggle) { current_state = SERVO; held_voltage = f.ch2; }
            break;
          case SERVO: default: {
            short err = f.ch2 - held_voltage;
            dac2_v = (short)(*servo_offset + (err >> GAIN_RIGHT_SHIFT));
            if (toggle) current_state = IDLE;
            break;
          }
        }

        axis_t o;
        o.data = (((uint32_t)dac2_v & 0xFFFF) << 16) | ((uint32_t)dac1_held & 0xFFFF);
        o.keep = 0xF; o.strb = 0xF; o.last = 0;
        dac_out.write(o);
    }
}

void dy_cavity_relocker_2(hls::stream<axis_t> &adc_in,
                          hls::stream<axis_t> &dac_out,
                          volatile bool *gpio_in,
                          volatile int *servo_offset,
                          volatile int *servo_arm) {
    #pragma HLS INTERFACE axis port=adc_in
    #pragma HLS INTERFACE axis port=dac_out
    #pragma HLS INTERFACE s_axilite port=servo_offset
    #pragma HLS INTERFACE s_axilite port=servo_arm
    #pragma HLS INTERFACE ap_none port=gpio_in
    #pragma HLS INTERFACE ap_ctrl_none port=return
    #pragma HLS DATAFLOW

    hls::stream<fast_t> ch2_s;   
    #pragma HLS STREAM variable=ch2_s depth=4
    hls::stream<slow_t> slow_s;  
    #pragma HLS STREAM variable=slow_s depth=4
    hls::stream<ctrl_t> ctrl_s;  
    #pragma HLS STREAM variable=ctrl_s depth=4

    front_end(adc_in, ch2_s, slow_s);
    controller_stage(slow_s, ctrl_s);
    back_end(ch2_s, ctrl_s, dac_out, gpio_in, servo_offset, servo_arm);
}