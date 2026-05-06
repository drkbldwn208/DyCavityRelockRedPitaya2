#include "dy_cavity_relocker_2.h"
#include "hinf_filter.hpp"

static void read_adc(hls::stream<axis_t> &adc_in,
                     hls::stream<short> &ch1_out,
                     hls::stream<short> &ch2_out)
{
  while(true) {
    #pragma HLS PIPELINE II=1
    axis_t val_adc = adc_in.read();
    ch1_out.write((short)(val_adc.data & 0xFFFF));         // Channel A
    ch2_out.write((short)((val_adc.data >> 16) & 0xFFFF)); // Channel B
    }
  }

static void decimate_by_4(hls::stream<short> &in, hls::stream<short> &out) {
  int32_t acc = 0;
  int count = 0;
  while(true) {
    #pragma HLS PIPELINE II=1
    acc += in.read();
    count++;

    if (count==4) {
      out.write((short)(acc >> 2)); // Average of 4 samples
      acc = 0;
      count = 0;
    }
  }
}

static void decimate_by_32(hls::stream<short> &in,
                            hls::stream<short> &out) {
  long acc = 0;
  short count = 0;
  while(true) {
    #pragma HLS PIPELINE II=1
    acc += in.read();
    count++;

    if (count == 32) {
      out.write((short)(acc >> 5)); // Average of 32 samples
      acc = 0;
      count = 0;
    }
  }
}

static void hinf_path(hls::stream<short> &ch1_in, 
                      hls::stream<short> &dac1_out) {
  while(true) {
    #pragma HLS PIPELINE II=1
    static HinfFilter controller;
    short averaged = ch1_in.read();
    sig_t u_in;
    
    u_in.range() = averaged;
    dac1_out.write((short)controller.process(u_in).range());
  }
}

static void servo_path(hls::stream<short> &ch2_in, 
                       hls::stream<short> &dac2_out,
                       bool gpio_in,
                       int servo_offset,
                       int servo_arm) {
  State current_state = IDLE;
  short held_voltage = 0;
  short error_signal = 0;
  
  while(true) {
    #pragma HLS PIPELINE II=1
    bool state_toggle = false;
    fsm_receiver(gpio_in, state_toggle);
    short val_adc2 = ch2_in.read();
    short out=0;

    switch (current_state) {
      case IDLE:
        out=servo_offset; // Output the servo offset to DAC channel 2 in IDLE state
        if (servo_arm && state_toggle) {
          current_state = SERVO;
          held_voltage = val_adc2; // Capture the current ADC value to hold in SERVO state
        }
        break;

      case SERVO:
        error_signal = val_adc2 - held_voltage;
        out = (servo_offset + (error_signal >> GAIN_RIGHT_SHIFT)); // Output the servo offset plus the error signal to DAC channel 2 in SERVO state
        
        if (state_toggle)
        {
          current_state = IDLE;
        }
        break;
    }
    dac2_out.write(out);
  }
}

static void hold_128(hls::stream<short> &in,
                   hls::stream<short> &out) {

  short current_val = 0;
  short new_val;
  while(true) {
    #pragma HLS PIPELINE II=1
    if (in.read_nb(new_val)) {
      current_val = new_val;
    }
    out.write(current_val);
  }
}


static void pack_dac(hls::stream<short> &dac1_in,
                     hls::stream<short> &dac2_in,
                     hls::stream<axis_t> &dac_out) {
  while(true) {
    #pragma HLS PIPELINE II=1
    axis_t val_dac;
    // short dac_1_out = dac1_in.read();
    // short dac_2_out = dac2_in.read();
    dac1_in.read(); // D/iscard the H-infinity path output for now
    dac2_in.read(); // Discard the servo path output for now
    short dac_1_out = 200;
    short dac_2_out = -200;
    val_dac.data = (((uint32_t)dac_2_out & 0xFFFF) << 16) | ((uint32_t)dac_1_out & 0xFFFF);
    val_dac.keep = 0xF;
    val_dac.strb = 0xF;
    val_dac.last = 0;
    
    dac_out.write(val_dac);
  }
}


void dy_cavity_relocker_2(hls::stream<axis_t> &adc_in,
                          hls::stream<axis_t> &dac_out,
                          bool gpio_in,
                          int servo_offset,
                          int servo_arm)
{
  #pragma HLS DATAFLOW
  #pragma HLS INTERFACE axis port = adc_in
  #pragma HLS INTERFACE axis port = dac_out
  #pragma HLS INTERFACE s_axilite port = servo_offset
  #pragma HLS INTERFACE s_axilite port = servo_arm
  #pragma HLS INTERFACE ap_none port = gpio_in
  #pragma HLS INTERFACE ap_ctrl_none port = return

  hls::stream<short> ch1, ch1_decim, ch1_decim_32, ch2, dac1_raw, dac1, dac2;
  #pragma HLS STREAM variable = ch2 depth = 256
  #pragma HLS STREAM variable = dac2 depth = 256
  
  read_adc(adc_in, ch1, ch2);
  decimate_by_4(ch1, ch1_decim);
  decimate_by_32(ch1_decim, ch1_decim_32);
  hinf_path(ch1_decim_32, dac1_raw);
  hold_128(dac1_raw, dac1);
  servo_path(ch2, dac2, gpio_in, servo_offset, servo_arm);
  pack_dac(dac1, dac2, dac_out);
}
