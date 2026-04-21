#include "dy_cavity_relocker_2.h"
#include "hinf_filter.hpp"

void dy_cavity_relocker_2(hls::stream<axis_t> &adc_in,
                          hls::stream<axis_t> &dac_out,
                          bool gpio_in,
                          int servo_offset,
                          int servo_arm)
{
  #pragma HLS PIPELINE II = 4
  #pragma HLS INTERFACE axis port = adc_in
  #pragma HLS INTERFACE axis port = dac_out
  #pragma HLS INTERFACE s_axilite port = servo_offset
  #pragma HLS INTERFACE s_axilite port = servo_arm
  #pragma HLS INTERFACE ap_none port = gpio_in
  #pragma HLS INTERFACE s_axilite port = return

  static State current_state = IDLE;
  bool state_toggle = false;

  static HinfFilter controller;

  static ap_uint<2> decim_cnt = 0;
  static int32_t adc1_acc = 0;
  static short hinf_out_hold = 0;

  fsm_receiver(gpio_in, state_toggle);

  axis_t val_adc = adc_in.read();
  short val_adc1 = (short)(val_adc.data & 0xFFFF);         // Channel A
  short val_adc2 = (short)((val_adc.data >> 16) & 0xFFFF); // Right shift reads the upper 16 bits, so this is channel B

  /*
  Channel 1 handling: x4 decimation and H-infinity filtering
  */
  adc1_acc += val_adc1;

  if (decim_cnt == 3) {
    short averaged_adc1 = (short)(adc1_acc >> 2); // Divide by 4 to get the average
    adc1_acc = 0;

    sig_t u_in;
    u_in.range() = averaged_adc1;

    sig_t y_out = controller.process(u_in);

    hinf_out_hold = (short)y_out.range(); // Hold the H-infinity output until the next decimated sample is ready

  }
  decim_cnt++;

  axis_t val_dac;
  short dac_1_out = hinf_out_hold; // Output the held H-infinity output to DAC channel 1
  
  /*
  Channel 2 handling: Cavity relocking servo logic based on the state machine
  */

  short dac_2_out = 0;

  static short held_voltage = 0;
  static short error_signal = 0;

  switch (current_state)
  {
  case IDLE:
    dac_2_out = servo_offset; // Output the servo offset to DAC channel 2 in IDLE state
    
    if (servo_arm)
    {
      if (state_toggle)
      {
        current_state = SERVO;
        held_voltage = val_adc2; // Capture the current ADC value to hold in SERVO state
      }
    }
    break;


  case SERVO:
    error_signal = val_adc2 - held_voltage;
    dac_2_out = servo_offset + (error_signal >> GAIN_RIGHT_SHIFT); // Output the servo offset plus the error signal to DAC channel 2 in SERVO state
    
    if (state_toggle)
    {
      current_state = IDLE;
    }
    break;
  }

  val_dac.data = (((uint32_t)dac_2_out & 0xFFFF) << 16) | ((uint32_t)dac_1_out & 0xFFFF);
  val_dac.keep = 0xF;
  val_dac.strb = 0xF;
  val_dac.last = 0;

  dac_out.write(val_dac);
}
