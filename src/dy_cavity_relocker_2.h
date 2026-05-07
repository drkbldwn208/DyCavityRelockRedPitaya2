#ifndef DY_CAVITY_RELOCKER_2_H
#define DY_CAVITY_RELOCKER_2_H

#include <stdint.h>
#include "ap_fixed.h"
#include "ap_axi_sdata.h"
#include "hls_stream.h"
#include "ap_int.h"

#include "fsm_receiver.h"

typedef ap_axiu<32, 0, 0, 0> axis_t;
typedef ap_axiu<16, 0, 0, 0> axis_t_16;

enum State
{
  IDLE,
  SERVO
};

#ifndef GAIN_RIGHT_SHIFT
  #define GAIN_RIGHT_SHIFT 4
#endif

void dy_cavity_relocker_2(hls::stream<axis_t> &adc_in,
                          hls::stream<axis_t> &dac_out,
                          volatile bool *gpio_in,
                          volatile int *servo_offset,
                          volatile int *servo_arm);

#endif