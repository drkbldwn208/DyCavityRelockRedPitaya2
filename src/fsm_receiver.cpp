#include "dy_cavity_relocker_2.h"

void fsm_receiver(bool gpio_in, bool &fsm_trigger_out)
{
#pragma HLS INLINE

  static bool sync1 = false;
  static bool sync2 = false;
  static bool prev_sync2 = false;

  sync2 = sync1;
  sync1 = gpio_in;

  if (sync2 != prev_sync2)
  {
    fsm_trigger_out = true;
  }
  else
  {
    fsm_trigger_out = false;
  }

  prev_sync2 = sync2;
}