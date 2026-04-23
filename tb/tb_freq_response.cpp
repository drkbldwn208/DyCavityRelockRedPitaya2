#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "../src/dy_cavity_relocker_2.h"

/*
 * Frequency-response testbench for dy_cavity_relocker_2
 *
 * ADC sample rate:     125     MHz
 * Total decimation:    128x    (4x CIC + 32x CIC)
 * Filter sample rate:  976.5625 kHz
 * Nyquist:             488.28  kHz
 *
 * Each call to dy_cavity_relocker_2 processes 128 ADC samples.
 */

static const int    FRAME_SIZE = 128;
static const double FS_ADC     = 125.0e6;
static const double FS_FILT    = FS_ADC / FRAME_SIZE;  // 976562.5 Hz
static const int    INPUT_AMPL = 50;

static const int    N_SETTLE   = 128 * 8000;    // 8000 filter samples (~8 ms)
static const int    N_MEASURE  = 128 * 16000;   // 16000 filter samples
static const int    N_TOTAL    = N_SETTLE + N_MEASURE;

static const double F_START    = 5.0e3;          // skip below 5 kHz (too slow to settle)
static const double F_STOP     = 400.0e3;
static const int    N_FREQS    = 40;

// ---- helpers ----
static axis_t make_adc_word(short ch1, short ch2)
{
    axis_t w;
    w.data = (((uint32_t)(unsigned short)ch2) << 16) |
              ((uint32_t)(unsigned short)ch1);
    w.keep = 0xF;
    w.strb = 0xF;
    w.last = 0;
    return w;
}

static short extract_dac_ch1(const axis_t &w) { return (short)(w.data & 0xFFFF); }
static short extract_dac_ch2(const axis_t &w) { return (short)((w.data >> 16) & 0xFFFF); }

struct MeasResult { double gain_db, phase_deg, peak_out, rms_out; };

static MeasResult measure_freq(double freq, bool verbose)
{
    double omega_adc  = 2.0 * M_PI * freq / FS_ADC;
    double omega_filt = 2.0 * M_PI * freq / FS_FILT;
    int n_calls = N_TOTAL / FRAME_SIZE;
    int n_filt_total  = n_calls;
    int n_filt_settle = N_SETTLE / FRAME_SIZE;

    hls::stream<axis_t> adc_in("adc_in");
    hls::stream<axis_t> dac_out("dac_out");

    short  *out_filt = new short[n_filt_total];
    double *in_filt  = new double[n_filt_total];

    for (int call = 0; call < n_calls; call++) {
        double acc_ref = 0;
        for (int k = 0; k < FRAME_SIZE; k++) {
            int n = call * FRAME_SIZE + k;
            double val = INPUT_AMPL * cos(omega_adc * n);
            short ch1_in = (short)round(val);
            acc_ref += ch1_in;
            adc_in.write(make_adc_word(ch1_in, 0));
        }
        in_filt[call] = acc_ref / FRAME_SIZE;

        dy_cavity_relocker_2(adc_in, dac_out, false, 0, 0);

        for (int k = 0; k < FRAME_SIZE; k++) {
            axis_t out = dac_out.read();
            if (k == 0)
                out_filt[call] = extract_dac_ch1(out);
        }
    }

    double in_sum_sq = 0, out_sum_sq = 0;
    double sum_sin = 0, sum_cos = 0;
    int n_meas = 0;

    for (int k = n_filt_settle; k < n_filt_total; k++) {
        double y = (double)out_filt[k];
        double x = in_filt[k];
        in_sum_sq  += x * x;
        out_sum_sq += y * y;
        sum_sin += y * sin(omega_filt * k);
        sum_cos += y * cos(omega_filt * k);
        n_meas++;
    }

    double rms_in  = sqrt(in_sum_sq / n_meas);
    double rms_out = sqrt(out_sum_sq / n_meas);
    double gain_linear = (rms_in > 1e-12) ? rms_out / rms_in : 0;
    double gain_db = (gain_linear > 1e-12) ? 20.0 * log10(gain_linear) : -200.0;
    double phase_deg = atan2(sum_cos, sum_sin) * 180.0 / M_PI;

    if (verbose) {
        printf("  [%.1f Hz] rms_in=%.1f rms_out=%.1f gain=%.4f (%.2f dB) phase=%.1f\n",
               freq, rms_in, rms_out, gain_linear, gain_db, phase_deg);
        printf("    First 8 decimated outputs: ");
        for (int i = n_filt_settle; i < n_filt_settle + 8 && i < n_filt_total; i++)
            printf("%6d ", out_filt[i]);
        printf("\n    First 8 decimated inputs:  ");
        for (int i = n_filt_settle; i < n_filt_settle + 8 && i < n_filt_total; i++)
            printf("%6.0f ", in_filt[i]);
        printf("\n");
    }

    delete[] out_filt;
    delete[] in_filt;
    return {gain_db, phase_deg, (double)rms_out, rms_out};
}

const int PLANT_N_SEC = 5;
const double PLANT_SOS[][5] = {
    {5.05347664e+00, 5.05347664e+00, 0.00000000e+00, 1.18577075e-02, 0.00000000e+00},
    {1.00000000e+00, -1.99072016e+00, 9.90790930e-01, -1.87530636e+00, 8.79222958e-01},
    {1.00000000e+00, -1.99846379e+00, 9.99782834e-01, -1.99838334e+00, 9.99699627e-01},
    {1.00000000e+00, -2.00282966e+00, 1.00302868e+00, -1.99952581e+00, 9.99792197e-01},
    {1.00000000e+00, -1.99874668e+00, 9.99490582e-01, -1.99926010e+00, 9.99944465e-01},
};

// Simple floating-point Biquad cascade to simulate the physical table
class PlantSimulator {
private:
    double w1[20] = {0};
    double w2[20] = {0};

public:
    double process(double x) {
        double y = x;
        for (int i = 0; i < PLANT_N_SEC; i++) {
            double b0 = PLANT_SOS[i][0], b1 = PLANT_SOS[i][1], b2 = PLANT_SOS[i][2];
            double a1 = PLANT_SOS[i][3], a2 = PLANT_SOS[i][4];

            // Direct Form II (Floating point has massive dynamic range, so DF-II is perfectly stable here)
            double w0 = y - a1 * w1[i] - a2 * w2[i];
            y = b0 * w0 + b1 * w1[i] + b2 * w2[i];
            
            w2[i] = w1[i];
            w1[i] = w0;
        }
        return y;
    }
};

int main()
{
    printf("=== dy_cavity_relocker_2 Testbench ===\n");
    printf("  ADC rate:    %.2f MHz\n", FS_ADC/1e6);
    printf("  Filter rate: %.4f kHz  (128x decimation)\n", FS_FILT/1e3);
    printf("  Nyquist:     %.2f kHz\n", FS_FILT/2e3);
    printf("  Frame size:  %d ADC samples per call\n\n", FRAME_SIZE);

    // ---- Test 1: DC step ----
    // printf("--- Test 1: DC step (input=50) ---\n");
    // {
    //     hls::stream<axis_t> adc_in("dc_in");
    //     hls::stream<axis_t> dac_out("dc_out");
    //     FILE *f = fopen("step_response.csv", "w");
    //     fprintf(f, "frame,input,output_ch1,output_ch2\n");

    //     short last_ch1 = 0;
    //     int n_step_calls = 200;
    //     for (int call = 0; call < n_step_calls; call++) {
    //         for (int k = 0; k < FRAME_SIZE; k++)
    //             adc_in.write(make_adc_word(50, 0));
    //         dy_cavity_relocker_2(adc_in, dac_out, false, 0, 0);
    //         short y1 = 0, y2 = 0;
    //         for (int k = 0; k < FRAME_SIZE; k++) {
    //             axis_t out = dac_out.read();
    //             if (k == 0) {
    //                 y1 = extract_dac_ch1(out);
    //                 y2 = extract_dac_ch2(out);
    //             }
    //         }
    //         fprintf(f, "%d,%d,%d,%d\n", call, 50, y1, y2);
    //         last_ch1 = y1;
    //         if (call < 10 || call % 50 == 0)
    //             printf("  frame=%3d  ch1=%6d  ch2=%6d\n", call, y1, y2);
    //     }
    //     fclose(f);
    //     printf("  Final ch1=%d (expect ~2000 for unity-gain DC)\n\n", last_ch1);
    // }

    // // ---- Test 2: Impulse ----
    // printf("--- Test 2: Impulse (50 at n=0) ---\n");
    // {
    //     hls::stream<axis_t> adc_in("imp_in");
    //     hls::stream<axis_t> dac_out("imp_out");
    //     FILE *f = fopen("impulse_response.csv", "w");
    //     fprintf(f, "frame,input,output_ch1\n");

    //     int nonzero = 0; short maxabs = 0;
    //     int n_imp_calls = 300;
    //     for (int call = 0; call < n_imp_calls; call++) {
    //         for (int k = 0; k < FRAME_SIZE; k++) {
    //             int n = call * FRAME_SIZE + k;
    //             short imp = (n == 0) ? 50 : 0;
    //             adc_in.write(make_adc_word(imp, 0));
    //         }
    //         dy_cavity_relocker_2(adc_in, dac_out, false, 0, 0);
    //         short frame_out = 0;
    //         for (int k = 0; k < FRAME_SIZE; k++) {
    //             axis_t out = dac_out.read();
    //             if (k == 0) frame_out = extract_dac_ch1(out);
    //         }
    //         fprintf(f, "%d,%d,%d\n", call, (call == 0) ? 50 : 0, frame_out);
    //         if (frame_out != 0) nonzero++;
    //         if (abs(frame_out) > maxabs) maxabs = abs(frame_out);
    //         if (call < 20)
    //             printf("  frame=%3d  in=%5d  out=%6d\n", call,
    //                    (call == 0) ? 50 : 0, frame_out);
    //     }
    //     fclose(f);
    //     printf("  Nonzero frames: %d/%d  Max|out|=%d\n", nonzero, n_imp_calls, maxabs);
    //     if (nonzero == 0) printf("  ERROR: no output from impulse!\n");
    //     printf("\n");
    // }

    // // ---- Test 3: Spot checks ----
    // printf("--- Test 3: Spot frequency checks ---\n");
    // double spots[] = {1e3, 5e3, 15e3, 50e3, 200e3};
    // for (int i = 0; i < 5; i++)
    //     measure_freq(spots[i], true);
    // printf("\n");
    #include "../src/hinf_filter.hpp"

// ... inside main() ...

    printf("--- Test 5a: FLOATING-POINT closed loop (verify control design) ---\n");
    {
        PlantSimulator plant;
        // Floating-point biquad cascade matching the controller's SOS
        PlantSimulator controller_float;  // We need a float version of the controller
        
        // Instead, use the EXACT same controller but in float:
        // Direct Form II, float64, same coefficients
        double cw1[20] = {0}, cw2[20] = {0};
        
        FILE *f = fopen("closed_loop_float.csv", "w");
        fprintf(f, "sample,error,effort\n");
        
        double disturbance = 50.0;
        double y = 0;  // plant output
        
        for (int k = 0; k < 5000; k++) {
            double e = disturbance - y;
            
            // Float controller using same SOS coefficients
            double u = e;
            for (int i = 0; i < HINF_N_SECTIONS; i++) {
                double b0 = (double)HINF_SOS[i].b0 / HINF_COEF_SCALE;
                double b1 = (double)HINF_SOS[i].b1 / HINF_COEF_SCALE;
                double b2 = (double)HINF_SOS[i].b2 / HINF_COEF_SCALE;
                double a1 = (double)HINF_SOS[i].a1 / HINF_COEF_SCALE;
                double a2 = (double)HINF_SOS[i].a2 / HINF_COEF_SCALE;
                
                double w0 = u - a1 * cw1[i] - a2 * cw2[i];
                u = b0 * w0 + b1 * cw1[i] + b2 * cw2[i];
                cw2[i] = cw1[i];
                cw1[i] = w0;
            }
            
            y = plant.process(u);
            fprintf(f, "%d,%.6f,%.6f\n", k, e, u);
            
            if (k < 20 || k % 500 == 0)
                printf("  k=%4d  error=%.4f  effort=%.4f  y=%.4f\n", k, e, u, y);
        }
        fclose(f);
    }\

      printf("--- Test 5b: FIXED-POINT closed loop ---\n");
    {
        HinfFilter lock_controller;
        PlantSimulator plant;
        
        FILE *f = fopen("closed_loop_fixed.csv", "w");
        fprintf(f, "sample,error,effort\n");
        
        double disturbance = 50.0;
        double y = 0;
        
        for (int k = 0; k < 5000; k++) {
            double e = disturbance - y;
            
            // Clip to Q1.15 range before loading
            short e_short = (short)round(fmax(-32768, fmin(32767, e)));
            sig_t u_in;
            u_in.range() = e_short;
            
            sig_t y_ctrl = lock_controller.process(u_in);
            short effort_short = (short)y_ctrl.range();
            
            // Controller output in Q1.15: effort_short/32768 is the actual value
            // Plant expects the same scaling as the synthesis
            double effort_float = (double)effort_short / 32768.0;
            y = plant.process(effort_float) * 32768.0;
            
            fprintf(f, "%d,%.6f,%d\n", k, e, effort_short);
            
            if (k < 20 || k % 500 == 0)
                printf("  k=%4d  error=%.4f  effort=%d  y=%.4f\n", k, e, effort_short, y);
        }
        fclose(f);
    }

    // ---- Test 4: Full sweep ----
    printf("--- Test 4: Frequency sweep (%d pts, %.0f Hz - %.0f kHz) ---\n",
           N_FREQS, F_START, F_STOP/1e3);
    FILE *csv = fopen("freq_response.csv", "w");
    fprintf(csv, "freq_hz,freq_mhz,gain_db,phase_deg,peak_out,rms_out\n");
    printf("%-14s %-12s %-12s %-12s\n", "Freq (kHz)", "Gain (dB)", "RMS out", "Phase");
    printf("------------------------------------------------------\n");

    for (int fi = 0; fi < N_FREQS; fi++) {
        double freq = F_START * pow(F_STOP / F_START, (double)fi / (N_FREQS - 1));
        MeasResult r = measure_freq(freq, false);
        printf("%-14.4f %-12.2f %-12.1f %-12.1f\n",
               freq/1e3, r.gain_db, r.rms_out, r.phase_deg);
        fprintf(csv, "%.6f,%.6f,%.2f,%.2f,%.1f,%.1f\n",
                freq, freq/1e6, r.gain_db, r.phase_deg, r.peak_out, r.rms_out);
    }
    fclose(csv);

    printf("\nCSVs: freq_response.csv, step_response.csv, impulse_response.csv\n");
    return 0;
}