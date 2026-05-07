#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "../src/dy_cavity_relocker_2.h"
#include "../src/hinf_filter.hpp"

/*
 * Testbench for HinfFilter (the controller inside dy_cavity_relocker_2).
 *
 * NOTE: This testbench does NOT call the top-level kernel.
 * The kernel uses while(true) + DATAFLOW + hls::stream and cannot run
 * in csim — the streaming model isn't simulated. Use cosim_design or
 * hardware to validate the full pipeline.
 *
 * What this DOES test (all in the post-decimation domain at FS_FILT):
 *   1  DC step                 → step_response.csv
 *   2  Impulse                 → impulse_response.csv
 *   3  Float closed loop       → closed_loop_float.csv
 *   4  Fixed-point closed loop → closed_loop_fixed.csv
 *   5  Frequency sweep         → freq_response.csv
 */

static const int    FRAME_SIZE = 128;
static const double FS_ADC     = 125.0e6;
static const double FS_FILT    = FS_ADC / FRAME_SIZE;   // 976562.5 Hz
static const int    INPUT_AMPL = 50;                    // raw Q1.15 ticks

// ----- frequency sweep settings (in filter-rate samples) -----
static const int    N_SETTLE   = 8000;
static const int    N_MEASURE  = 16000;
static const double F_START    = 5.0e3;
static const double F_STOP     = 400.0e3;
static const int    N_FREQS    = 40;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static inline sig_t short_to_sig(short s) {
    sig_t u; u.range() = s; return u;
}
static inline short sig_to_short(sig_t s) {
    return (short)s.range();
}

struct MeasResult { double gain_db, phase_deg, peak_out, rms_out; };

// Plant SOS — paste from ctrl.py "Plant SOS for tb_freq_response.cpp"
const int PLANT_N_SEC = 1;
const double PLANT_SOS[][5] = {
    {7.01024743e-02, 1.40204949e-01, 7.01024743e-02, -7.11018449e-01, -8.57165403e-03},
};

class PlantSimulator {
    double w1[10] = {0};
    double w2[10] = {0};
public:
    double process(double x) {
        double y = x;
        for (int i = 0; i < PLANT_N_SEC; i++) {
            double b0 = PLANT_SOS[i][0], b1 = PLANT_SOS[i][1], b2 = PLANT_SOS[i][2];
            double a1 = PLANT_SOS[i][3], a2 = PLANT_SOS[i][4];
            double w0 = y - a1 * w1[i] - a2 * w2[i];
            y = b0 * w0 + b1 * w1[i] + b2 * w2[i];
            w2[i] = w1[i];
            w1[i] = w0;
        }
        return y;
    }
};

// ---------------------------------------------------------------------------
// Frequency-sweep core: drives HinfFilter directly at FS_FILT
// ---------------------------------------------------------------------------
static MeasResult measure_freq(double freq, bool verbose)
{
    HinfFilter ctrl;
    double omega = 2.0 * M_PI * freq / FS_FILT;

    double in_ss = 0, out_ss = 0, sum_sin = 0, sum_cos = 0;
    int n_meas = 0;

    for (int n = 0; n < N_SETTLE + N_MEASURE; n++) {
        double x = (double)INPUT_AMPL * cos(omega * n);
        short  x_short = (short)round(x);
        sig_t  y_ctrl = ctrl.process(short_to_sig(x_short));
        short  y = sig_to_short(y_ctrl);

        if (n >= N_SETTLE) {
            in_ss   += x * x;
            out_ss  += (double)y * y;
            sum_sin += y * sin(omega * n);
            sum_cos += y * cos(omega * n);
            n_meas++;
        }
    }

    double rms_in      = sqrt(in_ss  / n_meas);
    double rms_out     = sqrt(out_ss / n_meas);
    double gain_lin    = (rms_in > 1e-12) ? rms_out / rms_in : 0;
    double gain_db     = (gain_lin > 1e-12) ? 20.0 * log10(gain_lin) : -200.0;
    double phase_deg   = atan2(sum_cos, sum_sin) * 180.0 / M_PI;

    if (verbose) {
        printf("  [%.1f Hz] rms_in=%.3f rms_out=%.3f gain=%.2f dB phase=%.1f\n",
               freq, rms_in, rms_out, gain_db, phase_deg);
    }
    return {gain_db, phase_deg, rms_out, rms_out};
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main()
{
    printf("=== HinfFilter testbench ===\n");
    printf("  Filter rate: %.4f kHz  (post 128x decim of %.2f MHz ADC)\n",
           FS_FILT / 1e3, FS_ADC / 1e6);
    printf("  Nyquist:     %.2f kHz\n\n", FS_FILT / 2e3);

    // --- Test 1: DC step ---
    printf("--- Test 1: DC step (input = %d) ---\n", INPUT_AMPL);
    {
        HinfFilter ctrl;
        FILE *f = fopen("step_response.csv", "w");
        fprintf(f, "sample,input,output\n");

        const int n_samples = 2000;
        short last_y = 0;
        for (int n = 0; n < n_samples; n++) {
            sig_t y_ctrl = ctrl.process(short_to_sig((short)INPUT_AMPL));
            short y = sig_to_short(y_ctrl);
            fprintf(f, "%d,%d,%d\n", n, INPUT_AMPL, y);
            last_y = y;
        }
        fclose(f);
        printf("  Final output = %d\n\n", last_y);
    }

    // --- Test 2: Impulse ---
    printf("--- Test 2: Impulse (amp=%d at n=0) ---\n", INPUT_AMPL);
    {
        HinfFilter ctrl;
        FILE *f = fopen("impulse_response.csv", "w");
        fprintf(f, "sample,input,output\n");

        const int n_samples = 3000;
        int    nonzero = 0;
        short  maxabs  = 0;
        for (int n = 0; n < n_samples; n++) {
            short x = (n == 0) ? (short)INPUT_AMPL : 0;
            sig_t y_ctrl = ctrl.process(short_to_sig(x));
            short y = sig_to_short(y_ctrl);
            fprintf(f, "%d,%d,%d\n", n, x, y);
            if (y != 0) nonzero++;
            if (abs(y) > maxabs) maxabs = abs(y);
        }
        fclose(f);
        printf("  Nonzero samples: %d/%d   max|out| = %d\n", nonzero, n_samples, maxabs);
        if (nonzero == 0) printf("  ERROR: no output from impulse!\n");
        printf("\n");
    }

    // --- Test 3: Floating-point closed loop (reference) ---
    printf("--- Test 3: Floating-point closed loop ---\n");
    {
        PlantSimulator plant;
        double cw1[20] = {0}, cw2[20] = {0};
        FILE *f = fopen("closed_loop_float.csv", "w");
        fprintf(f, "sample,error,effort\n");

        const int n_samples = 5000;
        const double disturbance = 50.0;
        double y = 0;
        for (int k = 0; k < n_samples; k++) {
            double e = disturbance - y;
            double u = e;
            for (int i = 0; i < HINF_N_SECTIONS; i++) {
                double b0 = (double)HINF_SOS[i].b0 / HINF_COEF_SCALE;
                double b1 = (double)HINF_SOS[i].b1 / HINF_COEF_SCALE;
                double b2 = (double)HINF_SOS[i].b2 / HINF_COEF_SCALE;
                double a1 = (double)HINF_SOS[i].a1 / HINF_COEF_SCALE;
                double a2 = (double)HINF_SOS[i].a2 / HINF_COEF_SCALE;
                double w0 = u - a1 * cw1[i] - a2 * cw2[i];
                u  = b0 * w0 + b1 * cw1[i] + b2 * cw2[i];
                cw2[i] = cw1[i];
                cw1[i] = w0;
            }
            y = plant.process(u);
            fprintf(f, "%d,%.6f,%.6f\n", k, e, u);
        }
        fclose(f);
        printf("  Final: error=%.3f  y=%.3f\n\n", disturbance - y, y);
    }

    // --- Test 4: Fixed-point closed loop (what the FPGA does) ---
    printf("--- Test 4: Fixed-point closed loop ---\n");
    {
        HinfFilter ctrl;
        PlantSimulator plant;
        FILE *f = fopen("closed_loop_fixed.csv", "w");
        fprintf(f, "sample,error,effort\n");

        const int n_samples = 5000;
        const double disturbance = 50.0;
        double y = 0;
        for (int k = 0; k < n_samples; k++) {
            double e = disturbance - y;
            short  e_short = (short)round(fmax(-32768.0, fmin(32767.0, e)));
            sig_t  y_ctrl  = ctrl.process(short_to_sig(e_short));
            short  effort  = sig_to_short(y_ctrl);
            y = plant.process((double)effort);
            fprintf(f, "%d,%.6f,%d\n", k, e, effort);
        }
        fclose(f);
        printf("  Final: error=%.3f  y=%.3f\n\n", disturbance - y, y);
    }

    // --- Test 5: Frequency sweep ---
    printf("--- Test 5: Frequency sweep (%d pts, %.0f Hz - %.0f kHz) ---\n",
           N_FREQS, F_START, F_STOP / 1e3);
    {
        FILE *csv = fopen("freq_response.csv", "w");
        fprintf(csv, "freq_hz,freq_mhz,gain_db,phase_deg,peak_out,rms_out\n");
        printf("%-14s %-12s %-12s %-12s\n", "Freq (kHz)", "Gain (dB)", "RMS out", "Phase");
        printf("------------------------------------------------------\n");
        for (int fi = 0; fi < N_FREQS; fi++) {
            double freq = F_START * pow(F_STOP / F_START, (double)fi / (N_FREQS - 1));
            MeasResult r = measure_freq(freq, false);
            printf("%-14.4f %-12.2f %-12.3f %-12.1f\n",
                   freq / 1e3, r.gain_db, r.rms_out, r.phase_deg);
            fprintf(csv, "%.6f,%.6f,%.2f,%.2f,%.3f,%.3f\n",
                    freq, freq / 1e6, r.gain_db, r.phase_deg, r.peak_out, r.rms_out);
        }
        fclose(csv);
    }

    printf("\nCSVs: step_response.csv, impulse_response.csv, "
           "closed_loop_float.csv, closed_loop_fixed.csv, freq_response.csv\n");
    return 0;
}