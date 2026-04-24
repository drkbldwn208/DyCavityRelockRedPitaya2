#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "../src/dy_cavity_relocker_2.h"
#include "../src/hinf_filter.hpp"

/*
 * Testbench for dy_cavity_relocker_2.
 *
 *   ADC 125 MHz  →  128x decim  →  filter at 976.5625 kHz
 *
 *   1  DC step               step_response.csv
 *   2  Impulse               impulse_response.csv
 *   3  Float closed loop     closed_loop_float.csv
 *   4  Fixed-point closed    closed_loop_fixed.csv
 *   5  Frequency sweep       freq_response.csv
 */

static const int    FRAME_SIZE = 128;
static const double FS_ADC     = 125.0e6;
static const double FS_FILT    = FS_ADC / FRAME_SIZE;   // 976562.5 Hz
static const int    INPUT_AMPL = 50;

// ----- frequency sweep settings -----
static const int    N_SETTLE   = 128 * 8000;
static const int    N_MEASURE  = 128 * 16000;
static const int    N_TOTAL    = N_SETTLE + N_MEASURE;
static const double F_START    = 5.0e3;
static const double F_STOP     = 400.0e3;
static const int    N_FREQS    = 40;

// ----- helpers -----
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
    int n_calls       = N_TOTAL / FRAME_SIZE;
    int n_filt_settle = N_SETTLE / FRAME_SIZE;

    hls::stream<axis_t> adc_in("adc_in");
    hls::stream<axis_t> dac_out("dac_out");

    short  *out_filt = new short[n_calls];
    double *in_filt  = new double[n_calls];

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
            if (k == 0) out_filt[call] = extract_dac_ch1(out);
        }
    }

    double in_sum_sq = 0, out_sum_sq = 0, sum_sin = 0, sum_cos = 0;
    int n_meas = 0;
    for (int k = n_filt_settle; k < n_calls; k++) {
        double y = (double)out_filt[k];
        double x = in_filt[k];
        in_sum_sq  += x * x;
        out_sum_sq += y * y;
        sum_sin    += y * sin(omega_filt * k);
        sum_cos    += y * cos(omega_filt * k);
        n_meas++;
    }

    double rms_in      = sqrt(in_sum_sq / n_meas);
    double rms_out     = sqrt(out_sum_sq / n_meas);
    double gain_linear = (rms_in > 1e-12) ? rms_out / rms_in : 0;
    double gain_db     = (gain_linear > 1e-12) ? 20.0 * log10(gain_linear) : -200.0;
    double phase_deg   = atan2(sum_cos, sum_sin) * 180.0 / M_PI;

    if (verbose) {
        printf("  [%.1f Hz] rms_in=%.1f rms_out=%.1f gain=%.2f dB phase=%.1f\n",
               freq, rms_in, rms_out, gain_db, phase_deg);
    }

    delete[] out_filt;
    delete[] in_filt;
    return {gain_db, phase_deg, (double)rms_out, rms_out};
}

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

int main()
{
    printf("=== dy_cavity_relocker_2 testbench ===\n");
    printf("  ADC rate:    %.2f MHz\n", FS_ADC / 1e6);
    printf("  Filter rate: %.4f kHz  (128x decim)\n", FS_FILT / 1e3);
    printf("  Nyquist:     %.2f kHz\n\n", FS_FILT / 2e3);

    printf("--- Test 1: DC step (input = %d) ---\n", INPUT_AMPL);
    {
        hls::stream<axis_t> adc_in("dc_in");
        hls::stream<axis_t> dac_out("dc_out");
        FILE *f = fopen("step_response.csv", "w");
        fprintf(f, "frame,input,output_ch1,output_ch2\n");

        const int n_frames = 200;
        short last_ch1 = 0;
        for (int call = 0; call < n_frames; call++) {
            for (int k = 0; k < FRAME_SIZE; k++)
                adc_in.write(make_adc_word(INPUT_AMPL, 0));
            dy_cavity_relocker_2(adc_in, dac_out, false, 0, 0);

            short y1 = 0, y2 = 0;
            for (int k = 0; k < FRAME_SIZE; k++) {
                axis_t out = dac_out.read();
                if (k == 0) { y1 = extract_dac_ch1(out); y2 = extract_dac_ch2(out); }
            }
            fprintf(f, "%d,%d,%d,%d\n", call, INPUT_AMPL, y1, y2);
            last_ch1 = y1;
        }
        fclose(f);
        printf("  Final ch1 = %d\n\n", last_ch1);
    }

    printf("--- Test 2: Impulse (amp=%d at n=0) ---\n", INPUT_AMPL);
    {
        hls::stream<axis_t> adc_in("imp_in");
        hls::stream<axis_t> dac_out("imp_out");
        FILE *f = fopen("impulse_response.csv", "w");
        fprintf(f, "frame,input,output_ch1\n");

        const int n_frames = 300;
        int nonzero = 0; short maxabs = 0;
        for (int call = 0; call < n_frames; call++) {
            for (int k = 0; k < FRAME_SIZE; k++) {
                int n = call * FRAME_SIZE + k;
                short imp = (n == 0) ? INPUT_AMPL : 0;
                adc_in.write(make_adc_word(imp, 0));
            }
            dy_cavity_relocker_2(adc_in, dac_out, false, 0, 0);

            short y1 = 0;
            for (int k = 0; k < FRAME_SIZE; k++) {
                axis_t out = dac_out.read();
                if (k == 0) y1 = extract_dac_ch1(out);
            }
            fprintf(f, "%d,%d,%d\n", call, (call == 0) ? INPUT_AMPL : 0, y1);
            if (y1 != 0) nonzero++;
            if (abs(y1) > maxabs) maxabs = abs(y1);
        }
        fclose(f);
        printf("  Nonzero frames: %d/%d   max|out| = %d\n", nonzero, n_frames, maxabs);
        if (nonzero == 0) printf("  ERROR: no output from impulse!\n");
        printf("\n");
    }

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

    printf("--- Test 4: Fixed-point closed loop ---\n");
    {
        HinfFilter  controller;
        PlantSimulator plant;
        FILE *f = fopen("closed_loop_fixed.csv", "w");
        fprintf(f, "sample,error,effort\n");

        const int n_samples = 5000;
        const double disturbance = 50.0;
        double y = 0;
        for (int k = 0; k < n_samples; k++) {
            double e = disturbance - y;
            short e_short = (short)round(fmax(-32768.0, fmin(32767.0, e)));
            sig_t u_in;
            u_in.range() = e_short;
            sig_t y_ctrl = controller.process(u_in);
            short effort = (short)y_ctrl.range();
            y = plant.process((double)effort);
            fprintf(f, "%d,%.6f,%d\n", k, e, effort);
        }
        fclose(f);
        printf("  Final: error=%.3f  y=%.3f\n\n", disturbance - y, y);
    }

    printf("--- Test 5: Frequency sweep (%d pts, %.0f Hz - %.0f kHz) ---\n",
           N_FREQS, F_START, F_STOP / 1e3);
    FILE *csv = fopen("freq_response.csv", "w");
    fprintf(csv, "freq_hz,freq_mhz,gain_db,phase_deg,peak_out,rms_out\n");
    printf("%-14s %-12s %-12s %-12s\n", "Freq (kHz)", "Gain (dB)", "RMS out", "Phase");
    printf("------------------------------------------------------\n");
    for (int fi = 0; fi < N_FREQS; fi++) {
        double freq = F_START * pow(F_STOP / F_START, (double)fi / (N_FREQS - 1));
        MeasResult r = measure_freq(freq, false);
        printf("%-14.4f %-12.2f %-12.1f %-12.1f\n",
               freq / 1e3, r.gain_db, r.rms_out, r.phase_deg);
        fprintf(csv, "%.6f,%.6f,%.2f,%.2f,%.1f,%.1f\n",
                freq, freq / 1e6, r.gain_db, r.phase_deg, r.peak_out, r.rms_out);
    }
    fclose(csv);

    printf("\nCSVs: step_response.csv, impulse_response.csv, "
           "closed_loop_float.csv, closed_loop_fixed.csv, freq_response.csv\n");
    return 0;
}
