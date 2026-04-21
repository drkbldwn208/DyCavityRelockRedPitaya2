// hinf_coeffs.h — auto-generated
// Q3.29, 7 SOS sections, fs=31.2500 MHz
// Quantization error: -15.0 dB
#ifndef HINF_COEFFS_H
#define HINF_COEFFS_H
#include <stdint.h>

#define HINF_N_SECTIONS     7
#define HINF_COEF_INT_BITS  3
#define HINF_COEF_FRAC_BITS 29
#define HINF_COEF_SCALE     536870912L

typedef struct { int32_t b0,b1,b2,a1,a2; } hinf_sos_t;

static const hinf_sos_t HINF_SOS[7] = {
    {  171823497, -340514798,  168697532,-1028454714,  492057499},  // section 0
    {  536870912,-1071585083,  534716354,-1072282190,  535411876},  // section 1
    {  536870912,-1073094441,  536223529,-1073633888,  536762982},  // section 2
    {  536870912,-1073736924,  536866702,-1073736619,  536866405},  // section 3
    {  536870912,-1073737748,  536867194,-1073739210,  536868661},  // section 4
    {  536870912,-1073740253,  536869482,-1073739719,  536868949},  // section 5
    {  536870912,   33294321, -503576591,       8193, -536861800}  // section 6
};

#endif