// hinf_coeffs.h — auto-generated
// Q3.29, 7 SOS sections, fs=0.9766 MHz
// Quantization error: -62.9 dB
#ifndef HINF_COEFFS_H
#define HINF_COEFFS_H
#include <stdint.h>

#define HINF_N_SECTIONS     7
#define HINF_COEF_INT_BITS  3
#define HINF_COEF_FRAC_BITS 29
#define HINF_COEF_SCALE     536870912L

typedef struct { int32_t b0,b1,b2,a1,a2; } hinf_sos_t;

static const hinf_sos_t HINF_SOS[7] = {
    {    8319172,  -15600946,    7314361, -317072051,          0},  // section 0
    {  536870912,    6366047,          0,    4462267, -524265377},  // section 1
    {  536870912, -824672835,  299793433,-1070314330,  533449522},  // section 2
    {  536870912,-1072900260,  536736259,-1072882749,  536726691},  // section 3
    {  536870912,-1073257134,  536752120,-1073298071,  536798815},  // section 4
    {  536870912,-1073551492,  536825045,-1073533881,  536808423},  // section 5
    {  536870912,  604202069,   67331157,    1374714, -535485817}  // section 6
};

#endif