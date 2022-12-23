#ifndef LUDLOWS_LIBFFT_H
#define LUDLOWS_LIBFFT_H

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/*
fast fourier transform using 2-base method
https://github.com/ludlows/libfft
complex signal x(n) with length of N is represented by an array of 2*N
array[0] is the real part of x(0)
array[1] is the complex part of x(0)
*/

struct fft_ctx_double {
    unsigned long n_fft;
    unsigned long * buffer_info;   // n_fft >> 1
    unsigned long * bit_swap;      // n_fft
    double * cos_sin_info;         // n_fft
    double * x;                    // n_fft << 1, real and image parts
};

typedef struct fft_ctx_double FFT_CTX;

void free_fft_ctx_double(FFT_CTX* pointer) {
    if (pointer == NULL) return;
    free(pointer->bit_swap);
    pointer->bit_swap = NULL;
    free(pointer->buffer_info);
    pointer->buffer_info = NULL;
    free(pointer->cos_sin_info);
    pointer->cos_sin_info = NULL;
    free(pointer->x);
    pointer->x = NULL;
}

FFT_CTX init_fft_ctx_double(unsigned long n) {
    FFT_CTX ctx;
    ctx.n_fft = n;
    ctx.buffer_info = NULL;
    ctx.cos_sin_info = NULL;
    ctx.bit_swap = NULL;
    ctx.x = NULL;

    void * p_ctx_buffer_info = malloc(sizeof(unsigned long) * (n >> 1));
    if (p_ctx_buffer_info == NULL) {
        printf("allocate memory of buffer_info failed!\n");
        free_fft_ctx_double(&ctx);
        return ctx;
    }
    ctx.buffer_info = (unsigned long *)p_ctx_buffer_info;
    void * p_ctx_cos_sin_info = malloc(sizeof(double) * n);
    if (p_ctx_cos_sin_info == NULL) {
        printf("allocate memory of cos_sin_info  failed!\n");
        free_fft_ctx_double(&ctx);
        return ctx;
    }
    ctx.cos_sin_info = (double*)p_ctx_cos_sin_info;
    void * p_ctx_bit_swap = malloc(sizeof(unsigned long) * n);
    if (p_ctx_bit_swap == NULL) {
        printf("allocate memory of bit_swap failed!\n");
        free_fft_ctx_double(&ctx);
        return ctx;
    }
    ctx.bit_swap = (unsigned long *)p_ctx_bit_swap;
    double * p_ctx_x = (double*)malloc(sizeof(double) * (n << 1));
    if (p_ctx_x == NULL) {
        printf("allocate memory of x failed!\n");
        free_fft_ctx_double(&ctx);
        return ctx;
    }
    ctx.x = p_ctx_x;
    unsigned long len = n << 1;
    for (unsigned long i = 0; i < len; i++) {
        p_ctx_x[i] = 0.0f;
    }
    return ctx;
}

/*
n_point fast fourier transform
n_point: integer , power of 2
x: length 2 * n_point, signal x(n) whose length is n_point only
x[0]: signal x(0)'s real part
x[1]: signal x(0)'s image part
*/
void fft_double_inplace(FFT_CTX* ctx_ptr) {
    if (ctx_ptr == NULL) {
        printf("Error! fft context is NULL.\n");
        return;
    }
    unsigned long n_point = ctx_ptr->n_fft;
    if (n_point <= 1) {
        printf("Error! nFFT <= 1.\n");
        return;
    }
    double * x = ctx_ptr->x;
    unsigned long * buffer_info = ctx_ptr->buffer_info;
    double * cos_sin_info = ctx_ptr->cos_sin_info;
    unsigned long * bit_swap = ctx_ptr->bit_swap;
    register double theta, index1_real, index2_real, index1_image, index2_image;
    unsigned long k = 0;
    // initialize cos_sin_info
    for (unsigned long i = 0, k = 0; i < (n_point >> 1); i++) {
        theta = 6.283185307179586f * i / n_point;
        cos_sin_info[k++] = cos(theta);
        cos_sin_info[k++] = sin(theta);
    }
    // initialize buffer_info
    buffer_info[0] = 0;
    k = n_point >> 2;
    unsigned long len = 1;
    while (k >= 1) {
        for (unsigned long i = 0; i < len; i++) {
            buffer_info[i + len] = buffer_info[i] + k;
        }
        len <<= 1;
        k >>= 1;
    }
    unsigned long index1 = 0;
    unsigned long index2 = 0;
    unsigned long n_cycle = 0;
    unsigned long step = 0;
    // start
    for (n_cycle = 1, step = n_point >> 1; n_cycle < n_point; n_cycle <<= 1, step >>= 1) { //log(N)
        index1 = 0;
        index2 = step << 1;

        for (unsigned long i = 0; i < n_cycle; i++) {
            unsigned long index_cos_sin_info = buffer_info[i] << 1;
            double real_part = cos_sin_info[index_cos_sin_info];
            double image_part = cos_sin_info[index_cos_sin_info + 1];
            for (unsigned long j = 0; j < step; j++) {
                index1_real = x[index1];
                index1_image = x[index1 + 1];
                index2_real = x[index2];
                index2_image = x[index2 + 1];
                x[index1++] = index1_real + real_part * index2_real + image_part * index2_image;
                x[index1++] = index1_image - image_part * index2_real + real_part * index2_image;
                x[index2++] = index1_real - real_part * index2_real - image_part * index2_image;
                x[index2++] = index1_image + image_part * index2_real - real_part * index2_image;
            }
            index1 = index2;
            index2 = index1 + (step << 1);
        }
    }
    n_cycle = n_point >> 1;
    for (unsigned long i = 0; i < n_cycle; i++) {
        bit_swap[i] = buffer_info[i] << 1;
        bit_swap[i + n_cycle] = 1 + bit_swap[i];
    }
    unsigned long s = 0;

    for (unsigned long i = 0; i < n_point; i++) {
        if (bit_swap[i] != i) {
            s = bit_swap[i];
            bit_swap[s] = s;
            index1 = i << 1;
            index2 = s << 1;
            index1_real = x[index1];
            x[index1++] = x[index2];
            x[index2++] = index1_real;
            index1_real = x[index1];
            x[index1] = x[index2];
            x[index2] = index1_real;
        }
    }
}

/*
n_point inverse fast fourier transform
n_point: integer , power of 2
x: length 2 * n_point, spectrum X(n) whose length is n_point only
x[0]: spectrum X(0)'s real part
x[1]: spectrum X(0)'s image part
*/
void ifft_double_inplace(FFT_CTX* ctx_ptr) {
    if (ctx_ptr == NULL) {
        printf("Error! fft context is NULL.\n");
        return;
    }
    unsigned long n_point = ctx_ptr->n_fft;
    if (n_point <= 1) {
        printf("Error! nFFT <= 1\n");
        return;
    }
    double * x = ctx_ptr->x;
    unsigned long * buffer_info = ctx_ptr->buffer_info;
    double * cos_sin_info = ctx_ptr->cos_sin_info;
    unsigned long * bit_swap = ctx_ptr->bit_swap;
    register double theta, index1_real, index2_real, index1_image, index2_image;
    unsigned long k = 0;
    // initialize cos_sin_info
    for (unsigned long i = 0, k = 0; i < (n_point >> 1); i++) {
        theta = 6.283185307179586f * i / n_point;
        cos_sin_info[k++] = cos(theta);
        cos_sin_info[k++] = sin(theta);
    }
    // initialize buffer_info
    buffer_info[0] = 0;
    k = n_point >> 2;
    unsigned long len = 1;
    while (k >= 1) {
        for (unsigned long i = 0; i < len; i++) {
            buffer_info[i + len] = buffer_info[i] + k;
        }
        len <<= 1;
        k >>= 1;
    }
    unsigned long index1 = 0;
    unsigned long index2 = 0;
    unsigned long n_cycle = 0;
    unsigned long step = 0;
    // start 
    for (n_cycle = 1, step = n_point >> 1; n_cycle < n_point; n_cycle <<= 1, step >>= 1) { //log(N)
        index1 = 0;
        index2 = step << 1;

        for (unsigned long i = 0; i < n_cycle; i++) {
            unsigned long index_cos_sin_info = buffer_info[i] << 1;
            double real_part = cos_sin_info[index_cos_sin_info];
            double image_part = cos_sin_info[index_cos_sin_info + 1];
            for (unsigned long j = 0; j < step; j++) {
                index1_real = x[index1];
                index1_image = x[index1 + 1];
                index2_real = x[index2];
                index2_image = x[index2 + 1];
                x[index1++] = index1_real + real_part * index2_real - image_part * index2_image;
                x[index1++] = index1_image + image_part * index2_real + real_part * index2_image;
                x[index2++] = index1_real - real_part * index2_real + image_part * index2_image;
                x[index2++] = index1_image - image_part * index2_real - real_part * index2_image;
            }
            index1 = index2;
            index2 = index1 + (step << 1);
        }
    }
    n_cycle = n_point >> 1;
    for (unsigned long i = 0; i < n_cycle; i++) {
        bit_swap[i] = buffer_info[i] << 1;
        bit_swap[i + n_cycle] = 1 + bit_swap[i];
    }
    unsigned long s = 0;

    for (unsigned long i = 0; i < n_point; i++) {
        if (bit_swap[i] != i) {
            s = bit_swap[i];
            bit_swap[s] = s;
            index1 = i << 1;
            index2 = s << 1;
            index1_real = x[index1];
            x[index1++] = x[index2];
            x[index2++] = index1_real;
            index1_real = x[index1];
            x[index1] = x[index2];
            x[index2] = index1_real;
        }
    }
    unsigned long N = ctx_ptr->n_fft;
    len = (N << 1);
    for (unsigned long i = 0; i < len; i++) { // normalization
        x[i] /= N;
    }
}
#endif // LUDLOWS_LIBFFT_H
