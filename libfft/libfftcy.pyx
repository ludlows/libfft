# distutils: language=c
#cython: language_level=3
# 2022-12-11
# https://www.github.com/ludlows/

cimport cython

cdef extern from "ludlows_libfft.h":
    cdef struct fft_ctx_double:
        unsigned long n_fft
        unsigned long * buffer_info
        unsigned long * bit_swap
        double * cos_sin_info
        double * x

    ctypedef fft_ctx_double FFT_CTX

    cdef FFT_CTX init_fft_ctx_double(unsigned long n)

    cdef void free_fft_ctx_double(FFT_CTX* pointer)

    cdef void fft_double_inplace(FFT_CTX* ctx_ptr)

    cdef void ifft_double_inplace(FFT_CTX* ctx_ptr)


cdef class FFTContext:
    cdef FFT_CTX _ctx
    cdef unsigned long _n_fft
    cdef double * _x
    def __cinit__(self, n_fft):
        if n_fft <= 1:
            raise ValueError("n_fft is too small!")
        if ((n_fft-1)&n_fft) != 0:
            raise ValueError("n_fft should equal power K of 2!")
        cdef FFT_CTX ctx
        ctx = init_fft_ctx_double(n_fft)
        if ctx.buffer_info == NULL:
            raise MemoryError("memory allocated failed")
        if ctx.bit_swap == NULL:
            raise MemoryError("memory allocated failed")
        if ctx.bit_swap == NULL:
            raise MemoryError("memory allocated failed")
        if ctx.x == NULL:
            raise MemoryError("memory allocated failed")
        self._n_fft = n_fft
        self._ctx = ctx
        self._x = ctx.x

    def __dealloc__(self):
        free_fft_ctx_double(&self._ctx)

    cpdef tuple fft(self, real_arr, image_arr):
        # check length
        cdef unsigned long arr_length = len(real_arr)
        cdef unsigned long n_fft = self._n_fft
        if arr_length > n_fft:
            raise ValueError("input is too long!")
        # copy inputs
        cdef double * x
        x = self._x
        cdef unsigned long i
        cdef unsigned long index = 0;
        for i in range(arr_length):
            x[index] = real_arr[i]
            index += 1
            x[index] = image_arr[i]
            index += 1
        # clean remain x to zeros
        while index < (n_fft << 1):
            x[index] = 0
            index += 1
        fft_double_inplace(&self._ctx)
        # copy outputs
        real_spec = []
        image_spec = []
        index = 0
        for i in range(n_fft):
            real_spec.append(x[index])
            index += 1
            image_spec.append(x[index])
            index += 1
        return (real_spec, image_spec)

    cpdef tuple ifft(self, real_arr, image_arr):
        cdef unsigned long arr_length = len(real_arr)
        cdef unsigned long n_fft = self._n_fft
        cdef double * x
        x = self._x
        cdef unsigned long i
        cdef unsigned long index = 0
        for i in range(arr_length):
            x[index] = real_arr[i]
            index += 1
            x[index] = image_arr[i]
            index += 1
        # clean remain x to zeros
        while index < (n_fft << 1):
            x[index] = 0
            index += 1
        ifft_double_inplace(&self._ctx)
        # copy outputs
        real_spec = []
        image_spec = []
        index = 0
        for i in range(n_fft):
            real_spec.append(x[index])
            index += 1
            image_spec.append(x[index])
            index += 1
        return (real_spec, image_spec)
