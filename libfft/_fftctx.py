from .libfftcy import FFTContext


class FFTCtx:
    """
    FFT (Fast Fourier Transform) Context
    """

    def __init__(self, n_fft):
        if n_fft < 2:
            raise ValueError("n_fft is too small")
        if (n_fft - 1) & n_fft != 0:
            raise ValueError("n_fft must be a power of 2!")
        self._n_fft = n_fft
        self._ctx = None

    def open(self):
        if self._ctx is None:
            self._ctx = FFTContext(self._n_fft)

    def close(self):
        if self._ctx is not None:
            del self._ctx
            self._ctx = None

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def fft(self, arr_list):
        """
        compute fft spectrum
        arr_list: signal x(n)
        return: [spectrum1, spectrum2, ...], with length of n_fft
        """
        if len(arr_list) > self._n_fft:
            raise ValueError("input list is too long!")
        real_in = [(1.0j * v).imag for v in arr_list]
        image_in = [(-1.0j * v).real for v in arr_list]
        real_out, image_out = self._ctx.fft(real_in, image_in)
        return [r+i*1j for r, i in zip(real_out, image_out)]

    def ifft(self, spec_list):
        """
        inverse fast fourier transform
        spec_list: spectrum X(K)
        return: [x(0), x(1), ...], with length of n_fft
        """
        if len(spec_list) > self._n_fft:
            raise ValueError("input list is too long!")
        real_in = [(1.0j * v).imag for v in spec_list]
        image_in = [(-1.0j * v).real for v in spec_list]
        real_out, image_out = self._ctx.ifft(real_in, image_in)
        return [r+i*1j for r, i in zip(real_out, image_out)]
