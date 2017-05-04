#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Kevin Copps'

import math
import numpy
import numexpr as ne


def frf(f_in: numpy.ndarray, x_response: numpy.ndarray, times: numpy.ndarray) -> (numpy.ndarray, numpy.ndarray):
    """
    Return the magnitude (real part) of the frequency response function of a single input/output system using 
    fast Fourier transforms.

    The frequency response function is a transfer function often used to express
    the structural response to an applied force as a function of frequency. Frequency response functions are complex
    functions, with real and imaginary components, magnitude and phase angle, respectively.

    In the frequency domain, the structural response X(ω) is expressed as the product of the frequency response
    function H(ω) and the input or applied force F(ω). The response X(ω) may be in terms of displacement, velocity, or
    acceleration.

    .. math:: X(ω) = H(ω)⋅F(ω)
        
    Args:

        f_in: input force values at the time values, with size equal to time.size
        x_response: response function values at the time values, with size equal to time.size
        times: time values
        
    Returns:
        (frequencies, magnitudes) - arrays of frequencies :math:`ω`, and magnitudes (real values) of the transfer 
        function :math:`H(ω)`, with sizes equal to ``times.size/2 + 1``

    Raises:
        ValueError: if ``x_response.size != num_times`` or ``f_in.size != num_times``
    """
    num_times = times.size
    if x_response.size != num_times or f_in.size != num_times:
        raise ValueError('input array parameters have different sizes')

    dt = times[1] - times[0]
    df = 1.0 / (float(num_times) * dt)
    m = int(math.floor(num_times / 2)) + 1

    x_response_f = numpy.fft.fft(x_response)[0:m] * dt
    f_in_f = numpy.fft.fft(f_in)[0:m] * dt

    # division using the numexpr package, avoiding divide by zero
    h_transfer = ne.evaluate("where(f_in_f != 0, real(abs(x_response_f / f_in_f)), 0.0)")

    frequency = numpy.linspace(0.0, float(m - 1) * df, m)

    return frequency, h_transfer
