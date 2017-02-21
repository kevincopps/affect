#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Kevin Copps'

import sys
import math
import numpy as np
import numexpr as ne

def frf(f_in, x_response, times):
    """
    Return the transfer function frequency response function of a single input, single output system. using fast Fourier
    transforms.

    The frequency response function is a transfer function. The frequency response function is used often to express
    the structural response to an applied force as a function of frequency. Frequency response functions are complex
    functions, with real and imaginary components.

    In the frequency domain, the structural response X(ω) is expressed as the product of the frequency response
    function H(ω) and the input or applied force F(ω). The response X(ω) may be in terms of displacement, velocity, or
    acceleration.

        X(ω) = H(ω)⋅F(ω)

    :param f_in: input force values at sampled times
    :type f_in: numpy.ndarray with size equal to times.size
    :param x_response: response function values at sampled times
    :type x_response: numpy.ndarray with size equal to times.size
    :param times: time values
    :type times: numpy.ndarray
    :return: frequency, h_transfer
    :rtype frequency: numpy array of frequencies ω, with size equal to times.size/2 + 1
    :rtype h_transfer: numpy array of the magnitude (real values) of transfer function H(ω)
    """
    num_times = times.size
    if x_response.size != num_times or f_in.size != num_times:
        raise ValueError('input array parameters have different sizes')

    dt = times[1] - times[0]
    df = 1.0 / (float(num_times) * dt)
    m = int(math.floor(num_times / 2)) + 1

    x_response_f = np.fft.fft(x_response)[0:m] * dt
    f_in_f = np.fft.fft(f_in)[0:m] * dt

    # division using the numexpr package, avoiding divide by zero
    h_transfer = ne.evaluate("where(f_in_f != 0, real(abs(x_response_f / f_in_f)), 0.0)")

    frequency = np.linspace(0.0, float(m-1)*df, m)

    return frequency, h_transfer
