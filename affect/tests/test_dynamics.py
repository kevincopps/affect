import numpy as np
from scipy.signal import argrelmax

from .. import dynamics
from .. import util


def test_frf(frf_data_dict):
    # e = edb_frf
    # times = e.globals.times()
    # print('    number of timesteps {}'.format(len(times)))
    # nodal = e.nodal
    # names = nodal.variable_names()
    # acceleration_index = names.index('AccZ')
    # force_index = names.index('AForceZ')

    util.print_function_starting()

    times = frf_data_dict['times']
    print(f'num_times = {times.size}')
    print(f'times = {times}')

    fz = frf_data_dict['force_z']
    print(f'force = {fz}')

    az = frf_data_dict['acceleration_z']
    print(f'acceleration = {az}')

    frequency, h_transfer = dynamics.frf(fz, az, times)  # get values of the frequency and FRF

    indices = argrelmax(h_transfer, order=2)  # get indices at the relative peaks of the FRF

    # Note that the return value 'indices' of argrelmax is a tuple even when data is 1 - dimensional.
    # Use the 0-th entry of the tuple.
    frequency_peaks = np.take(frequency, indices[0])
    frf_peaks = np.take(h_transfer, indices[0])

    # check if the first handful of peaks of frequency are close to what we expect
    util.print_array_info('frequency_peaks', frequency_peaks[:9])
    desired_frequency_peaks = np.array([451.0, 702.5, 1000.0, 1157.0, 1815.5, 1817.5, 2737.5, 2739.5, 3671.5])
    np.testing.assert_allclose(frequency_peaks[:9], desired_frequency_peaks, rtol=1e-3, atol=0.5, verbose=True)

    # and the peaks of the frf transfer function
    util.print_array_info('frf_peaks', frf_peaks[:6])
    desired_frf_peaks = np.array([1543.897, 3036.994, 2024.438, 481.288, 2.302227, 2.302244])
    np.testing.assert_allclose(frf_peaks[:6], desired_frf_peaks, rtol=1.5e-7, verbose=True)
