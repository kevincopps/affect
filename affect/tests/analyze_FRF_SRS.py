#!/usr/bin/env python

import sys
import math
import numpy as np
import numexpr as ne

sys.path.append("../affect")
import exodus
import dynamics


def frequency_response_function(path):

    with exodus.DatabaseFile(path) as e:

        times = e.globals.times()
        num_times = times.size

        names = e.nodal.variable_names()
        id_acceleration = names.index('AccZ')
        id_force = names.index('AForceZ')

        az = e.nodal.variable_at_times(id_acceleration, 0, 0, num_times)
        fz = e.nodal.variable_at_times(id_force, 1, 0, num_times)

    print "Finding frequency response function..."

    # dt = times[1] - times[0]
    # df = 1.0 / (float(num_times) * dt)
    #
    # num_times_2 = int(math.floor(num_times / 2))
    #
    # print "Computing fast fourier transform for acceleration..."
    # h_resp = np.fft.fft(az)[0:num_times_2+1] * dt
    #
    # print "Computing fast fourier transform for force..."
    # h_in = np.fft.fft(fz)[0:num_times_2+1] * dt
    #
    # print "Calculating Solution..."
    # #H = np.zeros((len(h_in)), dtype=np.double)
    # #idx = np.nonzero(h_in)  # only the non-zero entries
    # #H[idx] = np.abs(h_resp[idx] / h_in[idx])
    # H = ne.evaluate("where(h_in != 0, real(abs(h_resp / h_in)), 0.0)")
    # print H.size, H
    #
    # freq = np.linspace(0.0, float(num_times_2)*df, num_times_2+1)

    frequency, h_transfer = dynamics.frf(fz, az, times)

    # print "Finding peaks"
    # from scipy.signal import argrelmax
    # peaks = argrelmax(H)
    # for i,j in enumerate(peaks):
    #     print i, freq[j], H[j]

    '''
    print "Outputting results..."
    print 'Plotting...'
    from matplotlib.widgets import Cursor
    import matplotlib.pyplot as plt

    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in xrange(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)
    plt.figure(figsize=(12,9))
    plt.plot(frequency, h_transfer, ConsoleCode=tableau20[0])
    plt.yscale('log')
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # Make sure your axis ticks are large enough to be easily read.
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.xlabel("Frequency", fontsize=16)
    ##plt.savefig('FRF.png')
    plt.show('FRF.png')
    '''

    '''
    f1 = open('Freq.txt','w')
    newline = ''
    for line in freq.T: #turns freq into column vector
    for char in line[0][0]:
      if char != '[' or char != ']' or char != ' ': #formatting data
        newline = newline + char
    f1.write(str(newline)+'\n')
    newline = ''
    f1.close()
    '''

    '''
    f = open('h_resp.txt','w')
    for line in h_resp:
    f.write(str(line[0][0])+'\n')
    f.close()

    f = open('h_in.txt','w')
    for line in h_in:
    f.write(str(line[0][0])+'\n')
    f.close()
    '''


def main():
    base = "./SRS-FRF-example/model/1/"
    file = "p1f-out.h"
    path = base + file
    frequency_response_function(path)
  

if __name__== '__main__':
    main()
