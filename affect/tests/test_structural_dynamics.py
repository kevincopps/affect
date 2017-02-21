import sys
import numpy as np
sys.path.append("../affect")
import exodus
from nose.tools import raises, assert_equal

class TestStructuralDynamics():

    e = None  # database for most of the test functions

    def setUp(self):

        exodus.debug_messages(exodus.VERBOSE|exodus.DEBUG)
        base = "./SRS-FRF-example/model/1/"
        file = "p1f-out.h"
        path = base + file
        self.e = exodus.Database(path)
        print str(type(self.e))

    def tearDown(self):
        pass

    def test_frf(self):
        globals = self.e.globals
        times = globals.times()
        print 'number of timesteps {}'.format(len(times))
        nodal = self.e.nodal
        names = nodal.variable_names()
        acceleration_index = names.index('AccZ')
        force_index = names.index('AForceZ')

        a_stack = list()
        for t in xrange(len(times)):
            a_stack.append(nodal.variable(acceleration_index, t))
        a_matrix = np.hstack(a_stack)
        a_fft = np.fft.fft(a_matrix)
        print a_fft

        f_stack = list()
        for t in xrange(len(times)):
            f_stack.append(nodal.variable(force_index, t))
        f_matrix = np.hstack(f_stack)


        # f = fields['AccZ']
        # print f
        # acceleration = globals.field(f, node, 0, -1) # names.index('AccZ') + 1   #use 3
        # print acceleration
        #
        # f = fields['AForceZ']
        # print f
        # force = globals.field(f, node, 0, -1) # names.index('AForceZ')+1  #use 9
        # print force



# t=fexo.Time;
# az=fexo.NodalVars(3).Data(1,:)';
# fz=fexo.NodalVars(9).Data(2,:)';
# [freq,H]=find_frf(t,az,fz);
# figure,semilogy(freq,abs(H)), title('FRF')

# % To run the shock response spectra code do:
# fexo=exo_get('../model/1/p1f-out.h')
# t=fexo.Time;
# az=fexo.NodalVars(3).Data(1,:)';
# fz=fexo.NodalVars(9).Data(2,:)';
# sr=1/diff(t(1:2))
# itype=9
# fn=1:1:2000;  %logspace(log10(sr/1e4),log10(sr/4),50);
# [s,fn]=shspec(az(1:4096),fn,0.05,sr,itype);
# figure, plot(fn,s), title('SRS')