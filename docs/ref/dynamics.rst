=============================================
``dynamics`` --- Structural Dynamics Analysis
=============================================

:mod:`affect.dynamics`
======================

.. automodule:: affect.dynamics

.. currentmodule:: affect.dynamics

Summary
-------

This module contains useful functions and postprocessors concerned with analyzing the behavior of physical structures
when subjected to dynamic forces. This module is useful when the applied dynamic forces result in accelerations
high enough to excite the structure's natural frequency.

Dynamic analysis can be used to find dynamic displacements, time history, and modal analysis.

.. autosummary::

   frf

Frequency Response Function
---------------------------

Frequency response is the quantitative measure of the output spectrum of a system or device in response to a stimulus,
and is used to characterize the dynamics of the system. It is a measure of magnitude and phase of the output as a
function of frequency, in comparison to the input. The frequency response function is a transfer function used to
identify the resonant frequencies, damping and mode shapes of a physical structure.

.. math::
   :nowrap:

   \begin{gather*}
      \begin{split}
         & \text{Input Force}  \\
         & \quad F(ω)
      \end{split} & \longrightarrow &
      \begin{split}
         & \text{Transfer Function}\\
         & \qquad H(ω)
      \end{split} & \longrightarrow &
      \begin{split}
         & \text{Displacement Response}\\
         & \qquad\quad X(ω)
      \end{split}
   \end{gather*}

Here, :math:`F` is the input force as a function of the frequency :math:`\omega`, and :math:`H` is the transfer
function, while :math:`X` is the displacement (or velocity or acceleration) response function.
Each function is a complex function, with real and imaginary components, which may also be represented in terms of
magnitude and phase, and thus the functions are spectral functions. For sake of computation and simplicity, we consider
each to be a Fourier transform.

Thus, in the frequency domain, the structural response X(ω) is usually expressed as the product of the frequency
response function H(ω) and the input or applied force F(ω). Usually the response X(ω) may be in terms of displacement,
velocity, or acceleration.

.. math:: X(ω) = H(ω)⋅F(ω)

.. math:: H(ω) = \frac{X(ω)}{F(ω)}

Using a frequency response function, the following can be observed:

* Resonances - Peaks indicate the presence of the natural frequencies of the structure under test
* Damping - Damping is proportional to the width of the peaks. The wider the peak, the heavier the damping
* Mode Shape – The amplitude and phase of multiple FRFs acquired to a common reference on a structure are used to
  determine the mode shape

Nomenclature:
^^^^^^^^^^^^^

Various transfer functions are useful for measuring system response and these have common names:

====================  ===================================
Quantity              Name of Frequency Response Function
====================  ===================================
displacement / force  admittance, compliance, receptance
velocity / force      mobility
acceleration / force  accelerance, inertance
force / displacement  dynamic stiffness
force / velocity      mechanical impedance
force / acceleration  apparent mass, dynamic mass
====================  ===================================


Example:
^^^^^^^^

Examine the natural frequencies in the computational results of a structural model. Find the peak values of the
`accelerance` frequency response function, where the response is acceleration in the z-direction given an input
stimulus of force in the z-direction.

.. code-block:: python
   :linenos:
   :emphasize-lines: 12

   from affect.exodus import DatabaseFile
   from affect.dynamics import frf
   from scipy.signal import argrelmax

   with DatabaseFile('./SRS-FRF-example/model/1/p1f-out.h') as e:
       times = e.globals.times()
       num_times = times.size
       node_vars = e.nodal.variables
       az = e.nodal.variable_at_times(node_vars['AccZ'], 0, 0, num_times)
       fz = e.nodal.variable_at_times(node_vars['AForceZ'], 1, 0, num_times)

   frequency, h_transfer = frf(fz, az, times)
   peaks = argrelmax(h_transfer)
   for i, j in enumerate(peaks):
       print(i, frequency[j], h_transfer[j])

.. autofunction:: frf


