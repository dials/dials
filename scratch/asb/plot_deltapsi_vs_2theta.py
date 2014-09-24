from __future__ import division
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_experiments
import iotbx.phil
from scitbx.matrix import col
import math

"""
Script to show the deltapsi vs. 2theta plot.
Example usage:
libtbx.python plot_deltapsi_vs_2theta.py reflections.pickle experiments.json
"""

master_phil_scope = iotbx.phil.parse("""
""")

parser = OptionParser(
  phil=master_phil_scope,
  read_reflections=True,
  read_datablocks=True,
  read_experiments=True,
  check_format=False)

params, options = parser.parse_args(show_diff_phil=True)
experiments = flatten_experiments(params.input.experiments)
reflections = flatten_reflections(params.input.reflections)

beam = experiments[0].beam
s0 = col(beam.get_s0())

two_theta = [col(i).angle(s0, deg=True) for i in reflections[0]['s1']]

if 'delpsical.rad' in reflections[0]:
  delta_psi = [i*180/math.pi for i in reflections[0]['delpsical.rad']]
else:
  # need to calculate it manually. see algoirhtms/spot_prediction/stills_ray_predictor.h
  delta_psi = []
  crystal = experiments[0].crystal
  ub = crystal.get_A()
  unit_s0 = s0.normalize()
  for h in reflections[0]['miller_index']:
    # Calculate the reciprocal space vector and required unit vectors
    q = ub * h
    e1 = q.cross(unit_s0).normalize()
    c0 = unit_s0.cross(e1).normalize()

    # Calculate the vector rotated to the Ewald sphere
    qq = q.length_sq()
    wavelength = 1. / s0.length()
    a = 0.5 * qq * wavelength
    tmp = qq - a*a
    assert tmp > 0.0
    b = math.sqrt(tmp)
    r = -1.0 * a * unit_s0 + b * c0

    # Calculate delpsi value
    q0 = q.normalize()
    q1 = q0.cross(e1).normalize()
    delpsi = -1.0 * math.atan2(r.dot(q1), r.dot(q0))
    delta_psi.append(delpsi*180/math.pi)

from matplotlib import pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(111,aspect='auto')
plt.scatter(two_theta,delta_psi,c='blue',linewidth=0)
plt.show()
