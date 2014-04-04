from __future__ import division
def correct_intensity(experiment, reflections):
  from dials.util.command_line import Command
  from dials.array_family import flex
  Command.start('Performing LP-correction')
  lp = flex.double(
    [LP_calculations(experiment, s1)
     for s1 in reflections['s1']])
  reflections['lp'] = lp
  Command.end('Performed LP-correction on {0} reflections'.format(lp))
  return lp

def LP_calculations(experiment, s1):
  '''See Kabsch, J. Appl. Cryst 1988 21 916-924.'''

  tpl_n = experiment.beam.get_polarization_normal()
  tpl_s0 = experiment.beam.get_s0()
  tpl_m2 = experiment.goniometer.get_rotation_axis()
  tpl_s1 = s1
  p = experiment.beam.get_polarization_fraction()

  # FIXME hack for testing
  # p = 0.5

  from scitbx import matrix

  n = matrix.col(tpl_n)
  s0 = matrix.col(tpl_s0)
  u = matrix.col(tpl_m2)
  s = matrix.col(tpl_s1)

  L_f = abs(s.dot(u.cross(s0))) / (s.length() * s0.length())

  P_f = (1 - 2 * p) * (1 - (n.dot(s) / s.length()) ** 2.0) + \
        p * (1 + (s.dot(s0) / (s.length() * s0.length())) ** 2.0)

  return L_f / P_f
