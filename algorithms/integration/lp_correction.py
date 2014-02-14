from __future__ import division
def correct_intensity(sweep, crystal, reflections):
  from dials.util.command_line import Command
  from dials.array_family import flex
  Command.start('Performing LP-correction')
  s1 = reflections['s1']
  I = reflections['intensity.raw.value']
  varI = reflections['intensity.raw.variance']
  corrected_I = flex.double(len(I))
  corrected_varI = flex.double(len(I))
  for i in range(len(reflections)):
    lp = LP_calculations(sweep, crystal, s1[i])
    corrected_I[i] = I[i] * lp
    corrected_varI[i] = varI[i] * lp
  reflections['intensity.cor.value'] = corrected_I
  reflections['intensity.cor.variance'] = corrected_varI
  Command.end('Performed LP-correction on {0} reflections'.format(
    len(reflections)))
  return reflections

def LP_calculations(sweep, crystal, s1):
  '''See Kabsch, J. Appl. Cryst 1988 21 916-924.'''

  tpl_n = sweep.get_beam().get_polarization_normal()
  tpl_s0 = sweep.get_beam().get_s0()
  tpl_m2 = sweep.get_goniometer().get_rotation_axis()
  tpl_s1 = s1
  p = sweep.get_beam().get_polarization_fraction()

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

