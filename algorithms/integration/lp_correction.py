def correct_intensity(sweep, crystal, reflections):
    for ref in reflections:
        if ref.status == 0:
            LP_calculations(sweep, crystal, ref)
    return reflections

def LP_calculations(sweep, crystal, reflection):
    tpl_n = sweep.get_beam().get_polarization_normal()
    tpl_s0 = sweep.get_beam().get_s0()
    tpl_m2 = sweep.get_goniometer().get_rotation_axis()
    tpl_s1 = reflection.beam_vector
    p = sweep.get_beam().get_polarization_fraction()
    from scitbx import matrix
    #import math
    #import numpy as np

    n = matrix.col(tpl_n)
    s0 = matrix.col(tpl_s0)
    u = matrix.col(tpl_m2)
    s = matrix.col(tpl_s1)

    L_f = abs(s.dot(u.cross(s0))) / (s.length() * s0.length())

    P_f = (1 - 2 * p) * (1 - (n.dot(s) / s.length()) ** 2.0) + \
     p * (1 + (s.dot(s0) / (s.length() * s0.length())) ** 2.0)

    reflection.intensity = reflection.intensity * L_f * P_f

