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
    import numpy as np

    n = np.array(tpl_n)
    s0 = np.array(tpl_s0)
    u = np.array(tpl_m2)
    s = np.array(tpl_s1)

    L_f = np.abs(np.dot(s, np.cross(u, s0))) / (modulus_3D(s) * modulus_3D(s0))

    P_f = (1 - 2 * p) * (1 - (np.dot(n, s) / modulus_3D(s)) ** 2.0)          \
      + p * (1 + (np.dot(s, s0) / (modulus_3D(s) * modulus_3D(s0))) ** 2.0)

    reflection.intensity = reflection.intensity * L_f * P_f

def modulus_3D(vect):
    import numpy as np
    mold = np.sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2])
    return mold
