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


    print "______________________________________________________________________________"
    print "p  =", p
    print "tpl_n  =", tpl_n
    print "tpl_s0 =", tpl_s0
    print "tpl_m2 =", tpl_m2
    print "tpl_s1 =", tpl_s1
    print
    print "I(R) =", reflection.intensity
    import numpy as np

    n = np.array(tpl_n)
    s0 = np.array(tpl_s0)
    u = np.array(tpl_m2)
    s = np.array(tpl_s1)

    print "np vect(np_n ) =", n
    print "np vect(np_s0) =", s0
    print "np vect(np_m2) =", u
    print "np vect(np_s1) =", s
    print
    print "cross_prod(s0,s) =", np.cross(s0, s)

    L_f = abs(np.dot(s, np.cross(u, s0))) / np.dot(abs(s), abs(s0))
    print "Loretz factor =", L_f
    print



    print "______________________________________________________________________________"


    reflection.intensity = reflection.intensity / 2.0
