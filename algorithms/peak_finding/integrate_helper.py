from scitbx.array_family import flex
from dials.algorithms.peak_finding import model_2d, measure_2d_angl
from matplotlib import pyplot as plt
for int_ref_ang in range(30):
    print "_________________________________________"
    ref_ang = float(int_ref_ang) / 30.0

    ref2d = model_2d(30, 30, 2, 5, ref_ang, 25, 0.5)
    dat2d_ref = ref2d.as_numpy_array()
    for pr_line in dat2d_ref:
        print pr_line

    ang = measure_2d_angl(ref2d, 15, 15)
    print "ref_ang =", ref_ang
    print "ang =", ang

    new_ref2d = model_2d(30, 30, 2, 5, ang, 25, 0.5)
    dat2d_paint = new_ref2d.as_numpy_array()
    plt.imshow(dat2d_paint , interpolation = "nearest")
    plt.show()
