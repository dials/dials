from __future__ import absolute_import, division, print_function

import os
import six.moves.cPickle as pickle


def test_run(dials_regression):
    from dials.algorithms.background.simple import MosflmOutlierRejector
    from dials.algorithms.background.simple import Linear2dModeller

    # The directory path
    path = os.path.join(
        dials_regression, "integration_test_data", "i04-weak-data", "jmp_mosflm_test"
    )

    # The input files
    reflection_filename = os.path.join(path, "mosflm_reflections.pickle")
    shoebox_filename = os.path.join(path, "shoeboxes.pickle")
    fraction = 1.0
    n_sigma = 4
    outlier_rejector = MosflmOutlierRejector(fraction, n_sigma)
    linear_modeller = Linear2dModeller()

    from dials.array_family import flex
    from math import sqrt
    from dials.algorithms.shoebox import MaskCode

    print(shoebox_filename)
    # Read the data
    rtable = flex.reflection_table.from_pickle(reflection_filename)
    with open(shoebox_filename, "rb") as fh:
        shoeboxes, masks = pickle.load(fh)
    assert len(rtable) == len(shoeboxes)
    assert len(rtable) == len(masks)

    # Compute the background for each reflection and check against the values
    # read from the mosflm.lp file. Currently this fails for 1 strange
    # reflection whose pixel values in the mosflm file do not match those
    # extracted from the images.
    count = 0
    VAR1 = []
    VAR2 = []
    DIFF = []
    for i in range(len(rtable)):
        I = rtable[i]["intensity.sum.value"]
        Ivar = rtable[i]["intensity.sum.variance"]
        data = shoeboxes[i].as_double()
        mask = masks[i]
        try:
            assert len(data.all()) == 2
            assert len(mask.all()) == 2
            data.reshape(flex.grid(1, *data.all()))
            mask.reshape(flex.grid(1, *mask.all()))
            outlier_rejector(data, mask)
            mask2 = (mask.as_1d() & int(MaskCode.BackgroundUsed)) != 0
            mask2.reshape(flex.grid(*mask.all()))
            model = linear_modeller.create(data, mask2)
        except Exception:
            count += 1
            raise
        assert len(model.params()) == 3
        hy = data.all()[1] // 2
        hx = data.all()[2] // 2
        c1 = model.params()[0]
        a1 = model.params()[1]
        b1 = model.params()[2]
        c3 = c1 + a1 * (0.5 + hx) + b1 * (0.5 + hy)
        a2 = rtable[i]["background"][0]
        b2 = rtable[i]["background"][1]
        c2 = rtable[i]["background"][2]

        try:
            assert abs(a1 - b2) < 0.01
            assert abs(b1 + a2) < 0.01
            assert abs(c3 - c2) < 0.1
        except Exception:
            count += 1
            continue

        background = data.as_double()
        # hy = background.all()[1] // 2
        # hx = background.all()[2] // 2
        for jj in range(background.all()[1]):
            for ii in range(background.all()[2]):
                # x = ii - hx
                # y = jj - hy
                x = ii + 0.5
                y = jj + 0.5
                background[0, jj, ii] = a1 * x + b1 * y + c1

        # Test the summation results. Edge reflections use profile fitted
        # intensity in MOSFLM. Therefore ignore these. There also appears to be a
        # some discrepancy with very low <= 0 reflections where an extra 0.5 is
        # added. Not sure why this is so ignore these reflections as well.
        from dials.algorithms.integration.sum import integrate_by_summation

        intensity = integrate_by_summation(data.as_double(), background, mask)
        I2 = intensity.intensity()
        Ivar2 = intensity.variance()
        I1 = I
        Ivar1 = Ivar
        if mask.count(0) == 0 and mask.count(2) == 0 and I1 > 0:
            VAR1.append(sqrt(Ivar1))
            VAR2.append(sqrt(Ivar2))
            DIFF.append(sqrt(Ivar1) - sqrt(Ivar2))
            try:
                assert abs(I1 - I2) < 1.0
                assert abs(sqrt(Ivar1) - sqrt(Ivar2)) < 1.0
            except Exception:
                count += 1
                continue

    # Only 1 should fail
    assert count == 1
