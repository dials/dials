from __future__ import absolute_import, division, print_function

import itertools
import math
import random

from libtbx.phil import parse

master_phil = parse(
    """
nrefl = 0
  .type = int
shoebox_size {
  x = 10
    .type = int
  y = 10
    .type = int
  z = 10
    .type = int
}
spot_size {
  x = 1.0
    .type = float
  y = 1.0
    .type = float
  z = 1.0
    .type = float
}
spot_offset {
  x = 0.0
    .type = float
  y = 0.0
    .type = float
  z = 0.0
    .type = float
}
rs_node_size = 0.05
  .type = float
rs_window_size = 0.25
  .type = float
mask_nsigma = 3.0
  .type = float
integrated_data_file = None
  .type = path
integrated_data_file_scale = 0.0
  .type = float
counts = 0
  .type = int
background = 0
  .type = int
background_a = 0.0
  .type = float
background_b = 0.0
  .type = float
background_c = 0.0
  .type = float
background_d = 0.0
  .type = float
pixel_mask = *all static precise
  .type = choice
background_method = *xds mosflm
  .type = choice
integration_methpd = *xds mosflm
  .type = choice
output {
  over = None
    .type = path
  under = None
    .type = path
  all = None
    .type = path
}
rotation {
  axis {
    x = 0.0
      .type = float
    y = 0.0
      .type = float
    z = 0.0
      .type = float
  }

  angle = 0.0
    .type = float
}
"""
)


def random_background_plane2(sbox, a, b, c, d):
    """Draw values from Poisson distribution for each position where the mean for
    that distribition is equal to a + b * i + c * j + d * k where a, b, c, d are
    floating point values and i, j, k are the shoebox indices in directions x, y
    and z respectively."""

    from scitbx.random import variate, poisson_distribution

    dz, dy, dx = sbox.focus()

    if b == c == d == 0.0:
        g = variate(poisson_distribution(mean=a))
        for k in range(dz):
            for j in range(dy):
                for i in range(dx):
                    sbox[k, j, i] += next(g)
    else:
        for k in range(dz):
            for j in range(dy):
                for i in range(dx):
                    pixel = a + b * (i + 0.5) + c * (j + 0.5) + d * (k + 0.5)
                    g = variate(poisson_distribution(mean=pixel))
                    sbox[k, j, i] += next(g)
    return


def random_background_plane(sbox, a, b, c, d):
    """Draw values from Poisson distribution for each position where the mean for
    that distribition is equal to a + b * i + c * j + d * k where a, b, c, d are
    floating point values and i, j, k are the shoebox indices in directions x, y
    and z respectively."""

    from scitbx.random import variate, poisson_distribution

    dz, dy, dx = sbox.focus()

    if b == c == d == 0.0:
        g = variate(poisson_distribution(mean=a))
        for k in range(dz):
            for j in range(dy):
                for i in range(dx):
                    sbox[k, j, i] += next(g)
    else:
        for k in range(dz):
            for j in range(dy):
                for i in range(dx):
                    pixel = a + b * i + c * j + d * k
                    g = variate(poisson_distribution(mean=pixel))
                    sbox[k, j, i] += next(g)
    return


def simple_gaussian_spots(params):
    from dials.array_family import flex
    from scitbx import matrix

    r = params.rotation
    axis = matrix.col((r.axis.x, r.axis.y, r.axis.z))
    if axis.length() > 0:
        rotation = axis.axis_and_angle_as_r3_rotation_matrix(r.angle, deg=True)
    else:
        rotation = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

    # generate mask and peak values

    from dials.algorithms.shoebox import MaskCode

    mask_peak = MaskCode.Valid | MaskCode.Foreground
    mask_back = MaskCode.Valid | MaskCode.Background

    from dials.util.command_line import ProgressBar

    p = ProgressBar(title="Generating reflections")

    rlist = flex.reflection_table(params.nrefl)
    hkl = flex.miller_index(params.nrefl)
    s1 = flex.vec3_double(params.nrefl)
    xyzmm = flex.vec3_double(params.nrefl)
    xyzpx = flex.vec3_double(params.nrefl)
    panel = flex.size_t(params.nrefl)
    bbox = flex.int6(params.nrefl)

    for j in range(params.nrefl):
        p.update(j * 100.0 / params.nrefl)
        hkl[j] = (random.randint(0, 20), random.randint(0, 20), random.randint(0, 20))
        phi = 2 * math.pi * random.random()
        s1[j] = (0, 0, 0)
        xyzpx[j] = (0, 0, 0)
        xyzmm[j] = (0, 0, phi)
        panel[j] = 0
        bbox[j] = (
            0,
            params.shoebox_size.x,
            0,
            params.shoebox_size.y,
            0,
            params.shoebox_size.z,
        )

    p.finished("Generating %d reflections" % params.nrefl)
    intensity = flex.double(params.nrefl)
    shoebox = flex.shoebox(panel, bbox)
    shoebox.allocate_with_value(MaskCode.Valid)

    p = ProgressBar(title="Generating shoeboxes")

    for i in range(len(rlist)):

        p.update(i * 100.0 / params.nrefl)
        mask = shoebox[i].mask

        if params.pixel_mask == "precise":
            # flag everything as background: peak will me assigned later
            for j in range(len(mask)):
                mask[j] = mask_back
        elif params.pixel_mask == "all":
            # flag we have no idea what anything is
            mask_none = MaskCode.Valid | MaskCode.Foreground | MaskCode.Background
            for j in range(len(mask)):
                mask[j] = mask_none
        elif params.pixel_mask == "static":
            from scitbx.array_family import flex

            x0 = params.spot_offset.x + params.shoebox_size.x / 2
            y0 = params.spot_offset.x + params.shoebox_size.y / 2
            z0 = params.spot_offset.x + params.shoebox_size.z / 2
            sx = params.mask_nsigma * params.spot_size.x
            sy = params.mask_nsigma * params.spot_size.y
            sz = params.mask_nsigma * params.spot_size.z

            # The x, y, z indices
            z, y, x = zip(*itertools.product(*(range(n) for n in mask.all())))
            xyz = flex.vec3_double(flex.double(x), flex.double(y), flex.double(z))

            # Calculate SUM(((xj - xj0) / sxj)**2) for each element
            xyz0 = (x0, y0, z0)
            isxyz = (1.0 / sx, 1.0 / sy, 1.0 / sz)
            dxyz = sum(
                [
                    (x * isx) ** 2
                    for x, isx in zip(((xyz - xyz0) * rotation).parts(), isxyz)
                ]
            )

            # Set the mask values
            index = dxyz <= 1.0
            index.reshape(mask.accessor())
            mask.set_selected(index, MaskCode.Valid | MaskCode.Foreground)
            mask.set_selected(not index, MaskCode.Valid | MaskCode.Background)

        sbox = shoebox[i].data

        # reflection itself, including setting the peak region if we're doing that
        # FIXME use flex arrays to make the rotation bit more efficient as this is
        # now rather slow...

        counts_true = 0
        for j in range(params.counts):
            _x = random.gauss(0, params.spot_size.x)
            _y = random.gauss(0, params.spot_size.y)
            _z = random.gauss(0, params.spot_size.z)

            Rxyz = rotation * matrix.col((_x, _y, _z)).elems

            x = int(Rxyz[0] + params.spot_offset.x + params.shoebox_size.x / 2)
            y = int(Rxyz[1] + params.spot_offset.y + params.shoebox_size.y / 2)
            z = int(Rxyz[2] + params.spot_offset.z + params.shoebox_size.z / 2)

            if x < 0 or x >= params.shoebox_size.x:
                continue
            if y < 0 or y >= params.shoebox_size.y:
                continue
            if z < 0 or z >= params.shoebox_size.z:
                continue
            sbox[z, y, x] += 1
            counts_true += 1
            if params.pixel_mask == "precise":
                mask[z, y, x] = mask_peak

        intensity[i] = counts_true

        if params.background:
            # background:flat;
            for j in range(params.background * len(sbox)):
                x = random.randint(0, params.shoebox_size.x - 1)
                y = random.randint(0, params.shoebox_size.y - 1)
                z = random.randint(0, params.shoebox_size.z - 1)
                sbox[z, y, x] += 1
        else:
            # or inclined
            random_background_plane(
                sbox,
                params.background_a,
                params.background_b,
                params.background_c,
                params.background_d,
            )

    rlist["miller_index"] = hkl
    rlist["s1"] = s1
    rlist["xyzcal.px"] = xyzpx
    rlist["xyzcal.mm"] = xyzmm
    rlist["bbox"] = bbox
    rlist["panel"] = panel
    rlist["shoebox"] = shoebox
    rlist["intensity.sum.value"] = intensity
    p.finished("Generating %d shoeboxes" % params.nrefl)

    return rlist


# FIXME to do:
#
#  - generate the reflections in reciprocal space w.r.t. a reciprocal space
#    node; will require taking a position then calculating a prediction for it
#    (i.e. when it will be in reflecting position) then extrapolate this and
#    add it to the appropriate voxel
#  - implement this as a flex array of events then perform this mapping in C++
#    to the target space, then accumulate into the shoebox in the target space
#    (can this be done in C++ as well?)
#  - apply the detector geometric correction due to the depth of the pixel, will
#    depend on the actual detector geometry and so on
#  - include input of beam, crystal, detector, scan, goniometer models so that
#    the predictions are realistic
#  - worth noting that if I can make this work it could form the kernel of old
#    deadly - taking all of the dials models as input along with some profile
#    parameters and generating an image of a spot...


def background_xds(rlist):
    # FIXME
    # from dials.extensions import XdsBackgroundExt
    # background = XdsBackgroundExt(None, None)
    # background.compute_background(rlist)

    # from math import erf, sqrt
    # from dials.algorithms.shoebox import MaskCode
    # from dials.array_family import flex
    # nsigma = 3
    # shoebox = rlist['shoebox']
    # for i in range(len(shoebox)):
    # data = shoebox[i].data
    # mask = shoebox[i].mask
    # bkgrd = shoebox[i].background

    # fg = flex.bool([m & MaskCode.Foreground != 0 for m in mask])
    # bg = flex.bool([m & MaskCode.BackgroundUsed != 0 for m in mask])
    # n = fg.count(True)
    # m = bg.count(True)
    # Sp = flex.sum(data.select(fg))
    # Bp = flex.sum(data.select(bg)) / m

    # F = 1.0 - erf(nsigma / sqrt(2.0))
    # c = 1.0 / (1 + F * n / m)
    # B = c * (Bp - (F / m) * Sp)
    ##   print F, c, Bp, n, m, Sp, B
    # print Bp, Sp
    # for j in range(len(bkgrd)):
    # bkgrd[j] = B

    return


def background_inclined(rlist):
    from dials.algorithms.background import InclinedSubtractor

    background = InclinedSubtractor()
    background(None, rlist)


def integrate_3d_summation(rlist):
    from dials.algorithms.integration.sum import IntegrationAlgorithm

    integration = IntegrationAlgorithm()
    integration(rlist)


def main(params):
    # first generate the reflections - this could be called from elsewhere
    rlist = simple_gaussian_spots(params)
    correct_intensities = list(rlist["intensity.sum.value"])
    del rlist["intensity.sum.value"]

    # now integrate those reflections using code from elsewhere
    if params.background_method == "xds":
        background_xds(rlist)
    elif params.background_method == "mosflm":
        assert params.pixel_mask != "all"
        background_inclined(rlist)

    integrate_3d_summation(rlist)

    integrated_intensities = list(rlist["intensity.sum.value"])

    # now scan through the reflection list and find those where the integration
    # gave an apparently duff answer i.e. outside of 3 sigma from correct value

    from dials.array_family import flex

    overestimates = []
    underestimates = []

    for j, (c, i) in enumerate(zip(correct_intensities, integrated_intensities)):
        sigma = math.sqrt(c)
        if math.fabs(c - i) < 3 * sigma:
            continue
        if i > params.counts:
            overestimates.append(j)
        else:
            underestimates.append(j)

    print(
        "%d overestimates, %d underestimates"
        % (len(overestimates), len(underestimates))
    )

    overestimates = rlist.select(flex.size_t(overestimates))
    underestimates = rlist.select(flex.size_t(underestimates))
    # now pickle these, perhaps

    import six.moves.cPickle as pickle

    if params.output.under:
        with open(params.output.under, "wb") as fh:
            pickle.dump(underestimates, fh, pickle.HIGHEST_PROTOCOL)
    if params.output.over:
        with open(params.output.over, "wb") as fh:
            pickle.dump(overestimates, fh, pickle.HIGHEST_PROTOCOL)
    if params.output.all:
        with open(params.output.all, "wb") as fh:
            pickle.dump(rlist, fh, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    import sys

    from libtbx.phil import command_line

    cmd = command_line.argument_interpreter(master_params=master_phil)
    working_phil = cmd.process_and_fetch(args=sys.argv[1:])
    main(working_phil.extract())
