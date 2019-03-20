from __future__ import absolute_import, division, print_function

import libtbx.phil
from iotbx import mtz
from scitbx.array_family import flex

# LIBTBX_SET_DISPATCHER_NAME dev.dials.deviant

scope = libtbx.phil.parse(
    """
  hklin = None
    .type = path
  column = 'I'
    .type = str
"""
)


def main(args):
    work = scope
    for arg in args:
        work = work.fetch(libtbx.phil.parse(arg))

    p = work.extract()
    m = mtz.object(p.hklin)
    mad = m.as_miller_arrays_dict(merge_equivalents=False)

    z = flex.double()

    Icol = [mad[k] for k in mad if k[-1] == p.column][0]
    hkl = Icol.indices()

    mIcol = Icol.merge_equivalents().array()
    mhkl = mIcol.indices()

    for j, _mhkl in enumerate(mhkl):
        obs = Icol.data().select(hkl == _mhkl)
        ref = mIcol.data()[j]
        for o in obs:
            z.append((o - ref) / ref)
            # shopping
            # print((o - ref) / ref, obs.size())
            # if z > 5:
            #    print('%d %d %d' % _mhkl)

    zhist = flex.histogram(z, data_min=-10, data_max=10, n_slots=20)

    for c, v in zip(zhist.slot_centers(), zhist.slots()):
        print("%.2f %d" % (c, v))


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
