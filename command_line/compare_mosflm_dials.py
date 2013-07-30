def pull_reference(integrate_mtz):
    '''Generate reference data set from integrate.hkl, check out the calculated
    x, y and z centroids as well as the Miller indices as coordinates in some
    high dimensional space. Only consider measurements with meaningful
    centroids...'''

    # prepare input data as
    # sortmtz hklin thau_2_001.mtz hklout sorted.mtz << eof
    # H K L M/ISYM BATCH
    # eof
    #
    # scala hklin sorted.mtz hklout summed.mtz << eof
    # run 1 all
    # scales constant
    # output unmerged
    # sdcorrection noadjust norefine both 1 0 0
    # eof

    from iotbx import mtz

    m = mtz.object(integrate_mtz)

    hkl = m.extract_original_index_miller_indices()

    b0 = m.batches()[0]

    for c in m.columns():
        if c.label() == 'I':
            i = c.extract_valid_values()
        elif c.label() == 'SIGI':
            sigi = c.extract_valid_values()
        elif c.label() == 'XDET':
            xdet = c.extract_valid_values()
        elif c.label() == 'YDET':
            ydet = c.extract_valid_values()
        elif c.label() == 'ROT':
            rot = c.extract_valid_values()

    # extract x, y, z positions for this in image address space

    xyz = []

    dx = b0.detlm()[1]
    dz = b0.phiend() - b0.phistt()
    z0 = b0.phistt() / dz

    for x, y, r in zip(xdet, ydet, rot):
        _x = y
        _y = dx - x
        _z = z0 + r / dz

        xyz.append((_x, _y, _z))

    print 'Reference: %d observations' % len(hkl)
    return hkl, i, sigi, xyz

def integrate_mtz_to_unit_cell(integrate_mtz):
    '''Generate a cctbx unit_cell from an integrate_mtz file.'''

    from iotbx import mtz
    m = mtz.object(integrate_mtz)
    for c in m.crystals():
        return c.unit_cell()

    raise RuntimeError, 'unit cell not found'

def pull_calculated(integrate_pkl):
    from dials.model.data import ReflectionList
    import cPickle as pickle
    import math

    r_list = pickle.load(open(integrate_pkl, 'rb'))

    strong_reflections = []

    for r in r_list:
        if r.intensity > math.sqrt(r.intensity_variance):
            strong_reflections.append(r)

    del(r_list)

    hkl = []
    i = []
    sigi = []
    xyz = []

    for r in strong_reflections:
        if not r.is_valid():
            continue
        hkl.append(r.miller_index)
        i.append(r.intensity)
        sigi.append(math.sqrt(r.intensity_variance))
        x, y = r.image_coord_px
        z = r.frame_number
        xyz.append((x, y, z))

    print 'Computed: %d observations' % len(hkl)
    return hkl, i, sigi, xyz

def meansd(values):
    import math

    assert(len(values) > 3)

    mean = sum(values) / len(values)
    var = sum([(v - mean) * (v - mean) for v in values]) / (len(values) - 1)

    return mean, math.sqrt(var)

def cc(a, b):

    assert(len(a) == len(b))

    ma, sa = meansd(a)
    mb, sb = meansd(b)

    r = (1 / (len(a) - 1)) * sum([((a[j] - ma) / sa) * ((b[j] - mb) / sb)
                                  for j in range(len(a))])

    return r

def R(calc, obs, scale = None):

    import math

    assert(len(calc) == len(obs))

    if not scale:
        scale = sum(obs) / sum(calc)

    return sum([math.fabs(math.fabs(o) - math.fabs(scale * c)) \
                for c, o in zip(calc, obs)]) / \
                sum([math.fabs(o) for o in obs]), scale

def compare_chunks(integrate_mtz, integrate_pkl, crystal_json):

    from cctbx.array_family import flex
    from annlib_ext import AnnAdaptor as ann_adaptor

    # FIXME read crystal.json file

    uc = integrate_mtz_to_unit_cell(integrate_mtz)

    # FIXME read UB matrix from integrate.mtz


    # FIXME compute mapping from one to the other; compute reindex matrix


    xhkl, xi, xsigi, xxyz = pull_reference(integrate_mtz)
    dhkl, di, dsigi, dxyz = pull_calculated(integrate_pkl)

    reference = flex.double()
    query = flex.double()

    for xyz in xxyz:
        reference.append(xyz[0])
        reference.append(xyz[1])
        reference.append(xyz[2])

    for xyz in dxyz:
        query.append(xyz[0])
        query.append(xyz[1])
        query.append(xyz[2])

    # perform the match
    ann = ann_adaptor(data = reference, dim = 3, k = 1)
    ann.query(query)

    XDS = []
    DIALS = []
    HKL = []

    fout = open('ratio_xyz.dat', 'w')

    # perform the analysis
    for j, hkl in enumerate(dhkl):
        c = ann.nn[j]
        if hkl == xhkl[c]:
            XDS.append(xi[c])
            DIALS.append(di[j])
            HKL.append(hkl)
            fout.write('%d %d %d %f %f %f %f %f %f %f %f %f\n' % \
                       (hkl[0], hkl[1], hkl[2], uc.d(hkl),
                        xi[c], xxyz[c][0], xxyz[c][1], xxyz[c][2],
                        di[j], dxyz[j][0], dxyz[j][1], dxyz[j][2]))

    fout.close()

    # now compute resolution for every reflection - or at least each unique
    # Miller index...

    unique = set(HKL)

    resolutions = { }

    for hkl in unique:
        resolutions[hkl] = uc.d(hkl)

    # then resort the list in terms of resolution, then reverse it

    sort_me = []
    for hkl, xds, dials in zip(HKL, XDS, DIALS):
        sort_me.append((resolutions[hkl], xds, dials))

    sort_me.sort()
    sort_me.reverse()

    resolutions = [sm[0] for sm in sort_me]
    XDS = [sm[1] for sm in sort_me]
    DIALS = [sm[2] for sm in sort_me]

    # then extract the original observation structure

    print 'Paired %d observations' % len(XDS)

    scale = sum(XDS) / sum(DIALS)

    chunks = [(i, i + 1000) for i in range(0, len(XDS), 1000)]

    ccs = []
    rs = []
    ss = []

    for chunk in chunks:
        xds = XDS[chunk[0]:chunk[1]]
        dials = DIALS[chunk[0]:chunk[1]]
        resols = resolutions[chunk[0]:chunk[1]]

        if len(xds) < 100:
            break

        c = cc(dials, xds)
        r, s = R(dials, xds)
        print '%7d %4d %.3f %.3f %.3f %.3f %.3f' % (chunk[0], len(xds),
                                                    min(resols), max(resols),
                                                    c, r, s)
        ccs.append(c)
        rs.append(r)
        ss.append(s)

    chunks = [j for j in range(len(chunks))]

    # kludge - if we fall off

    chunks = chunks[:len(rs)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, ccs, label = 'CC')
    pyplot.plot(chunks, rs, label = 'R')
    pyplot.plot(chunks, ss, label = 'K')
    pyplot.legend()
    pyplot.savefig('plot.png')
    pyplot.close()

    return

##### FIXME derive REIDX matrix if appropriate from the unit cell vectors #####

if __name__ == '__main__':
    import sys
    compare_chunks(sys.argv[1], sys.argv[2])
