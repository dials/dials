def pull_reference(integrate_hkl):
    '''Generate reference data set from integrate.hkl - check out the calculated 
    x, y and z centroids as well as the Miller indices as coordinates in some
    high dimensional space. Only consider measurements with meaningful 
    centroids...'''

    hkl = []
    i = []
    sigi = []
    xyz = []

    for record in open(integrate_hkl):
        if record.startswith('!'):
            continue

        f_tokens = map(float, record.split())

        if f_tokens[12:15] == [0.0, 0.0, 0.0]:
            continue

        hkl.append(tuple(map(int, f_tokens[0:3])))
        i.append(f_tokens[3])
        sigi.append(f_tokens[4])
        xyz.append(tuple(f_tokens[5:8]))

    print 'Reference: %d observations' % len(hkl)
    return hkl, i, sigi, xyz

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

def compare(integrate_hkl, integrate_pkl):

    from cctbx.array_family import flex
    from annlib_ext import AnnAdaptor as ann_adaptor

    xhkl, xi, xsigi, xxyz = pull_reference(integrate_hkl)
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

    xds = []
    dials = []
    
    # perform the analysis
    for j, hkl in enumerate(dhkl):
        c = ann.nn[j]
        if hkl == xhkl[c]:
            xds.append(xi[c])
            dials.append(di[j])

    print 'Paired %d observations' % len(xds)
            
    c = cc(xds, dials)

    print 'Correlation coefficient: %.2f' % c

if __name__ == '__main__':
    import sys
    compare(sys.argv[1], sys.argv[2])
