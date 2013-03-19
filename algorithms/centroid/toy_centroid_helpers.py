
def covar_contrib(x1, x2, weights):
    '''Compute contribution to covariance matrix.'''

    covar = sum([_x1 * _x2 * _w for _x1, _x2, _w in zip(x1, x2, weights)]) / \
        sum(weights)

    return covar

def form_covariance_matrix(pixel_list, _f, _r, _c):
    '''Compute covariance matrix.'''

    # form the variance matrix

    pixels = [(f - _f, r - _r, c - _c, d) for f, r, c, d in pixel_list]

    data = ([pixel[0] for pixel in pixels],
            [pixel[1] for pixel in pixels],
            [pixel[2] for pixel in pixels],
            [pixel[3] for pixel in pixels])

    m_elems = []

    for i in range(3):
        for j in range(3):
            m_elems.append(covar_contrib(data[i], data[j], data[3]))

    from scitbx.array_family import flex
    from scitbx.linalg import eigensystem
    from scitbx import matrix

    m = flex.double(flex.grid(3, 3))
    for j in range(9):
        m[j] = m_elems[j]
    s = eigensystem.real_symmetric(m)
    values = s.values()
    vectors = s.vectors()

    for j in range(3):
        vec = tuple(vectors[j * 3:j * 3 + 3])
        print '%.3f %.3f %.3f %.3f' % (values[j], vec[0], vec[1], vec[2])

    return
