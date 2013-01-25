from __future__ import division

import sys
import math
import os
import itertools
import x2tbx

from cctbx.array_family import flex
from iotbx import mtz
from cctbx.miller import build_set
from cctbx.crystal import symmetry as crystal_symmetry
from scitbx import lbfgs
from libtbx.phil import parse

def nint(a):
    return int(round(a))

def poly_residual(xp, y, params):
    '''Compute the residual between the observations y[i] and sum_j
    params[j] x[i]^j. For efficiency, x[i]^j are pre-calculated in xp.'''

    r = 0.0

    n = len(params)
    c = len(y)

    e = flex.double([flex.sum(xp[j] * params) for j in range(c)])

    return flex.sum(flex.pow2(y - e))

def poly_gradients(xp, y, params):
    '''Compute the gradient of the residual w.r.t. the parameters, N.B.
    will be performed using a finite difference method. N.B. this should
    be trivial to do algebraicly.'''

    eps = 1.0e-6

    g = flex.double()

    n = len(params)

    for j in range(n):
        rs = []
        for signed_eps in [- eps, eps]:
            params_eps = params[:]
            params_eps[j] += signed_eps
            rs.append(poly_residual(xp, y, params_eps))
        g.append((rs[1] - rs[0]) / (2 * eps))

    return g

class poly_fitter:
    '''A class to do the polynomial fit. This will fit observations y
    at points x with a polynomial of order n.'''

    def __init__(self, points, values, order):
        self.x = flex.double([1.0 for j in range(order)])
        self._x = flex.double(points)
        self._y = flex.double(values)

        # precalculate x[j]^[0-(n - 1)] values

        self._xp = [flex.double([math.pow(x, j) for j in range(order)])
                    for x in self._x]

        return

    def refine(self):
        '''Actually perform the parameter refinement.'''

        return lbfgs.run(target_evaluator = self)

    def compute_functional_and_gradients(self):

        return poly_residual(self._xp, self._y, self.x), \
               poly_gradients(self._xp, self._y, self.x)

    def get_parameters(self):
        return list(self.x)

    def evaluate(self, x):
        '''Evaluate the resulting fit at point x.'''

        return sum([math.pow(x, k) * self.x[k] for k in range(len(self.x))])

def fit(x, y, order):
    '''Fit the values y(x) then return this fit. x, y should
    be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma.'''

    pf = poly_fitter(x, y, order)
    pf.refine()

    return [pf.evaluate(_x) for _x in x]

def log_fit(x, y, order):
    '''Fit the values log(y(x)) then return exp() to this fit. x, y should
    be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma.'''

    ly = [math.log(_y) for _y in y]

    pf = poly_fitter(x, ly, order)
    pf.refine()

    return [math.exp(pf.evaluate(_x)) for _x in x]

def log_inv_fit(x, y, order):
    '''Fit the values log(1 / y(x)) then return the inverse of this fit.
    x, y should be iterables, the order of the polynomial for the transformed
    fit needs to be specified. This will be useful for e.g. Rmerge.'''

    ly = [math.log(1.0 / _y) for _y in y]

    pf = poly_fitter(x, ly, order)
    pf.refine()

    return [(1.0 / math.exp(pf.evaluate(_x))) for _x in x]

def interpolate_value(x, y, t):
    '''Find the value of x: y(x) = t.'''

    if t > max(y) or t < min(y):
        raise RuntimeError, 't outside of [%f, %f]' % (min(y), max(y))

    for j in range(1, len(x)):
        x0 = x[j - 1]
        y0 = y[j - 1]

        x1 = x[j]
        y1 = y[j]

        if (y0 - t) * (y1 - t) < 0:
            return x0 + (t - y0) * (x1 - x0) / (y1 - y0)

def get_positive_values(x):
    '''Return a list of values v from x where v > 0.'''

    result = []

    for _x in x:
        if _x > 0:
            result.append(_x)
        else:
            return result

    return result

phil_defaults = '''
resolutionizer {
  rmerge = 0.0
    .type = float
  completeness = 0.0
    .type = float
  isigma = 1.0
    .type = float
  misigma = 2.0
    .type = float
  nbins = 100
    .type = int
}
'''

class resolutionizer:
    '''A class to calculate things from merging reflections.'''

    def __init__(self, args):

        working_phil = parse(phil_defaults)
        phil_args = []

        for arg in args[1:]:
            if os.path.exists(arg):
                working_phil = working_phil.fetch(
                    source = parse(open(arg).read()))
            else:
                phil_args.append(arg)

        for phil_arg in phil_args:
            interp = working_phil.command_line_argument_interpreter(
                home_scope = 'resolutionizer')
            more_phil = interp.process(phil_arg)
            working_phil = working_phil.fetch(source = more_phil)

        self._params = working_phil.extract().resolutionizer

        # look for useful information in the input file
        mtz_file = args[0]
        mtz_obj = mtz.object(mtz_file)

        i_data = None
        sigi_data = None

        mi = mtz_obj.extract_miller_indices()

        unit_cell = None

        for crystal in mtz_obj.crystals():
            unit_cell = crystal.unit_cell()
            for dataset in crystal.datasets():
                for column in dataset.columns():
                    if column.label() == 'I':
                        i_data = column.extract_values(
                            not_a_number_substitute = 0.0)
                    if column.label() == 'SIGI':
                        sigi_data = column.extract_values(
                            not_a_number_substitute = 0.0)

        assert(i_data)
        assert(sigi_data)

        self._unit_cell = unit_cell
        self._space_group = mtz_obj.space_group()

        self._r = x2tbx.ReflectionList()
        self._r.setup(mi, i_data, sigi_data)
        self._r.set_unit_cell(unit_cell.parameters())
        self._r.merge()

        return

    def calculate_resolution_ranges(self, nbins = 20):
        '''Calculate semi-useful resolution ranges for analysis.'''

        self._r.setup_resolution_shells(nbins)

        high = self._r.shell_high_limits()
        low = self._r.shell_low_limits()

        self._resolution_ranges = [(h, l) for h, l in zip(high, low)]
        self._hkl_ranges = [self._r.get_shell(j) for j in range(nbins)]

        return

    def get_resolution_bins(self):
        '''Return the reversed resolution limits - N.B. this is most
        important when considering resolution calculations, see
        resolution_completeness.'''

        return list(reversed(self._hkl_ranges)), \
               list(reversed(self._resolution_ranges))

    def calculate_completeness(self, resolution_bin):
        '''Calculate the completeness of observations in a named
        resolution bin.'''

        resolution_range = self._resolution_ranges[resolution_bin]
        hkl_list = self._hkl_ranges[resolution_bin]

        uc = self._unit_cell
        sg = self._space_group

        dmin = min(resolution_range)
        dmax = max(resolution_range)

        cs = crystal_symmetry(unit_cell = uc, space_group = sg)
        hkl_calc = [hkl for hkl in build_set(
            cs, False, d_min = dmin, d_max = dmax).indices()]

        # remove systematically absent reflections

        hkl_list = [hkl for hkl in
                    itertools.ifilterfalse(sg.is_sys_absent, hkl_list)]

        return float(len(hkl_list)) / float(len(hkl_calc))

    def calculate_rmerge(self, bin):
        '''Calculate the overall Rmerge for a given bin.'''

        # FIXME this is not efficient - should calculate all and cache...

        return self._r.rmerge_shells()[bin]

    def calculate_merged_isigma(self, bin):
        '''Calculate the average merged I/sigma.'''

        # FIXME this is not efficient - should calculate all and cache...

        return self._r.i_sigma_shells()[bin]

    def calculate_unmerged_isigma(self, bin):
        '''Calculate the average unmerged I/sigma.'''

        # FIXME this is not efficient - should calculate all and cache...

        return self._r.total_i_sigma_shells()[bin]

    def resolution_auto(self):
        '''Compute resolution limits based on the current self._params set.'''

        self.calculate_resolution_ranges(nbins = self._params.nbins)

        if self._params.rmerge:
            print 'Resolution rmerge:       %.2f' % \
                self.resolution_rmerge()

        if self._params.completeness:
            print 'Resolution completeness: %.2f' % \
                self.resolution_completeness()

        if self._params.isigma:
            print 'Resolution I/sig:        %.2f' % \
                self.resolution_unmerged_isigma()

        if self._params.misigma:
            print 'Resolution Mn(I/sig):    %.2f' % \
                self.resolution_merged_isigma()

        return

    def resolution_rmerge(self, limit = None):
        '''Compute a resolution limit where either rmerge = 1.0 (limit if
        set) or the full extent of the data. N.B. this fit is only meaningful
        for positive values.'''

        if limit is None:
            limit = self._params.rmerge

        bins, ranges = self.get_resolution_bins()

        if limit == 0.0:
            return ranges[-1][0]

        rmerge_s = get_positive_values(
            [self.calculate_rmerge(bin)
             for bin in reversed(range(len(bins)))])

        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(rmerge_s)]

        if limit == 0.0:
            return 1.0 / math.sqrt(max(s_s))

        if limit > max(rmerge_s):
            return 1.0 / math.sqrt(max(s_s))

        rmerge_f = log_inv_fit(s_s, rmerge_s, 6)

        try:
            r_rmerge = 1.0 / math.sqrt(interpolate_value(s_s, rmerge_f, limit))
        except: # intentional
            r_rmerge = 1.0 / math.sqrt(max(s_s))

        return r_rmerge

    def new_resolution_unmerged_isigma(self, limit = None):
        '''Compute a resolution limit where either I/sigma = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.isigma

        bins, ranges = self.get_resolution_bins()

        isigma_s = get_positive_values(
            [self.calculate_unmerged_isigma(bin)
             for bin in reversed(range(len(bins)))])

        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(isigma_s)]

        if min(isigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        for _l, s in enumerate(isigma_s):
            if s < limit:
                break

        if _l > 10 and _l < (len(isigma_s) - 10):
            start = _l - 10
            end = _l + 10
        elif _l <= 10:
            start = 0
            end = 20
        elif _l >= (len(isigma_s) - 10):
            start = -20
            end = -1

        _s_s = s_s[start:end]
        _isigma_s = isigma_s[start:end]

        _isigma_f = log_fit(_s_s, _isigma_s, 3)

        try:
            r_isigma = 1.0 / math.sqrt(interpolate_value(_s_s, _isigma_f,
                                                         limit))
        except: # intentional
            r_isigma = 1.0 / math.sqrt(max(_s_s))

        return r_isigma

    def resolution_unmerged_isigma(self, limit = None):
        '''Compute a resolution limit where either I/sigma = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.isigma

        bins, ranges = self.get_resolution_bins()

        isigma_s = get_positive_values(
            [self.calculate_unmerged_isigma(bin)
             for bin in reversed(range(len(bins)))])

        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(isigma_s)]

        if min(isigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        isigma_f = log_fit(s_s, isigma_s, 6)

        try:
            r_isigma = 1.0 / math.sqrt(interpolate_value(s_s, isigma_f, limit))
        except: # intentional
            r_isigma = 1.0 / math.sqrt(max(s_s))

        return r_isigma

    def new_resolution_merged_isigma(self, limit = None):
        '''Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.misigma

        bins, ranges = self.get_resolution_bins()

        misigma_s = get_positive_values(
            [self.calculate_merged_isigma(bin) for bin in range(len(bins))])
        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(misigma_s)]

        if min(misigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        for _l, s in enumerate(misigma_s):
            if s < limit:
                break

        if _l > 10 and _l < (len(misigma_s) - 10):
            start = _l - 10
            end = _l + 10
        elif _l <= 10:
            start = 0
            end = 20
        elif _l >= (len(misigma_s) - 10):
            start = -20
            end = -1

        _s_s = s_s[start:end]
        _misigma_s = misigma_s[start:end]

        _misigma_f = log_fit(_s_s, _misigma_s, 3)

        try:
            r_misigma = 1.0 / math.sqrt(interpolate_value(_s_s, _misigma_f,
                                                          limit))
        except: # intentional
            r_misigma = 1.0 / math.sqrt(max(_s_s))

        return r_misigma

    def resolution_merged_isigma(self, limit = None):
        '''Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.misigma

        bins, ranges = self.get_resolution_bins()

        misigma_s = get_positive_values(
            [self.calculate_merged_isigma(bin)
             for bin in reversed(range(len(bins)))])
        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(misigma_s)]

        if min(misigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        misigma_f = log_fit(s_s, misigma_s, 6)

        try:
            r_misigma = 1.0 / math.sqrt(
                interpolate_value(s_s, misigma_f, limit))
        except: # intentional
            r_misigma = 1.0 / math.sqrt(max(s_s))

        return r_misigma

    def resolution_completeness(self, limit = None):
        '''Compute a resolution limit where completeness < 0.5 (limit if
        set) or the full extent of the data. N.B. this completeness is
        with respect to the *maximum* completeness in a shell, to reflect
        triclinic cases.'''

        if limit is None:
            limit = self._params.completeness

        bins, ranges = self.get_resolution_bins()

        s_s = [1.0 / (r[0] * r[0]) for r in reversed(ranges)]

        if limit == 0.0:
            return 1.0 / math.sqrt(max(s_s))

        comp_s = [self.calculate_completeness(j) for j, bin in enumerate(
            reversed(bins))]

        if min(comp_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        comp_f = fit(s_s, comp_s, 6)

        rlimit = limit * max(comp_s)

        try:
            r_comp = 1.0 / math.sqrt(
                interpolate_value(s_s, comp_f, rlimit))
        except: # intentional
            r_comp = 1.0 / math.sqrt(max(s_s))

        return r_comp

if __name__ == '__main__':

    m = resolutionizer(sys.argv[1:])
    m.resolution_auto()
