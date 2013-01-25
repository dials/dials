#!/usr/bin/env cctbx.python
# Resolutionizer.py
#
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A toolkit component for merging symmetry related intensity measurements
# to calculate:
#
#  - Rmerge vs. batch, resolution
#  - Chi^2
#  - Multiplicity
#  - Unmerged I/sigma
#  - Z^2 for centric and acentric reflections
#  - Completeness
#
# The standalone (moving to C++) version...

import sys
import math
import os
import time
import itertools
import copy

from cctbx.array_family import flex
from iotbx import mtz
from cctbx.miller import build_set
from cctbx.miller import map_to_asu
from cctbx.crystal import symmetry as crystal_symmetry
from cctbx.sgtbx import rt_mx
from scitbx import lbfgs
from libtbx.phil import parse

def nint(a):
    return int(round(a))

class mtz_dataset:
    '''A class to represent the MTZ dataset in the hierarchy. This will
    be instantiated in the mtz_crystal class below, and contain:

     - a list of columns

    Maybe more things will be added.'''

    def __init__(self, iotbx_dataset):
        self._name = iotbx_dataset.name()
        self._column_table = { }
        for column in iotbx_dataset.columns():
            self._column_table[column.label()] = column

        return

    def get_column_names(self):
        return list(self._column_table)

    def get_column(self, column_label):
        return self._column_table[column_label]

    def get_column_values(self, column_label, nan_value = 0.0):
        return self._column_table[column_label].extract_values(
            not_a_number_substitute = nan_value)

class mtz_crystal:
    '''A class to represent the MTZ crystal in the hierarchy. This will
    be instantiated by the factories below.'''

    def __init__(self, iotbx_crystal):
        self._name = iotbx_crystal.name()
        self._unit_cell = iotbx_crystal.unit_cell()

        self._dataset_table = { }
        for dataset in iotbx_crystal.datasets():
            self._dataset_table[dataset.name()] = mtz_dataset(dataset)

        self._column_table = { }

        for dname in self._dataset_table:
            dataset = self._dataset_table[dname]
            for column_name in dataset.get_column_names():
                assert(not column_name in self._column_table)
                self._column_table[column_name] = dataset.get_column(
                    column_name)

        return

    def get_dataset_names(self):
        return list(self._dataset_table)

    def get_dataset(self, dataset_name):
        return self._dataset_table[dataset_name]

    def get_unit_cell_parameters(self):
        return tuple(self._unit_cell.parameters())

    def get_unit_cell(self):
        return self._unit_cell

    def get_column_names(self):
        return list(self._column_table)

    def get_column(self, column_label):
        return self._column_table[column_label]

    def get_column_values(self, column_label, nan_value = 0.0):
        return self._column_table[column_label].extract_values(
            not_a_number_substitute = nan_value)

class mtz_file:
    '''A class to represent the full MTZ file in the hierarchy - this
    will have a list of one or more crystals contained within it each
    with its own unit cell and datasets.'''

    # FIXME need to keep in mind MTZ batch headers - can I access these?
    # yes - through stuff like this:
    #
    #    for batch in mtz_obj.batches():
    #        for token in dir(batch):
    #            print token
    #        print batch.num()
    #
    # but need to decide what I am looking for... will bake this in with
    # a future update (it would be interesting to look at how the UB
    # matrices behave.)

    def __init__(self, hklin):
        mtz_obj = mtz.object(hklin)

        self._miller_indices = mtz_obj.extract_miller_indices()
        self._resolution_range = mtz_obj.max_min_resolution()
        self._space_group = mtz_obj.space_group()

        self._crystal_table = { }

        for crystal in mtz_obj.crystals():
            self._crystal_table[crystal.name()] = mtz_crystal(crystal)

        self._column_table = { }

        for xname in self._crystal_table:
            crystal = self._crystal_table[xname]
            for dname in crystal.get_dataset_names():
                dataset = crystal.get_dataset(dname)
                for column_name in dataset.get_column_names():
                    assert(not column_name in self._column_table)
                    self._column_table[column_name] = dataset.get_column(
                        column_name)

        return

    def get_crystal_names(self):
        return list(self._crystal_table)

    def get_crystal(self, crystal_name):
        return self._crystal_table[crystal_name]

    def get_unit_cell(self):
        '''Get the unit cell object from HKL_base for other calculations.'''

        return self.get_crystal('HKL_base').get_unit_cell()

    def get_space_group(self):
        return self._space_group

    def get_resolution_range(self):
        return self._resolution_range

    def get_symmetry_operations(self):
        return [smx for smx in self._space_group.smx()]

    def get_centring_operations(self):
        return [ltr for ltr in self._space_group.ltr()]

    def get_miller_indices(self):
        return self._miller_indices

    def get_column_names(self):
        return list(self._column_table)

    def get_column(self, column_label):
        return self._column_table[column_label]

    def get_column_values(self, column_label, nan_value = 0.0):
        return self._column_table[column_label].extract_values(
            not_a_number_substitute = nan_value)

def mtz_dump(hklin):
    '''An implementation of mtzdump using the above classes.'''

    mtz = mtz_file(hklin)

    print 'Reading file: %s' % hklin
    print 'Spacegroup: %s' % mtz.get_space_group().type(
         ).universal_hermann_mauguin_symbol()

    print 'Centring operations:'
    for cenop in mtz.get_centring_operations():
        print cenop

    print 'Symmetry operations:'
    for symop in mtz.get_symmetry_operations():
        print symop

    for xname in mtz.get_crystal_names():
        crystal = mtz.get_crystal(xname)
        print 'Crystal: %s' % xname
        print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f' % \
              crystal.get_unit_cell_parameters()

        for dname in crystal.get_dataset_names():
            dataset = crystal.get_dataset(dname)
            print 'Dataset: %s' % dname
            print 'Columns (with min / max)'
            for column in dataset.get_column_names():
                values = dataset.get_column_values(column)
                print '%20s %.4e %.4e' % (column, min(values), max(values))

    print 'All columns:'
    for column in mtz.get_column_names():
        print column

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

class unmerged_intensity:
    '''A class to represent and encapsulate the multiple observations of a
    given intensity defined in terms of the Miller index. It is assumed that
    these are compatible observations.'''

    def __init__(self):
        self._observations = []

        return

    def add(self, misym, i, sigi, b):
        self._observations.append((misym, i, sigi, b))
        return

    def apply_scale(self, s):
        '''Apply scale factor s.'''

        scaled_observations = []

        for o in self._observations:
            scaled_observations.append((o[0], s * o[1], s * o[2], o[3]))

        self._observations = scaled_observations

        return

    def get(self):
        return self._observations

    def merge(self):
        '''Merge the observations recorded so far, return Imean, SigImean.'''

        assert(self._observations)

        sum_wi = 0.0
        sum_w = 0.0

        for o in self._observations:
            i = o[1]
            w = 1.0 / (o[2] * o[2])
            sum_w += w
            sum_wi += w * i

        i_mean = sum_wi / sum_w
        sigi_mean = math.sqrt(1.0 / sum_w)

        return i_mean, sigi_mean

    def merge_anomalous(self):
        '''Merge the observations recorded so far, return Imean+, SigImean+
        &c. N.B. will return a 4-ple. +/- defined by M_ISYM record, which
        if odd = I+, even I-.'''

        assert(self._observations)

        sum_wi_p = 0.0
        sum_w_p = 0.0

        sum_wi_m = 0.0
        sum_w_m = 0.0

        for o in self._observations:
            i = o[1]
            w = 1.0 / (o[2] * o[2])
            if o[0] % 2:
                sum_w_p += w
                sum_wi_p += w * i
            else:
                sum_w_m += w
                sum_wi_m += w * i

        if sum_w_p:
            i_mean_p = sum_wi_p / sum_w_p
            sigi_mean_p = math.sqrt(1.0 / sum_w_p)
        else:
            i_mean_p = 0.0
            sigi_mean_p = 0.0

        if sum_w_m:
            i_mean_m = sum_wi_m / sum_w_m
            sigi_mean_m = math.sqrt(1.0 / sum_w_m)
        else:
            i_mean_m = 0.0
            sigi_mean_m = 0.0

        return i_mean_p, sigi_mean_p, i_mean_m, sigi_mean_m

    def multiplicity(self):
        return len(self._observations)

    def multiplicity_anomalous(self):
        def p(x): return x[0] % 2
        def m(x): return not x[0] % 2
        return len(filter(p, self._observations)), \
               len(filter(m, self._observations))

    def rmerge_contribution(self, i_mean):
        '''Calculate the contribution of this reflection to Rmerge.'''

        return sum([math.fabs(o[1] - i_mean) for o in self._observations])

    def rmerge_contribution_anomalous(self, i_mean_p, i_mean_m):
        '''Calculate the contribution of this reflection to Rmerge,
        separating anomalous pairs. Returns contributions for I+, I-.'''

        rmerge_p = 0.0
        rmerge_m = 0.0

        for o in self._observations:
            if o[0] % 2:
                rmerge_p += math.fabs(o[1] - i_mean_p)
            else:
                rmerge_m += math.fabs(o[1] - i_mean_m)

        return rmerge_p, rmerge_m

    def isigma_contribution(self):
        '''Calculate the contribution to the I/sigma. N.B. multiplicity!'''
        return sum(o[1] / o[2] for o in self._observations)

    def chisq_contribution(self, i_mean):
        '''Calculate the contribution to the reduced chi^2.'''

        return [(o[1] - i_mean) / o[2] for o in self._observations]

    def calculate_unmerged_di(self):
        '''Calculate unmerged dI values, returns a list of putative unmerged
        observations in same structure, though misym is replaced with dose
        difference...'''

        def isp(x): return x[0] % 2
        def ism(x): return not x[0] % 2

        i_p = filter(isp, self._observations)
        i_m = filter(ism, self._observations)

        dis = sorted([(math.fabs(m[3] - p[3]), ip, im) \
                      for ip, p in enumerate(i_p) \
                      for im, m in enumerate(i_m)])

        ip_used = []
        im_used = []

        result = unmerged_intensity()

        for d, ip, im in dis:
            if ip in ip_used:
                continue
            if im in im_used:
                continue

            ip_used.append(ip)
            im_used.append(im)

            # calculate intensity difference ...

            _ip = i_p[ip]
            _im = i_m[im]

            _i = _ip[1] - _im[1]
            _si = math.sqrt(_ip[2] * _ip[2] + _im[2] * _im[2])
            _b = 0.5 * (_ip[3] + _im[3])

            result.add(d, _i, _si, _b)

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

        hklin = args[0]

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
            
        self._mf = mtz_file(hklin)
        self._unmerged_reflections = { }
        self._merged_reflections = { }
        self._merged_reflections_anomalous = { }
        self._unmerged_di = { }

        all_columns = self._mf.get_column_names()

        assert('M_ISYM' in all_columns)
        assert('I' in all_columns)
        assert('SIGI' in all_columns)

        if 'DOSE' in all_columns and False:
            self._b_column = 'DOSE'
        elif 'BATCH' in all_columns:
            self._b_column = 'BATCH'
        else:
            raise RuntimeError, 'no baseline column (DOSE or BATCH) found'

        self._read_unmerged_reflections()
        self._merge_reflections()
        self._merge_reflections_anomalous()

        # self._calculate_unmerged_di()

        return

    def debug_info(self):
        '''Pull out some information for debugging, namely total intensity,
        number of reflections &c.'''

        n = 0
        I = 0.0

        for hkl in self._unmerged_reflections:
            for o in self._unmerged_reflections[hkl].get():
                n += 1
                I += o[1]

        return n, I

    def reload(self):
        '''Reload the reflection list &c.'''

        self._unmerged_reflections = { }
        self._merged_reflections = { }
        self._merged_reflections_anomalous = { }

        self._read_unmerged_reflections()
        self._merge_reflections()
        self._merge_reflections_anomalous()

        return

    def accumulate(self, other_r):
        '''Accumulate all of the measurements from another resolutionizer class
        instance.'''

        other_unmerged_reflections = other_r.get_unmerged_reflections()

        for hkl in other_unmerged_reflections:
            if not hkl in self._unmerged_reflections:
                self._unmerged_reflections[hkl] = unmerged_intensity()
            for observation in other_unmerged_reflections[hkl].get():
                m_isym, i, sigi, b = observation
                self._unmerged_reflections[hkl].add(m_isym, i, sigi, b)

        self._merge_reflections()
        self._merge_reflections_anomalous()

        return

    def _read_unmerged_reflections(self):
        '''Actually read the reflections in to memory.'''

        mi = self._mf.get_miller_indices()
        m_isym = self._mf.get_column_values('M_ISYM')
        i = self._mf.get_column_values('I')
        sigi = self._mf.get_column_values('SIGI')
        b = self._mf.get_column_values(self._b_column)

        for j in range(len(i)):
            hkl = mi[j]
            if not hkl in self._unmerged_reflections:
                self._unmerged_reflections[hkl] = unmerged_intensity()
            self._unmerged_reflections[hkl].add(
                m_isym[j], i[j], sigi[j], b[j])

        return

    def _merge_reflections(self):
        '''Merge the currently recorded unmerged reflections.'''

        for hkl in self._unmerged_reflections:
            self._merged_reflections[hkl] = self._unmerged_reflections[
                hkl].merge()

        return

    def _merge_reflections_anomalous(self):
        '''Merge the currently recorded unmerged reflections.'''

        for hkl in self._unmerged_reflections:
            self._merged_reflections_anomalous[
                hkl] = self._unmerged_reflections[hkl].merge_anomalous()

        return

    def _calculate_unmerged_di(self):
        '''Calculate a set of unmerged intensity differences.'''

        for hkl in self._unmerged_reflections:
            self._unmerged_di[hkl] = self._unmerged_reflections[
                hkl].calculate_unmerged_di()

        return

    def apply_kb(self, k, b):
        '''Apply kB scale factors to the recorded measurements, for all
        merged and unmerged observations.'''

        for hkl in self._merged_reflections:
            d = self.resolution(hkl)
            scale = k * math.exp(-1 * b / (d * d))
            i, sigi = self._merged_reflections[hkl]
            self._merged_reflections[hkl] = (i * scale, sigi * scale)

        for hkl in self._merged_reflections_anomalous:
            d = self.resolution(hkl)
            scale = k * math.exp(-1 * b / (d * d))
            ip, sigip, im, sigim = self._merged_reflections_anomalous[hkl]
            self._merged_reflections_anomalous[hkl] = (
                ip * scale, sigip * scale, im * scale, sigim * scale)

        for hkl in self._unmerged_reflections:
            d = self.resolution(hkl)
            scale = k * math.exp(-1 * b / (d * d))
            self._unmerged_reflections[hkl].apply_scale(scale)

        return

    def reindex(self, reindex_operation):
        '''Reindex the reflections by the given reindexing operation.'''

        R = rt_mx(reindex_operation).inverse()

        # first construct mapping table - native, anomalous

        map_native = { }

        hkls = flex.miller_index()

        for hkl in self._merged_reflections:
            Fhkl = R * hkl
            Rhkl = nint(Fhkl[0]), nint(Fhkl[1]), nint(Fhkl[2])
            hkls.append(Rhkl)

        map_to_asu(self._mf.get_space_group().type(), False, hkls)

        for j, hkl in enumerate(self._merged_reflections):
            map_native[hkl] = hkls[j]

        map_anomalous = { }

        hkls = flex.miller_index()

        for hkl in self._merged_reflections_anomalous:
            Fhkl = R * hkl
            Rhkl = nint(Fhkl[0]), nint(Fhkl[1]), nint(Fhkl[2])
            hkls.append(Rhkl)

        map_to_asu(self._mf.get_space_group().type(), True, hkls)

        for j, hkl in enumerate(self._merged_reflections_anomalous):
            map_anomalous[hkl] = hkls[j]

        # then remap the actual measurements

        merged_reflections = { }

        for hkl in self._merged_reflections:
            Rhkl = map_native[hkl]
            merged_reflections[Rhkl] = self._merged_reflections[hkl]

        self._merged_reflections = merged_reflections

        merged_reflections = { }

        for hkl in self._merged_reflections_anomalous:
            Rhkl = map_anomalous[hkl]
            merged_reflections[Rhkl] = self._merged_reflections_anomalous[hkl]

        self._merged_reflections_anomalous = merged_reflections

        unmerged_reflections = { }

        for hkl in self._unmerged_reflections:
            Rhkl = map_native[hkl]
            unmerged_reflections[Rhkl] = self._unmerged_reflections[hkl]

        self._unmerged_reflections = unmerged_reflections

        return

    def get_merged_reflections(self):
        return self._merged_reflections

    def get_unmerged_reflections(self):
        return self._unmerged_reflections

    def resolution(self, hkl):
        '''Compute the resolution corresponding to this miller index.'''

        return self._mf.get_unit_cell().d(hkl)

    def calculate_resolution_ranges(self, nbins = 20):
        '''Calculate semi-useful resolution ranges for analysis.'''

        miller_indices = list(self._merged_reflections)
        uc = self._mf.get_unit_cell()

        d_mi = []

        for mi in miller_indices:
            d = uc.d(mi)
            d_mi.append((d, mi))

        d_mi.sort()

        chunk_size = int(round(float(len(d_mi)) / nbins))

        hkl_ranges = []
        resolution_ranges = []

        for chunk in [d_mi[i:i + chunk_size] \
                      for i in range(0, len(d_mi), chunk_size)]:
            mi = [c[1] for c in chunk]
            d = [c[0] for c in chunk]
            hkl_ranges.append(mi)
            resolution_ranges.append((min(d), max(d)))

        # stitch together the two low res bins

        self._hkl_ranges = hkl_ranges[:-1]
        for mi in hkl_ranges[-1]:
            self._hkl_ranges[-1].append(mi)

        self._resolution_ranges = resolution_ranges[:-1]
        self._resolution_ranges[-1] = (self._resolution_ranges[-1][0],
                                       resolution_ranges[-1][1])

        return

    def get_resolution_bins(self):
        '''Return the reversed resolution limits - N.B. this is most
        important when considering resolution calculations, see
        resolution_completeness.'''

        return list(reversed(self._hkl_ranges)), \
               list(reversed(self._resolution_ranges))

    def apply_resolution_limit(self, dmin):
        '''Remove reflections with resolution < dmin.'''

        delete = []
        for hkl in self._merged_reflections:
            if self.resolution(hkl) < dmin:
                delete.append(hkl)

        for hkl in delete:
            del(self._merged_reflections[hkl])

        delete = []
        for hkl in self._merged_reflections_anomalous:
            if self.resolution(hkl) < dmin:
                delete.append(hkl)

        for hkl in delete:
            del(self._merged_reflections_anomalous[hkl])

        delete = []
        for hkl in self._unmerged_reflections:
            if self.resolution(hkl) < dmin:
                delete.append(hkl)

        for hkl in delete:
            del(self._unmerged_reflections[hkl])

        return

    def calculate_completeness(self, resolution_bin = None):
        '''Calculate the completeness of observations in a named
        resolution bin.'''

        if resolution_bin is None:
            resolution_range = self._mf.get_resolution_range()
            hkl_list = list(self._merged_reflections)
        else:
            resolution_range = self._resolution_ranges[resolution_bin]
            hkl_list = self._hkl_ranges[resolution_bin]

        uc = self._mf.get_unit_cell()
        sg = self._mf.get_space_group()

        dmin = min(resolution_range)
        dmax = max(resolution_range)

        cs = crystal_symmetry(unit_cell = uc, space_group = sg)
        hkl_calc = [hkl for hkl in build_set(
            cs, False, d_min = dmin, d_max = dmax).indices()]

        # remove systematically absent reflections

        hkl_list = [hkl for hkl in
                    itertools.ifilterfalse(sg.is_sys_absent, hkl_list)]

        return float(len(hkl_list)) / float(len(hkl_calc))

    def calculate_rmerge(self, hkl_list = None):
        '''Calculate the overall Rmerge.'''

        t = 0.0
        b = 0.0

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        for hkl in hkl_list:
            # if we have only one observation, do not include in the
            # rmerge calculations
            if self._unmerged_reflections[hkl].multiplicity() > 1:
                i_mean = self._merged_reflections[hkl][0]
                t += self._unmerged_reflections[hkl].rmerge_contribution(
                    i_mean)
                b += self._unmerged_reflections[hkl].multiplicity() * i_mean

        if not b:
            return 0.0

        return t / b

    def calculate_rmerge_anomalous(self, hkl_list = None):
        '''Calculate the overall Rmerge, separating anomalous pairs.'''

        t = 0.0
        b = 0.0

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        for hkl in hkl_list:
            is_pm = self._merged_reflections_anomalous[hkl]
            i_mean_p, i_mean_m = is_pm[0], is_pm[2]
            r_pm = self._unmerged_reflections[
                hkl].rmerge_contribution_anomalous(i_mean_p, i_mean_m)
            t += r_pm[0] + r_pm[1]

            # if we have only one observation, do not include in the
            # rmerge calculations

            mult_p, mult_m = self._unmerged_reflections[
                hkl].multiplicity_anomalous()

            if mult_p == 1:
                mult_p = 0
            if mult_m == 1:
                mult_m = 0

            b += mult_p * i_mean_p + mult_m * i_mean_m

        if not b:
            return 0.0

        return t / b

    def calculate_chisq(self, hkl_list = None):
        '''Calculate the overall erzatz chi^2.'''

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        deltas = []

        for hkl in hkl_list:
            i_mean = self._merged_reflections[hkl][0]
            for d in self._unmerged_reflections[hkl].chisq_contribution(
                i_mean):
                deltas.append(d)

        mean = sum(deltas) / len(deltas)
        var = sum([(d - mean) * (d - mean) for d in deltas]) / len(deltas)

        return mean, math.sqrt(var)

    def calculate_multiplicity(self, hkl_list = None):
        '''Calculate the overall average multiplicity.'''

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        multiplicity = [float(self._unmerged_reflections[hkl].multiplicity()) \
                        for hkl in hkl_list]

        return sum(multiplicity) / len(multiplicity)

    def calculate_merged_isigma(self, hkl_list = None):
        '''Calculate the average merged I/sigma.'''

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        isigma_values = [self._merged_reflections[hkl][0] / \
                         self._merged_reflections[hkl][1] \
                         for hkl in hkl_list]

        return sum(isigma_values) / len(isigma_values)

    def calculate_unmerged_isigma(self, hkl_list = None):
        '''Calculate the average unmerged I/sigma.'''

        if not hkl_list:
            hkl_list = list(self._unmerged_reflections)

        return sum([self._unmerged_reflections[hkl].isigma_contribution() \
                    for hkl in hkl_list]) / \
                    sum([self._unmerged_reflections[hkl].multiplicity() \
                         for hkl in hkl_list])

    def calculate_z2(self, hkl_list = None):
        '''Calculate average Z^2 values, where Z = I/<I> in the bin,
        from the merged observations. Now also separate centric and
        acentric reflections.'''

        if not hkl_list:
            hkl_list = list(self._merged_reflections)

        # separate centric and acentric reflections

        sg = self._mf.get_space_group()

        hkl_centric = [hkl for hkl in
                       itertools.ifilter(sg.is_centric, hkl_list)]
        hkl_acentric = [hkl for hkl in
                        itertools.ifilterfalse(sg.is_centric, hkl_list)]

        i_s = [self._merged_reflections[hkl][0] for hkl in hkl_centric]
        mean_i = sum(i_s) / len(i_s)

        z_s = [i / mean_i for i in i_s]
        z_centric = sum([z * z for z in z_s]) / len(z_s)

        i_s = [self._merged_reflections[hkl][0] for hkl in hkl_acentric]
        mean_i = sum(i_s) / len(i_s)

        z_s = [i / mean_i for i in i_s]
        z_acentric = sum([z * z for z in z_s]) / len(z_s)

        return z_centric, z_acentric

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


            
    def resolution_rmerge(self, limit = None, log = None):
        '''Compute a resolution limit where either rmerge = 1.0 (limit if
        set) or the full extent of the data. N.B. this fit is only meaningful
        for positive values.'''

        if limit is None:
            limit = self._params.rmerge

        bins, ranges = self.get_resolution_bins()

        if limit == 0.0:
            return ranges[-1][0]

        rmerge_s = get_positive_values(
            [self.calculate_rmerge(bin) for bin in bins])

        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(rmerge_s)]

        if limit == 0.0:
            return 1.0 / math.sqrt(max(s_s))

        if limit > max(rmerge_s):
            return 1.0 / math.sqrt(max(s_s))

        rmerge_f = log_inv_fit(s_s, rmerge_s, 6)

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(s_s):
                d = 1.0 / math.sqrt(s)
                o = rmerge_s[j]
                m = rmerge_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_rmerge = 1.0 / math.sqrt(interpolate_value(s_s, rmerge_f, limit))
        except:
            r_rmerge = 1.0 / math.sqrt(max(s_s))

        return r_rmerge

    def new_resolution_unmerged_isigma(self, limit = None, log = None):
        '''Compute a resolution limit where either I/sigma = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.isigma

        bins, ranges = self.get_resolution_bins()

        isigma_s = get_positive_values(
            [self.calculate_unmerged_isigma(bin) for bin in bins])

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

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(_s_s):
                d = 1.0 / math.sqrt(s)
                o = _isigma_s[j]
                m = _isigma_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_isigma = 1.0 / math.sqrt(interpolate_value(_s_s, _isigma_f,
                                                         limit))
        except:
            r_isigma = 1.0 / math.sqrt(max(_s_s))

        return r_isigma

    def resolution_unmerged_isigma(self, limit = None, log = None):
        '''Compute a resolution limit where either I/sigma = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.isigma

        bins, ranges = self.get_resolution_bins()

        isigma_s = get_positive_values(
            [self.calculate_unmerged_isigma(bin) for bin in bins])

        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(isigma_s)]

        if min(isigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        isigma_f = log_fit(s_s, isigma_s, 6)

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(s_s):
                d = 1.0 / math.sqrt(s)
                o = isigma_s[j]
                m = isigma_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_isigma = 1.0 / math.sqrt(interpolate_value(s_s, isigma_f, limit))
        except:
            r_isigma = 1.0 / math.sqrt(max(s_s))

        return r_isigma

    def new_resolution_merged_isigma(self, limit = None, log = None):
        '''Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.misigma

        bins, ranges = self.get_resolution_bins()

        misigma_s = get_positive_values(
            [self.calculate_merged_isigma(bin) for bin in bins])
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

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(_s_s):
                d = 1.0 / math.sqrt(s)
                o = _misigma_s[j]
                m = _misigma_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_misigma = 1.0 / math.sqrt(interpolate_value(_s_s, _misigma_f,
                                                          limit))
        except:
            r_misigma = 1.0 / math.sqrt(max(_s_s))

        return r_misigma

    def resolution_merged_isigma(self, limit = None, log = None):
        '''Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
        set) or the full extent of the data.'''

        if limit is None:
            limit = self._params.misigma

        bins, ranges = self.get_resolution_bins()

        misigma_s = get_positive_values(
            [self.calculate_merged_isigma(bin) for bin in bins])
        s_s = [1.0 / (r[0] * r[0]) for r in ranges][:len(misigma_s)]

        if min(misigma_s) > limit:
            return 1.0 / math.sqrt(max(s_s))

        misigma_f = log_fit(s_s, misigma_s, 6)

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(s_s):
                d = 1.0 / math.sqrt(s)
                o = misigma_s[j]
                m = misigma_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_misigma = 1.0 / math.sqrt(
                interpolate_value(s_s, misigma_f, limit))
        except:
            r_misigma = 1.0 / math.sqrt(max(s_s))

        return r_misigma

    def resolution_completeness(self, limit = None, log = None):
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

        if log:
            fout = open(log, 'w')
            for j, s in enumerate(s_s):
                d = 1.0 / math.sqrt(s)
                o = comp_s[j]
                m = comp_f[j]
                fout.write('%f %f %f %f\n' % (s, d, o, m))
            fout.close()

        try:
            r_comp = 1.0 / math.sqrt(
                interpolate_value(s_s, comp_f, rlimit))
        except:
            r_comp = 1.0 / math.sqrt(max(s_s))

        return r_comp

if __name__ == '__main__':

    m = resolutionizer(sys.argv[1:])
    m.resolution_auto()
