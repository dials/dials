#
# statistics.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.array_family import flex
from dials.array_family.flex import Binner


def flex_ios(val, var):
  ''' Compute I/sigma or return zero for each element. '''
  assert(var.all_ge(0))
  assert(len(val) == len(var))
  result = flex.double(len(val),0)
  indices = flex.size_t(range(len(val))).select(var > 0)
  val = val.select(indices)
  var = var.select(indices)
  assert(var.all_gt(0))
  result.set_selected(indices, val / flex.sqrt(var))
  return result


class ImageSummary(object):
  ''' A class to produce statistics per image. '''

  def __init__(self, data, experiment):
    ''' Compute stats. '''

    # Check some table columns
    assert("flags" in data)
    assert("bbox" in data)
    assert("partiality" in data)
    assert("intensity.sum.value" in data)
    assert("intensity.sum.variance" in data)

    # Get the array range
    try:
      array_range = experiment.imageset.get_array_range()
    except:
      array_range = (0, len(experiment.imageset))

    # Get arrays for each frame
    data = data.select(data.get_flags(data.flags.integrated, all=False))
    data.split_partials()
    frames = data['bbox'].parts()[4]

    # Create the binner with the bins per image
    binner = Binner(flex.int(range(
      array_range[0],
      array_range[1]+1)).as_double())

    # Get the bins
    self.bins = binner.bins()

    # Get full and partial counts
    full = data['partiality'] > 0.997300203937
    bin_indexer = binner.indexer(frames.as_double())
    self.full = bin_indexer.sum(full.as_double())
    self.part = bin_indexer.sum((~full).as_double())

    # Get the mean background values
    try:
      i_flg = data.get_flags(data.flags.integrated, all=False)
      i_bg = data['background.mean'].select(i_flg)
      bin_indexer = binner.indexer(frames.select(i_flg).as_double())
      self.i_bg = bin_indexer.mean(i_bg)
    except RuntimeError:
      self.i_bg = flex.double(len(self.bins), 0)

    # Get stuff from table for summation
    i_sum_flg = data.get_flags(data.flags.integrated_sum)
    i_sum_val = data['intensity.sum.value'].select(i_sum_flg)
    i_sum_var = data['intensity.sum.variance'].select(i_sum_flg)
    ios_sum = flex_ios(i_sum_val, i_sum_var)
    bin_indexer = binner.indexer(frames.select(i_sum_flg).as_double())
    self.num_sum = bin_indexer.count()
    self.ios_sum = bin_indexer.mean(ios_sum)

    # Get stuff from table for profile fitting
    try:
      i_prf_flg = data.get_flags(data.flags.integrated_prf)
      i_prf_val = data['intensity.prf.value'].select(i_prf_flg)
      i_prf_var = data['intensity.prf.variance'].select(i_prf_flg)
      cor_prf = data['profile.correlation'].select(i_prf_flg)
      ios_prf = flex_ios(i_prf_val, i_prf_var)
      bin_indexer = binner.indexer(frames.select(i_prf_flg).as_double())
      self.num_prf = bin_indexer.count()
      self.ios_prf = bin_indexer.mean(ios_prf)
      self.cor_prf = bin_indexer.mean(cor_prf)
    except RuntimeError:
      self.num_prf = flex.size_t(len(self.bins), 0)
      self.ios_prf = flex.double(len(self.bins), 0)
      self.cor_prf = flex.double(len(self.bins), 0)

  def __len__(self):
    ''' The number of bins. '''
    assert(len(self.bins) > 1)
    return len(self.bins)-1

  def table(self):
    ''' Produce a table of results. '''
    from libtbx.table_utils import format as table
    rows = [["Image",
             "# full",
             "# part",
             "# sum",
             "# prf",
             "<Ibg>",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)",
             "<CC prf>"]]
    for i in range(len(self)):
      rows.append([
        '%d' % self.bins[i],
        '%d' % self.full[i],
        '%d' % self.part[i],
        '%d' % self.num_sum[i],
        '%d' % self.num_prf[i],
        '%.2f' % self.i_bg[i],
        '%.2f' % self.ios_sum[i],
        '%.2f' % self.ios_prf[i],
        '%.2f' % self.cor_prf[i]])
    return table(rows, has_header=True, justify='right', prefix=' ')


class ResolutionSummary(object):
  ''' A class to produce statistics in resolution shells. '''

  def __init__(self, data, experiment, nbins=20):
    ''' Compute the statistics. '''
    from cctbx import miller
    from cctbx import crystal

    # Check some table columns
    assert("flags" in data)
    assert("d" in data)
    assert("intensity.sum.value" in data)
    assert("intensity.sum.variance" in data)

    # Select integrated reflections
    data = data.select(data.get_flags(data.flags.integrated, all=False))

    # Create the crystal symmetry object
    cs = crystal.symmetry(
      space_group=experiment.crystal.get_space_group(),
      unit_cell=experiment.crystal.get_unit_cell())

    # Create the array of bins
    ms = miller.set(cs, data['miller_index'])
    ms.setup_binner(n_bins=nbins)
    binner = ms.binner()
    brange = list(binner.range_used())
    bins = [binner.bin_d_range(brange[0])[0]]
    for i in brange:
      bins.append(binner.bin_d_range(i)[1])
    bins = flex.double(reversed(bins))

    # Create the binner
    binner = Binner(bins)

    # Get the bins
    self.bins = binner.bins()

    # Get full and partial counts
    full = data['partiality'] > 0.997300203937
    over = data.get_flags(data.flags.overloaded)
    ice = data.get_flags(data.flags.in_powder_ring)
    bin_indexer = binner.indexer(data['d'])
    self.num_full = bin_indexer.sum(full.as_double())
    self.num_part = bin_indexer.sum((~full).as_double())
    self.num_over = bin_indexer.sum(over.as_double())
    self.num_ice  = bin_indexer.sum(ice.as_double())

    # Get the mean background values
    try:
      i_flg = data.get_flags(data.flags.integrated, all=False)
      i_bg = data['background.mean'].select(i_flg)
      bin_indexer = binner.indexer(data['d'].select(i_flg).as_double())
      self.i_bg = bin_indexer.mean(i_bg)
    except RuntimeError:
      self.i_bg = flex.double(len(self.bins), 0)

    # Get stuff from table for summation
    i_sum_flg = data.get_flags(data.flags.integrated_sum)
    i_sum_val = data['intensity.sum.value'].select(i_sum_flg)
    i_sum_var = data['intensity.sum.variance'].select(i_sum_flg)
    ios_sum = flex_ios(i_sum_val, i_sum_var)
    bin_indexer = binner.indexer(data['d'].select(i_sum_flg))
    self.num_sum = bin_indexer.count()
    self.ios_sum = bin_indexer.mean(ios_sum)

    # Get stuff from table for profile fitting
    try:
      i_prf_flg = data.get_flags(data.flags.integrated_prf)
      i_prf_val = data['intensity.prf.value'].select(i_prf_flg)
      i_prf_var = data['intensity.prf.variance'].select(i_prf_flg)
      cor_prf = data['profile.correlation'].select(i_prf_flg)
      ios_prf = flex_ios(i_prf_val, i_prf_var)
      bin_indexer = binner.indexer(data['d'].select(i_prf_flg))
      self.num_prf = bin_indexer.count()
      self.ios_prf = bin_indexer.mean(ios_prf)
      self.cor_prf = bin_indexer.mean(cor_prf)
    except Exception:
      self.num_prf = flex.size_t(len(bins)-1, 0)
      self.ios_prf = flex.double(len(bins)-1, 0)
      self.cor_prf = flex.double(len(bins)-1, 0)

  def __len__(self):
    ''' The number of bins. '''
    assert(len(self.bins) > 1)
    return len(self.bins)-1

  def table(self):
    ''' Produce a table. '''
    from libtbx.table_utils import format as table
    rows = [["d min",
             "d max",
             "# full",
             "# part",
             "# over",
             "# ice",
             "# sum",
             "# prf",
             "<Ibg>",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)",
             "<CC prf>"]]
    for i in range(len(self)):
      rows.append([
        '%.2f' % self.bins[i],
        '%.2f' % self.bins[i+1],
        '%d'   % self.num_full[i],
        '%d'   % self.num_part[i],
        '%d'   % self.num_over[i],
        '%d'   % self.num_ice[i],
        '%d'   % self.num_sum[i],
        '%d'   % self.num_prf[i],
        '%.2f' % self.i_bg[i],
        '%.2f' % self.ios_sum[i],
        '%.2f' % self.ios_prf[i],
        '%.2f' % self.cor_prf[i]])
    return table(rows, has_header=True, justify='right', prefix=' ')


class WholeSummary(object):
  ''' A class to produce statistics for the whole dataset. '''

  def __init__(self, data, experiment):
    ''' Compute the results. '''

    # Check columns are there
    assert("partiality" in data)
    assert("intensity.sum.value" in data)
    assert("intensity.sum.variance" in data)

    # Compute some flag stuff
    full = data['partiality'] > 0.997300203937
    over = data.get_flags(data.flags.overloaded)
    ice = data.get_flags(data.flags.in_powder_ring)
    self.num_full = full.count(True)
    self.num_part = full.count(False)
    self.num_over = over.count(True)
    self.num_ice = ice.count(True)

    # Get the mean background values
    try:
      i_flg = data.get_flags(data.flags.integrated, all=False)
      self.i_bg = flex.mean(data['background.mean'].select(i_flg))
    except RuntimeError:
      self.i_bg = 0

    # Compute for summation
    flags_sum = data.get_flags(data.flags.integrated_sum)
    I_sum_val = data['intensity.sum.value'].select(flags_sum)
    I_sum_var = data['intensity.sum.variance'].select(flags_sum)
    self.ios_sum = flex.mean(flex_ios(I_sum_val, I_sum_var))
    self.num_sum = flags_sum.count(True)

    # Compute for profile fitting
    try:
      flags_prf = data.get_flags(data.flags.integrated_prf)
      I_prf_val = data['intensity.prf.value'].select(flags_prf)
      I_prf_var = data['intensity.prf.variance'].select(flags_prf)
      self.ios_prf = flex.mean(flex_ios(I_prf_val, I_prf_var))
      self.num_prf = flags_prf.count(True)
      self.cor_prf = flex.mean(data['profile.correlation'])
    except Exception:
        self.ios_prf = 0.0
        self.num_prf = 0
        self.cor_prf = 0

  def table(self):
    ''' Produce a table of results. '''
    from libtbx.table_utils import format as table
    rows = [["Number fully recorded",                 '%d'   % self.num_full],
            ["Number partially recorded",             '%d'   % self.num_part],
            ["Number with overloaded pixels",         '%d'   % self.num_over],
            ["Number in powder rings",                '%d'   % self.num_ice],
            ["Number processed with summation",       '%d'   % self.num_sum],
            ["Number processed with profile fitting", '%d'   % self.num_prf],
            ["<Ibg>",                                 '%.2f' % self.i_bg],
            ["<I/sigI> (summation)",                  '%.2f' % self.ios_sum],
            ["<I/sigI> (profile fitting)",            '%.2f' % self.ios_prf],
            ["<CC prf>",                              '%.2f' % self.cor_prf]]
    return table(rows, justify='left', prefix=' ')


class Summary(object):
  ''' A class to present a summary of integration results. '''

  def __init__(self, index, data, experiment):
    ''' Initialise. '''
    self._index = index
    self._image_summary = ImageSummary(data, experiment)
    self._resolution_summary = ResolutionSummary(data, experiment)
    self._whole_summary = WholeSummary(data, experiment)

  def __str__(self):
    ''' Return as a string. '''
    from dials.util.command_line import heading
    img_summary = self._image_summary.table()
    res_summary = self._resolution_summary.table()
    who_summary = self._whole_summary.table()
    return (
      '%s\n'
      '\n'
      '%s\n'
      '\n'
      ' Summary of integration results as a function of image number'
      '\n%s\n\n'
      ' Summary of integration results binned by resolution'
      '\n%s\n\n'
      ' Summary of integration results for the whole dataset'
      '\n%s\n'
    ) % ('=' * 80,
         heading('Summary of integration results for experiment %d' %
                 self._index),
         img_summary,
         res_summary,
         who_summary)


def statistics(data, experiments):
  ''' Return some simple statistics for a reflection table. '''
  tables = data.split_by_experiment_id()
  assert(len(tables) == len(experiments))
  summaries = []
  for index, (table, experiment) in enumerate(zip(tables, experiments)):
    summaries.append(Summary(index, table, experiment))
  return summaries
