from __future__ import division
from cctbx.array_family import flex

class TestForOverlaps(object):
  from dials.algorithms.shoebox import MaskCode

  code_fgd = MaskCode.Foreground | MaskCode.Valid
  @staticmethod
  def is_fgd(code):
    return (code & TestForOverlaps.code_fgd) == TestForOverlaps.code_fgd

  code_bgd = MaskCode.Background | MaskCode.Valid
  @staticmethod
  def is_bgd(code):
    return (code & TestForOverlaps.code_bgd) == TestForOverlaps.code_bgd

  code_overlap = code_fgd | code_bgd
  @staticmethod
  def is_overlap(code):
    return (code & TestForOverlaps.code_overlap) == TestForOverlaps.code_overlap

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser, Importer, flatten_experiments, flatten_reflections
    import libtbx.load_env

    # The script usage
    usage = "usage: %s experiments.json reflections.pickle" \
            % libtbx.env.dispatcher_name

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      read_experiments=True,
      read_reflections=True)

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    self.experiments = flatten_experiments(params.input.experiments)
    self.reflections = flatten_reflections(params.input.reflections)

  def run(self):
    for expt, refl in zip(self.experiments, self.reflections):
      det = expt.detector
      size_fast, size_slow = det[0].get_image_size()
      mask_array = flex.size_t(size_fast*size_slow)
      for obs in refl:
        shoebox = obs['shoebox']
        fast_coords = xrange(shoebox.xsize())
        slow_coords = xrange(shoebox.ysize())
        for f, s in zip(fast_coords, slow_coords):
          f_abs = f + shoebox.bbox[0] # relative to detector
          s_abs = s + shoebox.bbox[2] # relative to detector
          posn = f_abs + s_abs*size_fast # position in mask_array
          posn_in_shoebox = f + shoebox.xsize()*s # position in shoebox
          if self.is_fgd(shoebox.mask[posn_in_shoebox]) and self.is_fgd(mask_array[posn]):
            print "Overlapping foreground found at indexed position (%d, %d), " % (f_abs, s_abs) \
            + "observed centroid (%d, %d)" % (obs['xyzcal.px'][0], obs['xyzcal.px'][1])
            assert False
          try:
            mask_array[posn] |= shoebox.mask[posn_in_shoebox]
          except IndexError:
            continue
      for i in xrange(len(mask_array)):
        this_code = mask_array[i]
        if self.is_overlap(this_code):
          f = i % shoebox.xsize()
          s = i // shoebox.xsize()
          print "Overlapping foreground and background found at (%d, %d)" % (f, s)
          assert False
    print "OK"

class Test(TestForOverlaps):

  def __init__(self):
    from dials.util.options import Importer, flatten_experiments, flatten_reflections
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # test data
    import os
    refl_path = os.path.join(dials_regression,
      'integration_test_data', 'stills_PSII', 'idx-20161021225550223_integrated.pickle')
    expt_path = os.path.join(dials_regression,
      'integration_test_data', 'stills_PSII', 'idx-20161021225550223_refined_experiments.json')

    importer = Importer([refl_path, refl_path, expt_path, expt_path], read_experiments=True, read_reflections=True, check_format=False)

    self.reflections = flatten_reflections(importer.reflections)
    self.experiments = flatten_experiments(importer.experiments)

    from dials.algorithms.integration.overlaps_filter import OverlapsFilterMultiExpt
    overlaps_filter = OverlapsFilterMultiExpt(self.reflections[0], self.experiments)
    overlaps_filter.remove_foreground_foreground_overlaps()
    overlaps_filter.remove_foreground_background_overlaps()
    self.reflections = [overlaps_filter.refl]

if __name__ == '__main__':
  test = Test()
  test.run()
