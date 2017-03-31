#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.extract_background_signal

from __future__ import absolute_import, division

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

help_message = '''
%s experiments.json reflections.pickle
''' % libtbx.env.dispatcher_name


from libtbx.phil import parse
phil_scope = parse('''
  padding = 0
    .type = int(value_min=0)
    .help = "Add padding around shoebox"
''')

class Script(object):
  ''' Class to run the script. '''

  def __init__(self):
    ''' Initialise the script. '''

    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiment.json reflections.pickle" % \
      libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True,
      read_reflections=True,
      read_datablocks=True)

  def run(self):
    ''' Extract the shoeboxes. '''
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments
    from dials.util.options import flatten_datablocks
    from dials.util import log
    from dials.array_family import flex
    from libtbx.utils import Sorry

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Configure logging
    log.config()

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Get the data
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    if not any([experiments, datablocks, reflections]):
      self.parser.print_help()
      exit(0)
    elif experiments and datablocks:
      raise Sorry('Both experiment list and datablocks set')
    elif len(experiments) > 1:
      raise Sorry('More than 1 experiment set')
    elif len(datablocks) > 1:
      raise Sorry('More than 1 datablock set')
    elif len(experiments) == 1:
      imageset = experiments[0].imageset
    elif len(datablocks) == 1:
      imagesets = datablocks[0].extract_imagesets()
      if len(imagesets) != 1:
        raise Sorry('Need 1 imageset, got %d' % len(imagesets))
      imageset = imagesets[0]
    if len(reflections) != 1:
      raise Sorry('Need 1 reflection table, got %d' % len(reflections))
    else:
      reflections = reflections[0]

    # Check the reflections contain the necessary stuff
    assert("bbox" in reflections)
    assert("panel" in reflections)

    # Get some models
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    frame0, frame1 = scan.get_array_range()

    # Add some padding but limit to image volume
    if params.padding > 0:
      logger.info('Adding %d pixels as padding' % params.padding)
      x0, x1, y0, y1, z0, z1 = reflections['bbox'].parts()
      x0 -= params.padding
      x1 += params.padding
      y0 -= params.padding
      y1 += params.padding
      panel = reflections['panel']
      for i in range(len(reflections)):
        width, height = detector[panel[i]].get_image_size()
        if x0[i] < 0: x0[i] = 0
        if x1[i] > width: x1[i] = width
        if y0[i] < 0: y0[i] = 0
        if y1[i] > height: y1[i] = height
      reflections['bbox'] = flex.int6(x0, x1, y0, y1, z0, z1)

    # Now iterate through the images, masking the shoebox pixels as necessary,
    # then model the background & compute intensity statistics

    x0, x1, y0, y1, z0, z1 = reflections['bbox'].parts()

    for frame in range(frame0, frame1):
      data = imageset.get_raw_data(frame)[0]
      mask = (data < 0)
      data = data.as_double()
      sel = (z0 <= frame) & (z1 >= frame)
      _x0 = x0.select(sel)
      _x1 = x1.select(sel)
      _y0 = y0.select(sel)
      _y1 = y1.select(sel)
      for __x0, __x1, __y0, __y1 in zip(_x0, _x1, _y0, _y1):
        for y in range(__y0, __y1):
          for x in range(__x0, __x1):
            data[y, x] = 0.0
            mask[y, x] = True
      imask = (~mask).as_1d().as_int()
      imask.reshape(data.accessor())

      from dials.algorithms.image.filter import summed_area

      # smooth over scale corresponding to max shoebox size
      d = max(flex.max(_x1 - _x0), flex.max(_y1 - _y0))

      summed_background = summed_area(data, (d,d))
      summed_mask = summed_area(imask, (d,d))
      mean_background = (summed_background /
                         summed_mask.as_double())
      data.as_1d().set_selected(mask.as_1d(), mean_background.as_1d())
      print flex.sum(data.select(data.as_1d() > 0))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
