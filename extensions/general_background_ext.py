#!/usr/bin/env python
#
# general_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class GeneralBackgroundExt(BackgroundIface):
  ''' An extension class implementing XDS background subtraction. '''

  name = 'general'

  phil = '''
    outlier
      .help = "Outlier rejection prior to background fit"
    {
      algorithm = null *nsigma truncated normal
        .help = "The outlier rejection algorithm."
        .type = choice

      nsigma
        .help = "Parameters for nsigma outlier rejector"
      {
        lower = 3
          .help = "Lower n sigma"
          .type = float
        upper = 3
          .help = "Upper n sigma"
          .type = float
      }

      truncated
        .help = "Parameters for truncated outlier rejector"
      {
        lower = 0.01
          .help = "Lower bound"
          .type = float
        upper = 0.01
          .help = "Upper bound"
          .type = float
      }

      normal
        .help = "Parameters for normal outlier rejector"
      {
        min_pixels = 10
          .help = "The minimum number of pixels to use in calculating the"
                  "background intensity."
          .type = int
      }
    }

    model
      .help = "Background model"
    {
      algorithm = constant2d *constant3d linear2d linear3d
        .help = "The choice of background model"
        .type = choice
    }

  '''

  default=True

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import Creator
    from dials.algorithms.background import TruncatedOutlierRejector
    from dials.algorithms.background import NSigmaOutlierRejector
    from dials.algorithms.background import NormalOutlierRejector
    from dials.algorithms.background import Constant2dModeller
    from dials.algorithms.background import Constant3dModeller
    from dials.algorithms.background import Linear2dModeller
    from dials.algorithms.background import Linear3dModeller

    def select_modeller():
      model = params.integration.background.general.model
      if model.algorithm == 'constant2d':
        return Constant2dModeller()
      elif model.algorithm == 'constant3d':
        return Constant3dModeller()
      elif model.algorithm == 'linear2d':
        return Linear2dModeller()
      elif model.algorithm == 'linear3d':
        return Linear3dModeller()
      raise RuntimeError("Unexpected background model %s" % model.algorithm)

    def select_rejector():
      outlier = params.integration.background.general.outlier
      if outlier.algorithm == 'null':
        return None
      elif outlier.algorithm == 'truncated':
        return TruncatedOutlierRejector(
          outlier.truncated.lower,
          outlier.truncated.upper)
      elif outlier.algorithm == 'nsigma':
        return NSigmaOutlierRejector(
          outlier.nsigma.lower,
          outlier.nsigma.upper)
      elif outlier.algorithm == 'normal':
        return NormalOutlierRejector(
          outlier.normal.min_pixels)
      raise RuntimeError("Unexpected outlier rejector %s" % outlier.algorithm)

    modeller = select_modeller()
    rejector = select_rejector()
    self._creator = Creator(modeller, rejector)

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    from dials.util.command_line import Command
    from dials.array_family import flex

    # Do the background subtraction
    Command.start('Calculating reflection background')
    reflections['background.mse'] = flex.double(len(reflections))
    success = self._creator(
      reflections['shoebox'],
      reflections['background.mse'])
    reflections['background.mean'] = reflections['shoebox'].mean_background()
    reflections.del_selected(success != True)
    Command.end('Calculated {0} background values'.format(len(reflections)))
