#!/usr/bin/env python
#
# dials.merge_spot_lists.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class MergeSpotLists(object):
  ''' Class to merge spot lists '''

  def __init__(self):
    pass

  def __call__(self, spotlist):
    ''' Merge the spot lists

    Params:
        spotlist The spotlist

    Returns:
        The list of spot shoeboxes

    '''
    from dials.util.command_line import ProgressBar, Command
    from dials.algorithms.image.connected_components import LabelPixels3d
    from dials.algorithms.shoebox import MaskCode
    from dials.array_family import flex
    from dials.model.data import ReflectionList

    # Construct the pixel labeller
    Command.start('Extracting pixels from spot list')
    maxx = max([s.bounding_box[1] for s in spotlist])
    maxy = max([s.bounding_box[3] for s in spotlist])
    maxz = max([s.bounding_box[5] for s in spotlist])
    label = LabelPixels3d((maxz, maxy, maxx))
    for s in spotlist:
      assert(s.panel_number == 0)
      bbox = s.bounding_box
      sbox = s.shoebox
      mask = s.shoebox_mask
      coords = []
      values = []
      for k in range(bbox[5] - bbox[4]):
        for j in range(bbox[3] - bbox[2]):
          for i in range(bbox[1] - bbox[0]):
            if mask[k,j,i] & MaskCode.Valid:
              coords.append((i+bbox[0], j+bbox[2], k+bbox[4]))
              values.append(int(sbox[k,j,i]))
      label.add_pixels(flex.int(values), flex.vec3_int(coords))
    Command.end('Extracted pixels from spot list')

    # Extract the shoeboxes
    Command.start('Extracting spots from pixels')
    shoeboxes = flex.shoebox(label, i)
    Command.end('Extracted {0} spots from pixels'.format(len(shoeboxes)))

    # Calculate the spot centroids
    Command.start('Calculating {0} spot centroids'.format(len(shoeboxes)))
    centroid = shoeboxes.centroid_valid()
    Command.end('Calculated {0} spot centroids'.format(len(shoeboxes)))

    # Calculate the spot intensities
    Command.start('Calculating {0} spot intensities'.format(len(shoeboxes)))
    intensity = shoeboxes.summed_intensity_valid()
    Command.end('Calculated {0} spot intensities'.format(len(shoeboxes)))

    # Create the observations
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    # Return as a reflection list
    return ReflectionList(observed, shoeboxes)


class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage = "usage: %prog [options] [param.phil] "\
            "{reflection1.file [reflection2.file ...]}"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
        '-o', '--output-filename',
        dest = 'output_filename',
        type = 'string', default = 'merged.pickle',
        help = 'Set the filename for matched reflections.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.util.command_line import Importer, Command
    from dials.model.serialize import dump

    # Try importing the command line arguments
    importer = Importer(args)
    if len(importer.reflections) == 0:
      self.config().print_help()
      return
    reflections = importer.reflections[0]
    for rlist in importer.reflections:
      reflections.extend(rlist)

    merge = MergeSpotLists()
    spotlist = merge(reflections)

    # Save the reflections to file
    Command.start('Saving {0} spot list to {1}'.format(
        len(spotlist), options.output_filename))
    dump.reflections(spotlist, options.output_filename)
    Command.end('Saved {0} spot list to {1}'.format(
        len(spotlist), options.output_filename))




if __name__ == '__main__':
  script = Script()
  script.run()
