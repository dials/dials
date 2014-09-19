#
# algorithm.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class IntegrationAlgorithm(object):
  ''' A class to perform profile fitting '''

  def __init__(self,
               experiments,
               profile_model,
               grid_size=5,
               debug=False):
    '''Initialise algorithm.'''
    assert(len(experiments) == len(profile_model))
    self._experiments = experiments
    self._profile_model = profile_model
    self._grid_size = grid_size
    self._debug = debug

  def __call__(self, reflections):
    '''Process the reflections.

    Params:
        reflections The reflections to process

    Returns:
        The list of integrated reflections

    '''
    from dials.algorithms.integration.fit_image import ImageSpaceProfileFitting
    from dials.algorithms.integration.fit_image import Spec
    from dials.algorithms.integration.interface import job_id
    from dials.array_family import flex
    from dials.util.command_line import Command

    # Get the flags
    flags = flex.reflection_table.flags

    # Create the algorithm
    algorithm = ImageSpaceProfileFitting(self._grid_size)

    # Add the specs
    for experiment, model in zip(self._experiments, self._profile_model):
      algorithm.add(Spec(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        model.delta_b(deg=False),
        model.delta_m(deg=False)))

    # Perform the integration
    num = reflections.get_flags(flags.dont_integrate).count(False)
    Command.start('Integrating %d reflections with profile fitting' % num)
    profiles = algorithm.execute(reflections)

    # Print the number integrated
    num = reflections.get_flags(flags.integrated_prf).count(True)
    Command.end('Integrated %d reflections with profile fitting' % num)

    # Output the reference profiles
    if self._debug:
      import cPickle as pickle
      filename = 'debug_%d.pickle' % job_id()
      print 'Saving debugging information to %s' % filename
      reference = [profiles.data(i) for i in range(len(profiles))]
      rprofiles = []
      for r in reflections:
        rprofiles.append(profiles.get(
          r['id'],
          r['panel'],
          r['s1'],
          r['xyzcal.mm'][2],
          r['bbox']))
      output = {
        'reflections' : reflections,
        'experiments' : self._experiments,
        'profile_model' : self._profile_model,
        'reference' : reference,
        'profiles' : rprofiles
      }
      with open(filename, 'wb') as outfile:
        pickle.dump(output, outfile, protocol=pickle.HIGHEST_PROTOCOL)



    # import numpy
    # numpy.set_printoptions(threshold='nan')
    # from dials.util import pprint
    # print pprint.profile3d(profiles.data(0))
    # for i in range(len(profiles)):
    #   data = profiles.data(i)
    #   vmax = flex.max(data)
    #   data = data.as_numpy_array()
    #   from matplotlib import pylab
    #   for k in range(data.shape[0]):
    #     pylab.subplot(1, data.shape[0],  k+1)
    #     pylab.imshow(data[k,:,:], interpolation='none', vmin=0, vmax=vmax)
    #   pylab.show()

    # Return the reflections
    return reflections
