from __future__ import division
from dials.framework import toplevel

class SpotFinderInterface(toplevel.Interface):

  def __init__(self, sweep, crystal, params):

    from dials.framework.plugin import Registry
    from dials.interfaces.standard_interfaces import SpotFinderThreshold
    from dials.interfaces.standard_interfaces import SpotFinderFilter

    # Save stuff internally
    self.sweep = sweep

    # Configure the thresholding algorithm
    threshold_algorithm = Registry.factory(SpotFinderThreshold)(
        params.spotfinder.threshold.algorithm,
        sweep, crystal, params)

    # Configure the filtering algorithm
    filters = []
    for name in params.spotfinder.filter.algorithm:
      filters.append(Registry.factory(SpotFinderFilter)(
          name, sweep, crystal, params))
    filter_algorithm = FilterRunner(filters)

    # Configure the spot finding algorithm
    self.spotfinder = SpotFinder(
        threshold_algorithm=threshold_algorithm,
        filter_algorithm=filter_algorithm,
        scan_range=params.spotfinder.scan_range)

  def prepare(self):
    self.prepared = True

  def process(self):
    self.result = self.spotfinder(self.sweep)
    self.processed = True

  def finish(self):
    self.finished = True


class SpotFinderFactory(object):

  @staticmethod
  def create(sweep, crystal, params):
    return SpotFinderInterface(sweep, crystal, params)
