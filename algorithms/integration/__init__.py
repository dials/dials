from __future__ import absolute_import, division, print_function

from tabulate import tabulate

import boost.python
from dials_algorithms_integration_ext import *
from dials_algorithms_integration_kapton_ext import *
from dials.algorithms.integration.integrator_stills import *


class TimingInfo(object):
    """
    A class to contain timing info.
    """

    def __init__(self):
        self.read = 0
        self.extract = 0
        self.initialize = 0
        self.process = 0
        self.finalize = 0
        self.total = 0
        self.user = 0

    def __str__(self):
        """ Convert to string. """
        rows = [
            [description, "%.2f seconds" % value]
            for description, value in (
                ["Read time", self.read],
                ["Extract time", self.extract],
                ["Pre-process time", self.initialize],
                ["Process time", self.process],
                ["Post-process time", self.finalize],
                ["Total time", self.total],
                ["User time", self.user],
            )
            if value
        ]
        return tabulate(rows, tablefmt="rst")
