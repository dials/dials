from __future__ import absolute_import, division, print_function

import collections

from dials.util import tabulate
from dials_algorithms_integration_ext import *  # noqa: F403; lgtm
from dials_algorithms_integration_kapton_ext import *  # noqa: F403; lgtm

__all__ = (  # noqa: F405
    "Corrections",
    "CorrectionsMulti",
    "lp_correction",
    "qe_correction",
    "Result",
)

Result = collections.namedtuple(
    "Result",
    "index, reflections, data, read_time, extract_time, process_time, total_time",
)
#        :param index: The processing job index
#        :param reflections: The processed reflections
#        :param data: Other processed data


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
        """Convert to string."""
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
        return tabulate(rows)
