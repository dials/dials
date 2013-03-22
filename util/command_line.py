#
# command_line.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def parse_range_list_string(string):
    """Parse a string in the following ways:
    string: 1, 2, 3        -> [1, 2, 3]
    string: 1 - 6          -> [1, 2, 3, 4, 5, 6]
    string: 1 - 6, 7, 8, 9 -> [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    items = string.split(',')
    for i in range(len(items)):
        items[i] = items[i].split("-")
        if len(items[i]) == 1:
            items[i] = [int(items[i][0])]
        elif len(items[i]) == 2:
            items[i] = range(int(items[i][0]), int(items[i][1]) + 1)
        else:
            raise SyntaxError
    items = [item for sublist in items for item in sublist]
    return set(items)

def interactive_console(namespace):
    """ Enter an interactive console session. """
    try:
        from IPython import embed
        embed(user_ns = namespace)
    except ImportError:
        print "IPython not available"


class ProgressBarTimer:
    """ A simple timer for the progress bar. """

    def __init__(self):
        """ Init the progress bar timer. """
        from time import time
        self._last_time = time()
        self._last_perc = 0
        self._update_period = 1
        self._n_seconds_left = -1

    def update(self, percent):
        """ Update the timer. """
        from time import time

        # Get the current time diff between last time
        curr_time = time()
        diff_time = curr_time - self._last_time

        # Only update after certain period or at 100%
        if percent < 0: percent = 0
        if percent > 100: percent = 100
        if diff_time >= self._update_period or percent >= 100:

            # Check the difference in percentage and calculate
            # number of seconds remaining
            diff_perc = percent - self._last_perc
            if (diff_perc == 0):
                self._n_seconds_left = 0
            else:
                self._n_seconds_left = diff_time * (100 - percent) / diff_perc

        # Return number of seconds
        return self._n_seconds_left

class ProgressBar:
    """ A command line progress bar. """

    def __init__(self, spinner=True, bar=True, estimate_time=True, length=50):
        """ Init the progress bar parameters. """

        # Set the parameters
        self._spinner = spinner
        self._estimate_time = estimate_time
        self._bar = bar
        self._length = length

        self._timer = ProgressBarTimer()

        # Print 0 percent
        self.update(0)

    def update(self, fpercent):
        """ Update the progress bar with a percentage. """
        import sys
        from math import ceil

        # Get integer percentage
        percent = int(fpercent)
        if percent < 0: percent = 0
        if percent > 100: percent = 100

        # Add a percentage counter
        progress_str = '\r'
        progress_str += '{0: >3}%'.format(percent)

        # Add a spinner
        if self._spinner:
            progress_str += ' '
            progress_str += '[ {0} ]'.format('-\|/'[percent % 4])

        # Add a bar
        if self._bar:
            n_char = int(percent * self._length / 100)
            n_space = self._length - n_char
            progress_str += ' '
            progress_str += '[ {0}>{1} ]'.format('=' * n_char, ' ' * n_space)

        # Add a timer
        if self._estimate_time:
            n_seconds_left = self._timer.update(fpercent)
            if n_seconds_left < 0:
                n_seconds_left = '?'
            else:
                n_seconds_left = int(ceil(n_seconds_left))
            time_str = '{0: >5}'.format(n_seconds_left)
            progress_str += ' '
            progress_str += 'remaining: {0}s'.format(time_str)

        # Print progress string to stdout
        sys.stdout.write(progress_str)
        sys.stdout.flush()

    def finished(self):
        """ The progress bar is finished. """
        self.update(100)
        print ""
