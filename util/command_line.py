#
# command_line.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import time
import sys

from dials.util import debug_console


def parse_range_list_string(string):
    """Parse a string in the following ways:
    string: 1, 2, 3        -> [1, 2, 3]
    string: 1 - 6          -> [1, 2, 3, 4, 5, 6]
    string: 1 - 6, 7, 8, 9 -> [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    items = string.split(",")
    for i in range(len(items)):
        items[i] = items[i].split("-")
        if len(items[i]) == 1:
            items[i] = [int(items[i][0])]
        elif len(items[i]) == 2:
            items[i] = list(range(int(items[i][0]), int(items[i][1]) + 1))
        else:
            raise SyntaxError
    items = [item for sublist in items for item in sublist]
    return set(items)


interactive_console = debug_console


class ProgressBarTimer:
    """ A simple timer for the progress bar. """

    def __init__(self):
        """ Init the progress bar timer. """
        self._last_time = time.time()
        self._last_perc = 0
        self._update_period = 0.5
        self._n_seconds_left = -1

    def get_elapsed_time(self):
        return time.time() - self._last_time

    def update(self, percent):
        """ Update the timer. """
        # Get the current time diff between last time
        curr_time = time.time()
        diff_time = curr_time - self._last_time

        # Only update after certain period or at 100%
        if percent < 0:
            percent = 0
        if percent > 100:
            percent = 100
        if diff_time >= self._update_period or percent >= 100:

            # Check the difference in percentage and calculate
            # number of seconds remaining
            diff_perc = percent - self._last_perc
            if diff_perc == 0:
                self._n_seconds_left = 0
            else:
                self._n_seconds_left = diff_time * (100 - percent) / diff_perc

        # Return number of seconds
        return self._n_seconds_left


class ProgressBar:
    """ A command line progress bar. """

    def __init__(
        self,
        title=None,
        spinner=True,
        bar=True,
        estimate_time=True,
        indent=0,
        length=80,
    ):
        """ Init the progress bar parameters. """

        # Set the parameters
        self._title = title
        self._indent = indent
        self._spinner = spinner
        self._estimate_time = estimate_time
        self._bar = bar
        self._length = length

        self._timer = ProgressBarTimer()

        # Print 0 percent
        self.update(0)

    def update(self, fpercent):
        """ Update the progress bar with a percentage. """
        from math import ceil

        # do not update if not a tty
        if not sys.stdout.isatty():
            return

        # Get integer percentage
        percent = int(fpercent)
        if percent < 0:
            percent = 0
        if percent > 100:
            percent = 100

        # Add a percentage counter
        right_str = ""
        left_str = ""
        if sys.stdout.isatty():
            left_str = "\r"
        left_str += " " * self._indent

        # Add a title if given
        if self._title:
            left_str += self._title + ": "

        left_str += "{: >3}%".format(percent)

        # Add a spinner
        if self._spinner:
            left_str += " "
            left_str += "[ {} ]".format(r"-\|/"[percent % 4])

        # Add a timer
        if self._estimate_time:
            n_seconds_left = self._timer.update(fpercent)
            if n_seconds_left < 0:
                n_seconds_left = "?"
            else:
                n_seconds_left = int(ceil(n_seconds_left))
            right_str = " " + "est: {}s".format(n_seconds_left) + right_str

        # Add a bar
        if self._bar:
            bar_length = self._length - (len(left_str) + len(right_str)) - 5
            n_char = int(percent * bar_length / 100)
            n_space = bar_length - n_char
            left_str += " "
            left_str += "[ {0}>{1} ]".format("=" * n_char, " " * n_space)

        # Append strings
        progress_str = left_str + right_str

        # Print progress string to stdout
        sys.stdout.write(progress_str)
        sys.stdout.flush()

    def finished(self, string=None):
        """ The progress bar is finished. """
        if string:
            self._title = string
        else:
            string = ""

        """ Print the 'end of comand' string."""
        if self._estimate_time:
            # Get the time string
            time_string = "{:.2f}s".format(self._timer.get_elapsed_time())

            # Truncate the string
            max_length = self._length - self._indent - len(time_string) - 1
            string = string[:max_length]

            # Add an indent and a load of dots and then the time string
            dot_length = 1 + max_length - len(string)
            string = (" " * self._indent) + string
            string = string + "." * (dot_length)
            string = string + time_string

        else:

            # Truncate the string
            max_length = self._length - self._indent
            string = string[:max_length]

            # Add a load of dots
            dot_length = max_length - len(string)
            string = (" " * self._indent) + string
            string = string + "." * (dot_length)

        # Write the string to stdout
        if sys.stdout.isatty():
            string = "\r" + string + "\n"
        else:
            string = string + "\n"
        sys.stdout.write(string)
        sys.stdout.flush()


class Command(object):
    """Class to nicely print out a command with timing info."""

    # Variables available in class methods
    indent = 0
    max_length = 80
    print_time = True

    @classmethod
    def start(cls, string):
        """ Print the 'start command' string."""
        # from termcolor import colored

        # Get the command start time
        cls._start_time = time.time()

        # do not output if not a tty
        if not sys.stdout.isatty():
            return

        # Truncate the string to the maximum length
        max_length = cls.max_length - cls.indent - 3
        string = string[:max_length]
        string = (" " * cls.indent) + string + "..."

        # Write the string to stdout
        sys.stdout.write(string)
        sys.stdout.flush()

    @classmethod
    def end(cls, string):
        """ Print the 'end of command' string."""
        # from termcolor import colored

        # Check if we want to print the time or not
        if cls.print_time:

            # Get the time string
            time_string = "{:.2f}s".format(time.time() - cls._start_time)

            # Truncate the string
            max_length = cls.max_length - cls.indent - len(time_string) - 1
            string = string[:max_length]

            # Add an indent and a load of dots and then the time string
            dot_length = 1 + max_length - len(string)
            string = (" " * cls.indent) + string
            string = string + "." * (dot_length)
            string = string + time_string

        else:

            # Truncate the string
            max_length = cls.max_length - cls.indent
            string = string[:max_length]

            # Add a load of dots
            dot_length = max_length - len(string)
            string = (" " * cls.indent) + string
            string = string + "." * (dot_length)

        # Write the string to stdout
        if sys.stdout.isatty():
            string = "\r" + string + "\n"
        else:
            string = string + "\n"
        sys.stdout.write(string)
        sys.stdout.flush()


try:
    import termcolor
except ImportError:
    termcolor = None


def coloured(text, *args, **kwargs):
    if not sys.stdout.isatty() or termcolor is None:
        return text
    return termcolor.colored(text, *args, **kwargs)


def heading(text):
    return coloured(text, attrs=["bold"])


if __name__ == "__main__":
    p = ProgressBar()

    for j in range(100):
        p.update(j)
        time.sleep(0.05)

    p.finished()

    Command.start("Starting to do a command")
    time.sleep(1)
    Command.end("Ending the command")
