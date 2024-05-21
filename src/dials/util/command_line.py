from __future__ import annotations

import argparse
import glob
import sys
import time

from dials.util import debug_console

interactive_console = debug_console


from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.phil import parse

from dials.array_family.flex import reflection_table


def standard_scope():
    """Define the standard scope as a function so every instance is unique,
    not global."""
    return parse(
        """
    input {
      experiments = None
        .type = path
        .multiple = True
      reflections = None
        .type = path
        .multiple = True
      }

      output {
      experiments = None
        .type = path
        .multiple = True
      reflections = None
        .type = path
        .multiple = True
      }
    """
    )


class OptionParser:
    """Thin wrapper around argparse and phil to get DIALS arguments
    parsed without actually loading all the data files."""

    @staticmethod
    def guess_input_file_type(filenames):
        """Guess if input filenames are experiment-like, reflection-like or
        param-like, or optionally image like. Return dictionary."""

        experiment_file_endings = ".expt".split()
        reflection_file_endings = ".refl".split()
        parameter_file_endings = ".params .phil".split()
        image_file_endings = ".h5 .nxs .cbf .img".split()
        image_compression = ["", ".gz", ".bz2"]

        experiment_files = [
            filename
            for filename in filenames
            if any(filename.endswith(exten) for exten in experiment_file_endings)
        ]
        reflection_files = [
            filename
            for filename in filenames
            if any(filename.endswith(exten) for exten in reflection_file_endings)
        ]
        parameter_files = [
            filename
            for filename in filenames
            if any(filename.endswith(exten) for exten in parameter_file_endings)
        ]
        image_files = [
            filename
            for filename in filenames
            if any(
                filename.endswith(exten + compression)
                for exten in image_file_endings
                for compression in image_compression
            )
        ]

        unknown = [
            filename
            for filename in filenames
            if not any(
                (
                    filename in experiment_files,
                    filename in reflection_files,
                    filename in parameter_files,
                    filename in image_files,
                )
            )
        ]

        return {
            "experiments": experiment_files,
            "reflections": reflection_files,
            "parameters": parameter_files,
            "images": image_files,
            "unknown": unknown,
        }

    def __init__(self, phil=None):
        self._input_experiments = None
        self._input_reflections = None

        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-c",
            "--show-config",
            action="store_true",
            dest="show_config",
            help="Show the configuration parameters.",
        )
        parser.add_argument(
            "-a",
            "--attributes-level",
            type=int,
            default=0,
            dest="attributes_level",
            help="Set the attributes level for showing configuration parameters",
        )
        parser.add_argument(
            "-e",
            "--expert-level",
            type=int,
            default=0,
            dest="expert_level",
            help="Set the expert level for showing configuration parameters",
        )

        parser.add_argument("phil", nargs="*")
        args = parser.parse_args()

        standard = standard_scope()

        if phil:
            standard.adopt_scope(phil)
        clai = standard.command_line_argument_interpreter()

        phil_params = [arg for arg in args.phil if "=" in arg]

        other = OptionParser.guess_input_file_type(
            [arg for arg in args.phil if "=" not in arg]
        )

        for expt in other["experiments"]:
            standard = standard.fetch(clai.process(f"input.experiments={expt}"))

        for refl in other["reflections"]:
            standard = standard.fetch(clai.process(f"input.reflections={refl}"))

        for image in other["images"]:
            standard = standard.fetch(clai.process(f"input.experiments={image}"))

        self._phil = standard.fetch(clai.process_and_fetch(phil_params))
        self._params = self._phil.extract()

        if args.phil == [] or hasattr(args, "show_config") and args.show_config:
            print(
                "Showing configuration parameters with:\n"
                "  attributes_level = %d\n"
                "  expert_level = %d\n" % (args.attributes_level, args.expert_level)
            )
            print(
                self._phil.as_str(
                    expert_level=args.expert_level,
                    attributes_level=args.attributes_level,
                )
            )
            sys.exit(0)

    def __repr__(self):
        return self._phil.format(python_object=self._params).as_str()

    def params(self):
        return self._params

    def input_experiments(self):
        return sum(map(glob.glob, self._params.input.experiments), [])

    def input_reflections(self):
        return sum(map(glob.glob, self._params.input.reflections), [])

    def output_experiments(self):
        return self._params.output.experiments

    def output_reflections(self):
        return self._params.output.reflections

    @staticmethod
    def read_experiments(experiment_filenames, check_format=True):
        inputs = OptionParser.guess_input_file_type(experiment_filenames)
        experiments = ExperimentListFactory.from_filenames(inputs["images"])
        for experimentlist in [
            ExperimentListFactory.from_json_file(filename, check_format=check_format)
            for filename in inputs["experiments"]
        ]:
            experiments.extend(experimentlist)
        return experiments

    @staticmethod
    def read_reflections(reflection_filenames):
        return [
            reflection_table.from_file(filename) for filename in reflection_filenames
        ]


class ProgressBarTimer:
    """A simple timer for the progress bar."""

    def __init__(self):
        """Init the progress bar timer."""
        self._last_time = time.time()
        self._last_perc = 0
        self._update_period = 0.5
        self._n_seconds_left = -1

    def get_elapsed_time(self):
        return time.time() - self._last_time

    def update(self, percent):
        """Update the timer."""
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
    """A command line progress bar."""

    def __init__(
        self,
        title=None,
        spinner=True,
        bar=True,
        estimate_time=True,
        indent=0,
        length=80,
    ):
        """Init the progress bar parameters."""

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
        """Update the progress bar with a percentage."""
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

        left_str += f"{percent: >3}%"

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
            right_str = " " + f"est: {n_seconds_left}s" + right_str

        # Add a bar
        if self._bar:
            bar_length = self._length - (len(left_str) + len(right_str)) - 5
            n_char = int(percent * bar_length / 100)
            n_space = bar_length - n_char
            left_str += " "
            left_str += f"[ {'=' * n_char}>{' ' * n_space} ]"

        # Append strings
        progress_str = left_str + right_str

        # Print progress string to stdout
        sys.stdout.write(progress_str)
        sys.stdout.flush()

    def finished(self, string=None):
        """The progress bar is finished."""
        if string:
            self._title = string
        else:
            string = ""
        """ Print the 'end of comand' string."""
        if self._estimate_time:
            # Get the time string
            time_string = f"{self._timer.get_elapsed_time():.2f}s"

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


class Command:
    """Class to nicely print out a command with timing info."""

    # Variables available in class methods
    indent = 0
    max_length = 80
    print_time = True

    @classmethod
    def start(cls, string):
        """Print the 'start command' string."""
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
        """Print the 'end of command' string."""
        # from termcolor import colored

        # Check if we want to print the time or not
        if cls.print_time:

            # Get the time string
            time_string = f"{time.time() - cls._start_time:.2f}s"

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
