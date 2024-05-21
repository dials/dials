from __future__ import annotations

import argparse
import copy
import itertools
import logging
import os
import pickle
import sys
import traceback
import warnings
from collections import defaultdict, namedtuple
from glob import glob

import libtbx.phil
from dxtbx.model import ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.util import get_url_scheme
from orderedset import OrderedSet

from dials.array_family import flex
from dials.util import Sorry
from dials.util.multi_dataset_handling import (
    renumber_table_id_columns,
    sort_tables_to_experiments_order,
)
from dials.util.phil import FilenameDataWrapper


class InvalidPhilError(ValueError):
    pass


logger = logging.getLogger(__name__)

geometry_phil_scope = libtbx.phil.parse(
    """
geometry
  .help = "Allow overrides of experimental geometry"
  .expert_level = 2
{
  include scope dxtbx.model.beam.beam_phil_scope
  include scope dxtbx.model.detector.detector_phil_scope
  include scope dxtbx.model.goniometer.goniometer_phil_scope
  include scope dxtbx.model.scan.scan_phil_scope

  convert_stills_to_sequences = False
    .type = bool
    .help = "When overriding the scan, convert stills into sequences"
    .short_caption = "Convert stills into sequences"

  convert_sequences_to_stills = False
    .type = bool
    .help = "When overriding the scan, convert sequences into stills"
    .short_caption = "Convert sequences into stills"
}
""",
    process_includes=True,
)


class PhilCommandParser:
    """A class to parse phil parameters from positional arguments"""

    def __init__(
        self,
        phil=None,
        read_experiments=False,
        read_reflections=False,
        read_experiments_from_images=False,
        check_format=True,
    ):
        """
        Initialise the parser.

        :param phil: The phil scope
        :param read_experiments: Try to read the experiments
        :param read_reflections: Try to read the reflections
        :param read_experiments_from_images: Try to read the experiments from images
        :param check_format: Check the format when reading images
        """
        from dials.util.phil import parse

        # Set the system phil scope
        if phil is None:
            self._system_phil = parse("")
        else:
            self._system_phil = copy.deepcopy(phil)

        # Set the flags
        self._read_experiments = read_experiments
        self._read_reflections = read_reflections
        self._read_experiments_from_images = read_experiments_from_images
        self._check_format = check_format

        # Adopt the input scope
        input_phil_scope = self._generate_input_scope()
        if input_phil_scope is not None:
            self.system_phil.adopt_scope(input_phil_scope)

        # Set the working phil scope
        self._phil = self.system_phil.fetch(source=parse(""))

    @property
    def phil(self):
        """
        Get the phil object

        :return: The phil scope
        """
        return self._phil

    @property
    def system_phil(self):
        """
        Get the system phil.

        :return: The system phil scope
        """
        return self._system_phil

    @property
    def diff_phil(self):
        """
        Get the diff phil.

        :return: The difference phil scope
        """
        return self.system_phil.fetch_diff(source=self.phil)

    def parse_args(
        self, args, verbose=False, return_unhandled=False, quick_parse=False
    ):
        """
        Parse the command line arguments.

        :param args: The input arguments
        :param verbose: Print verbose output
        :param return_unhandled: True/False also return unhandled arguments
        :param quick_parse: Return as fast as possible and without reading any data,
                            ignoring class constructor options.
        :return: The options and parameters and (optionally) unhandled arguments
        """
        from dxtbx.model.experiment_list import (
            BeamComparison,
            DetectorComparison,
            GoniometerComparison,
        )

        from dials.util.phil import parse

        # Parse the command line phil parameters
        user_phils = []
        experiment_files = []
        reflection_table_files = []
        unhandled = []
        interpreter = self.system_phil.command_line_argument_interpreter()

        def _is_a_phil_file(filename):
            return filename.endswith((".phil", ".param", ".params", ".eff", ".def"))

        def _is_an_experiment_file(filename):
            return filename.endswith(".expt")

        def _is_a_refl_file(filename: str) -> bool:
            return filename.endswith(".refl")

        for arg in args:
            if (
                _is_a_phil_file(arg)
                and os.path.isfile(arg)
                and os.path.getsize(arg) > 0
            ):
                try:
                    user_phils.append(parse(file_name=arg))
                except Exception:
                    if return_unhandled:
                        unhandled.append(arg)
                    else:
                        raise
            elif (
                _is_an_experiment_file(arg)
                and os.path.isfile(arg)
                and os.path.getsize(arg) > 0
            ):
                experiment_files.append(arg)

            elif (
                _is_a_refl_file(arg)
                and os.path.isfile(arg)
                and os.path.getsize(arg) > 0
            ):
                reflection_table_files.append(arg)

            # Treat "has a schema" as "looks like a URL (not phil)
            elif "=" in arg and not get_url_scheme(arg):
                try:
                    user_phils.append(interpreter.process_arg(arg=arg))
                except Exception:
                    if return_unhandled:
                        unhandled.append(arg)
                    else:
                        raise
            else:
                unhandled.append(arg)

        # Fetch the phil parameters
        self._phil, unused = self.system_phil.fetch(
            sources=user_phils, track_unused_definitions=True
        )

        # Print if bad definitions
        if len(unused) > 0:
            msg = [item.object.as_str().strip() for item in unused]
            msg = "\n".join(["  %s" % line for line in msg])
            raise RuntimeError(f"The following definitions were not recognised\n{msg}")

        # Extract the parameters
        try:
            params = self._phil.extract()
        except Exception as e:
            raise InvalidPhilError(e)

        # Stop at this point if quick_parse is set. A second pass may be needed.
        if quick_parse:
            return params, unhandled

        try:
            load_models = params.load_models
        except AttributeError:
            load_models = True

        # Add the cached arguments
        for obj in importer.experiments:
            params.input.experiments.append(obj)
        for obj in importer.reflections:
            params.input.reflections.append(obj)

        # Convert to phil
        self._phil = self.system_phil.format(python_object=params)

        return params, importer.unhandled

    def _generate_input_scope(self):
        """
        Generate the required input scope.

        :return: The input phil scope
        """
        from dials.util.phil import parse

        # Create the input scope
        require_input_scope = (
            self._read_experiments
            or self._read_reflections
            or self._read_experiments_from_images
        )
        if not require_input_scope:
            return None
        input_phil_scope = parse("input {}")
        main_scope = input_phil_scope.get_without_substitution("input")
        assert len(main_scope) == 1
        main_scope = main_scope[0]

        # Return the input scope
        return input_phil_scope


class ArgumentParserBase(argparse.ArgumentParser):
    """The base class for the option parser."""

    def __init__(self, config_options=False, sort_options=False, **kwargs):
        """
        Initialise the class.

        :param config_options: True/False show configuration options
        :param sort_options: True/False show argument sorting options
        """

        # Initialise the option parser
        super().__init__(add_help=False, **kwargs)

        # Add an option to show configuration parameters
        if config_options:
            self.add_argument(
                "-c",
                "--show-config",
                action="store_true",
                default=False,
                dest="show_config",
                help="Show the configuration parameters.",
            )
            self.add_argument(
                "-a",
                "--attributes-level",
                default=0,
                type=int,
                dest="attributes_level",
                help="Set the attributes level for showing configuration parameters",
            )
            self.add_argument(
                "-e",
                "--expert-level",
                type=int,
                default=0,
                dest="expert_level",
                help="Set the expert level for showing configuration parameters",
            )
            self.add_argument(
                "--export-autocomplete-hints",
                action="store_true",
                default=False,
                dest="export_autocomplete_hints",
                help=argparse.SUPPRESS,
            )

        # Add an option to sort
        if sort_options:
            self.add_argument(
                "-s",
                "--sort",
                action="store_true",
                dest="sort",
                default=False,
                help="Sort the arguments",
            )

        # Set a help parameter
        self.add_argument(
            "-h",
            "--help",
            action="count",
            default=0,
            dest="help",
            help="Show this help message and exit. Can be specified multiple times to increase verbosity.",
        )

        # Set a verbosity parameter
        self.add_argument(
            "-v",
            action="count",
            default=0,
            dest="verbose",
            help="Increase verbosity (can be specified multiple times)",
        )

        # Add an option for PHIL file to parse - PHIL files passed as
        # positional arguments are also read but this allows the user to
        # explicitly specify STDIN
        self.add_argument(
            "--phil",
            action="append",
            metavar="FILE",
            help="PHIL files to read. Pass '-' for STDIN. Can be specified multiple times, but duplicates ignored.",
        )

    def parse_known_args(self, args=None, quick_parse=False):
        """
        Parse the command line arguments and get system configuration.

        :param args: The arguments to parse.
        :returns: The options and phil parameters
        """

        # Parse the command line arguments, this will separate out
        # options (e.g. -o, --option) and positional arguments, in
        # which phil options will be included.
        options, args = super().parse_known_args(args=args)

        # Read any argument-specified PHIL file. Ignore duplicates.
        if options.phil:
            for philfile in OrderedSet(options.phil):
                # Should we read STDIN?
                if philfile == "-":
                    lines = sys.stdin.readlines()
                else:
                    # Otherwise, assume we've been given a path
                    with open(philfile) as phil_input:
                        lines = phil_input.readlines()
                # Add these to the unparsed argument list
                args.extend(l.strip() for l in lines)

        # Maybe sort the data
        if hasattr(options, "sort") and options.sort:
            args = sorted(args)

        # Return the parameters
        return options, args

    def format_epilog(self, formatter):
        """Don't do formatting on epilog."""
        if self.epilog is None:
            return ""
        return self.epilog


class ArgumentParser(ArgumentParserBase):
    """A class to parse command line options and get the system configuration.
    The class extends argparse.ArgumentParser to include the reading of phil
    parameters."""

    def __init__(
        self,
        phil=None,
        read_experiments=False,
        read_reflections=False,
        read_experiments_from_images=False,
        check_format=True,
        sort_options=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        **kwargs,
    ):
        """
        Initialise the class.

        :param phil: The phil scope
        :param read_experiments: Try to read the experiments
        :param read_reflections: Try to read the reflections
        :param read_experiments_from_images: Try to read the experiments from images
        :param check_format: Check the format when reading images
        :param sort_options: Show argument sorting options
        """

        # Create the phil parser
        self._phil_parser = PhilCommandParser(
            phil=phil,
            read_experiments=read_experiments,
            read_reflections=read_reflections,
            read_experiments_from_images=read_experiments_from_images,
            check_format=check_format,
        )

        # Initialise the option parser
        super().__init__(
            sort_options=sort_options,
            config_options=self.system_phil.as_str() != "",
            formatter_class=formatter_class,
            **kwargs,
        )

    def parse_args(
        self,
        args=None,
        show_diff_phil=False,
        return_unhandled=False,
        ignore_unhandled=False,
        quick_parse=False,
    ):
        """
        Parse the command line arguments and get system configuration.

        :param args: The input arguments
        :param show_diff_phil: True/False Print the diff phil
        :param return_unhandled: True/False return unhandled arguments
        :param ignore_unhandled: True/False ignore unhandled arguments
                                  if return_unhandled is False
        :param quick_parse: Return as fast as possible and without reading any data,
                            ignoring class constructor options
        :return: The options and phil parameters
        """

        # Parse the command line arguments, this will separate out
        # options (e.g. -o, --option) and positional arguments, in
        # which phil options will be included.
        options, args = super().parse_known_args(args=args, quick_parse=quick_parse)

        if options.help:
            self.print_help()

        # Show config
        if hasattr(options, "show_config") and options.show_config:
            show_config = True
            attributes_level = options.attributes_level
            expert_level = options.expert_level
        elif options.help:
            show_config = True
            attributes_level = options.verbose
            expert_level = options.help - 1
        else:
            show_config = False

        if show_config:
            print(
                "Showing configuration parameters with:\n"
                f"  attributes_level = {attributes_level}\n"
                f"  expert_level = {expert_level}\n"
            )
            print(
                self.phil.as_str(
                    expert_level=expert_level,
                    attributes_level=attributes_level,
                )
            )
            exit(0)

        if (
            hasattr(options, "export_autocomplete_hints")
            and options.export_autocomplete_hints
        ):
            self._export_autocomplete_hints()
            exit(0)

        # Parse the phil parameters
        try:
            params, args = self._phil_parser.parse_args(
                args,
                options.verbose > 0,
                return_unhandled=return_unhandled,
                quick_parse=quick_parse,
            )
        except InvalidPhilError as e:
            self.error(message=f"Invalid phil parameter: {e}")

        # Print the diff phil
        if show_diff_phil:
            diff_phil_str = self.diff_phil.as_str()
            if diff_phil_str != "":
                print("The following parameters have been modified:\n")
                print(diff_phil_str)

        # Return the parameters
        if return_unhandled:
            return params, options, args
        elif len(args) > 0 and not quick_parse:
            # Handle printing any messages to diagnose unhandled arguments
            msg = self._warn_about_unhandled_args(args, verbosity=options.verbose)

            if ignore_unhandled:
                print(msg)
            else:
                raise Sorry(msg)
        return params, options

    def _warn_about_unhandled_args(self, unhandled, verbosity=0):
        """
        Generate any messages about unhandled arguments.

        This separates errors by validation/non-validation related, and only
        gives the user validation information if there's no other reason (or
        asked for verbose output).

        :param unhandled: List of unhandled arguments
        :param verbosity: The output verbosity determined during parsing
        :returns:         A formatted information string
        """
        msg = []
        msg.append("Unable to handle the following arguments:")

        # If we have any detailed information about why any of the
        # arguments weren't processed, give this to the user
        for arg in [x for x in unhandled if x in self._phil_parser.handling_errors]:
            # Split the reasons for unhandling into validation, non-validation
            non_valid = [
                x for x in self._phil_parser.handling_errors[arg] if not x.validation
            ]
            valid = [x for x in self._phil_parser.handling_errors[arg] if x.validation]

            # If we have non-validation-related errors, these are more important
            for _, err in itertools.groupby(non_valid, key=lambda x: x.message):
                err = list(err)
                # Grouping the errors by message lets us avoid repeating messages
                if len(err) > 1:
                    msg.append(
                        '  "{}" failed repeatedly during processing:\n{}\n'.format(
                            arg, "    " + err[0].message
                        )
                    )
                elif isinstance(err[0].exception, Sorry):
                    msg.append(
                        '  "{}" failed during {} processing:\n    {}\n'.format(
                            arg, err[0].type, err[0].message
                        )
                    )
                else:
                    msg.append(
                        '  "{}" failed during {} processing:\n{}\n'.format(
                            arg,
                            err[0].type,
                            "\n".join(
                                "    " + x for x in err[0].traceback.splitlines()
                            ),
                        )
                    )

            # Otherwise (or if asked for verbosity), list the validation errors
            if valid and (not non_valid or verbosity > 1):
                msg.append(
                    "  {} did not appear to conform to any{} expected format:".format(
                        arg, " other" if non_valid else ""
                    )
                )
                slen = max(len(x.type) for x in valid)
                for err in valid:
                    msg.append(f"    - {f'{err.type}:'.ljust(slen + 1)} {err.message}")
        # The others
        for arg in [x for x in unhandled if x not in self._phil_parser.handling_errors]:
            msg.append("  " + str(arg))
        return "\n".join(msg)

    @property
    def phil(self):
        """
        Get the phil object

        :returns: The phil scope
        """
        return self._phil_parser.phil

    @property
    def system_phil(self):
        """
        Get the system phil.

        :returns: The system phil scope
        """
        return self._phil_parser.system_phil

    @property
    def diff_phil(self):
        """
        Get the diff phil.

        :returns: The diff phil scope
        """
        return self._phil_parser.diff_phil

    def _strip_rst_markup(self, text):
        """
        Strip rst markup

        :param text: The text to strip
        :return: The stripped text
        """
        return text.replace("::", ":")

    def format_help(self):
        """
        Format the help string

        :param formatter: The formatter to use
        :return: The formatted help text
        """
        result = super().format_help()
        return self._strip_rst_markup(result)

    def _export_autocomplete_hints(self):
        # complete list of all parameters
        parameter_list = []
        # short name -> full name expansion for unique names
        parameter_expansion_list = {}
        # full name -> list of flags for choice parameters
        parameter_choice_list = {}

        for d in self.phil.all_definitions():
            # Create complete list of all parameters
            parameter_list.append(d.path)

            # Create expansion (alias) list for unique names that are not expert commands
            if (
                d.object.name not in parameter_expansion_list
                and d.object.expert_level is None
                and d.parent.expert_level is None
            ):
                parameter_expansion_list[d.object.name] = d.path
            else:
                parameter_expansion_list[d.object.name] = None

            # Extract parameter choice lists
            if d.object.type.phil_type == "choice":
                parameter_choice_list[d.path] = [
                    w[1:] if w.startswith("*") else w
                    for w in (str(x) for x in d.object.words)
                ]
            elif d.object.type.phil_type == "bool":
                parameter_choice_list[d.path] = ["true", "false"]

        def construct_completion_tree(paths):
            """Construct a tree of parameters, grouped by common prefixes"""

            # Split parameter paths at '.' character
            paths = [p.split(".", 1) for p in paths]

            # Identify all names that are directly on this level
            # or represent parameter groups with a common prefix
            top_elements = {f"{x[0]}{'=' if len(x) == 1 else '.'}" for x in paths}

            # Partition all names that are further down the tree by their prefix
            subpaths = {}
            for p in paths:
                if len(p) > 1:
                    if p[0] not in subpaths:
                        subpaths[p[0]] = []
                    subpaths[p[0]].append(p[1])

            # If there are prefixes with only one name beneath them, put them on the top level
            for s in list(subpaths.keys()):
                if len(subpaths[s]) == 1:
                    top_elements.remove(f"{s}.")
                    top_elements.add(f"{s}.{subpaths[s][0]}=")
                    del subpaths[s]

            result = {"": top_elements}
            # Revursively process each group
            for n, x in subpaths.items():
                result[n] = construct_completion_tree(x)

            return result

        print("function _dials_autocomplete_flags ()")
        print("{")
        print(' case "$1" in')
        for p in parameter_choice_list:
            print(f"\n  {p})")
            print(
                '   _dials_autocomplete_values="%s";;'
                % " ".join(parameter_choice_list[p])
            )
        print("\n  *)")
        print('    _dials_autocomplete_values="";;')
        print(" esac")
        print("}")

        print("function _dials_autocomplete_expansion ()")
        print("{")
        print(' case "$1" in')
        for p, exp in parameter_expansion_list.items():
            if exp is not None:
                print(f"\n  {p}=)")
                print(f'   _dials_autocomplete_values="{exp}=";;')
        print("\n  *)")
        print('    _dials_autocomplete_values="";;')
        print(" esac")
        print("}")

        tree = construct_completion_tree(parameter_list)

        def _tree_to_bash(prefix, tree):
            for subkey in tree:
                if subkey != "":
                    _tree_to_bash(prefix + subkey + ".", tree[subkey])
                    print(f"\n  {prefix + subkey + '.'}*)")
                    print(
                        '    _dials_autocomplete_values="%s";;'
                        % " ".join(
                            sorted(
                                [prefix + subkey + "." + x for x in tree[subkey][""]]
                            )
                        )
                    )

        print("function _dials_autocomplete_hints ()")
        print("{")
        print(' case "$1" in')
        _tree_to_bash("", tree)

        toplevelset = tree[""] | {
            p + "=" for p, exp in parameter_expansion_list.items() if exp is not None
        }

        print("\n  *)")
        print('    _dials_autocomplete_values="%s";;' % " ".join(sorted(toplevelset)))
        print(" esac")
        print("}")


class OptionParser(ArgumentParser):
    def __init__(self, *args, **kwargs):
        # Backwards compatibility 2021-11-10; 2022-04-06
        # Remove after Dec 2022
        warnings.warn(
            "OptionParser is deprecated, use ArgumentParser instead",
            UserWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)


def flatten_reflections(filename_object_list):
    """
    Flatten a list of reflections tables

    A check is also made for the 'id' values in the reflection tables, which are
    renumbered from 0..n-1 to avoid clashes. The experiment_identifiers dict is
    also updated if present in the input tables.

    :param filename_object_list: The parameter item
    :return: The flattened reflection table
    """
    tables = [o.data for o in filename_object_list]
    if len(tables) > 1:
        tables = renumber_table_id_columns(tables)
    return tables


def flatten_experiments(filename_object_list):
    """
    Flatten a list of experiment lists

    :param filename_object_list: The parameter item
    :return: The flattened experiment lists
    """

    result = ExperimentList()
    for o in filename_object_list:
        result.extend(o.data)
    return result


def reflections_and_experiments_from_files(
    reflection_file_object_list, experiment_file_object_list
):
    """Extract reflection tables and an experiment list from the files.
    If experiment identifiers are set, the order of the reflection tables is
    changed to match the order of experiments.
    """
    tables = flatten_reflections(reflection_file_object_list)

    experiments = flatten_experiments(experiment_file_object_list)

    if tables and experiments:
        tables = sort_tables_to_experiments_order(tables, experiments)

    return tables, experiments
