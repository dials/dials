"""
This module defines a general cross validation function that can be used
with any valid CrossValidator. To use in a command line program, the
phil_scope should be included. Provided here is a description of the options,
using the example of running cross_validation in dials.scale

General description:
cross_validate runs the script defined by the CrossValidator, running each
option in turn, using a free set to score the model - the results are printed
in a table and the model with the lowest free set rmsd is indicated. For each
option, the analysis will be repeated nfolds times, with a different free set
chosen each time, and the final rmsds averaged. For full k-fold cross-validation,
nfolds should be set to 100/free_set_percentage, which would be nfolds=10 for
the default free_set_percentage=10.0.

Two different modes are currently supported, controlled by cross_validation_mode=;
1) cross_validation_mode=single
   dials.scale is run nfolds times for the user specified dials.scale options
2) cross_validation_mode=multi
   optimise a dials.scale parameter, specified by parameter= .
   parameter_values must also be specified as a string of space separated values,
   unless the dials.scale parameter type is bool.

Therefore one must choose:
  cross_validation_mode= (single or multi)
  parameter=             (a supported command line option of the script run by
                          CrossValidator, optional if cross_validation_mode=single)
  parameter_values=      (values to test, only optional if parameter= selectes a
                          boolean command-line parameter)

For example
cross_validation_mode=multi parameter=absorption_term
cross_validation_mode=multi parameter=decay_interval parameter_values="5.0 10.0 15.0"
cross_validation_mode=multi parameter=model parameter_values="array physical"
"""

from __future__ import absolute_import, division, print_function

import logging
import itertools
import time

from libtbx import phil
import six

logger = logging.getLogger("dials")

phil_scope = phil.parse(
    """
  cross_validation {
    cross_validation_mode = multi single
      .type = choice
      .help = "Choose the cross validation running mode, for a full description"
              "see the module docstring. Choice is used for testing a parameter"
              "that can only have discreet values (a choice or bool phil parameter)."
              "Variable is used for testing a parameter that can have a float or"
              "int value (that is also not a 'choice' type). Single just performs"
              "cross validation on one parameter configuration."
      .expert_level = 2
    parameter = None
      .type = str
      .help = "Optimise a command-line parameter. parameter_values must also be"
              "specified, unless the parameter is a True/False option."
      .expert_level = 2
    parameter_values = None
      .type = strings
      .help = "Parameter values to compare, entered as a string of space"
              "separated values."
      .expert_level = 2
    nfolds = 1
      .type = int(value_min=1)
      .help = "Number of cross-validation folds to perform. If nfolds > 1, the"
              "minimisation for each option is repeated nfolds times, with an"
              "incremental offset for the free set. The max number of folds"
              "allowed is 1/free_set_percentage; if set greater than this then"
              "the repetition will finish afer 1/free_set_percentage folds."
      .expert_level = 2
  }
"""
)


def cross_validate(params, cross_validator):
    """Run cross validation script."""

    start_time = time.time()
    free_set_percentage = cross_validator.get_free_set_percentage(params)
    options_dict = {}

    if params.cross_validation.cross_validation_mode == "single":
        # just run the setup nfolds times
        cross_validator.create_results_dict(n_options=1)
        for n in range(params.cross_validation.nfolds):
            if n < 100.0 / free_set_percentage:
                params = cross_validator.set_free_set_offset(params, n)
                cross_validator.run_script(params, config_no=0)

    elif params.cross_validation.cross_validation_mode == "multi":
        # run each option nfolds times
        if params.cross_validation.parameter is None:
            raise ValueError(
                "parameter= must be set to specify what command line option should be optimised"
            )

        choice = params.cross_validation.parameter
        # #TODO extract allowed values to allow checking of user input

        # inspect the phil scope to see what the parameter type is e.g bool, int
        typ = cross_validator.get_parameter_type(choice)

        if typ == "bool" and not params.cross_validation.parameter_values:
            # values not specified, implied that should test both True and False
            options_dict[choice] = [True, False]
        else:
            if not params.cross_validation.parameter_values:
                raise ValueError(
                    "parameter_values= must be set to specify what options should be tested"
                )
            options_dict[choice] = []
            if typ == "bool":
                if (
                    "true" in params.cross_validation.parameter_values
                    or "True" in params.cross_validation.parameter_values
                ):
                    options_dict[choice].append(True)
                if (
                    "false" in params.cross_validation.parameter_values
                    or "False" in params.cross_validation.parameter_values
                ):
                    options_dict[choice].append(False)
            elif typ == "choice":
                for option in params.cross_validation.parameter_values:
                    options_dict[choice].append(option)
            elif typ == "int":
                for value in params.cross_validation.parameter_values:
                    options_dict[choice].append(int(value))
            elif typ == "float":
                for value in params.cross_validation.parameter_values:
                    options_dict[choice].append(float(value))
            else:
                raise ValueError("Error in interpreting parameter and parameter_values")

        # this code below should work for more than one parameter to be optimised,
        # but one cannot specify this yet from the command line
        keys, values = zip(*options_dict.items())

        cross_validator.create_results_dict(len(values[0]))
        cross_validator.set_results_dict_configuration(keys, values)

        for i, v in enumerate(itertools.product(*values)):
            e = dict(zip(keys, v))
            for k, val in six.iteritems(e):
                params = cross_validator.set_parameter(params, k, val)
            for n in range(params.cross_validation.nfolds):
                if n < 100.0 / free_set_percentage:
                    params = cross_validator.set_free_set_offset(params, n)
                    cross_validator.run_script(params, config_no=i)

    else:
        raise ValueError("Error in interpreting mode and options.")

    st = cross_validator.interpret_results()
    logger.info("Summary of the cross validation analysis: \n %s", st.format())

    finish_time = time.time()
    logger.info(
        "\nCross-validation finished.\nTotal time taken: {:.4f}s ".format(
            finish_time - start_time
        )
    )
    logger.info("\n" + "=" * 80 + "\n")
