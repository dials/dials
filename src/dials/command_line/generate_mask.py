"""
Mask images to remove unwanted pixels.

This program is used to generate mask to specify which pixels should be
considered "invalid" during spot finding and integration. It provides a few
options to create simple masks using the detector trusted range, or from
simple shapes or by setting different resolution ranges.

Examples::

  dials.generate_mask models.expt border=5

  dials.generate_mask models.expt \\
    untrusted.rectangle=50,100,50,100 \\
    untrusted.circle=200,200,100

  dials.generate_mask models.expt d_max=2.00
"""


from __future__ import annotations

import logging
import os.path
import pickle
from typing import List, Optional, Tuple

import libtbx.phil as phil
from dxtbx.format.image import ImageBool
from dxtbx.model.experiment_list import ExperimentList
from scitbx.array_family import flex

import dials.util
import dials.util.log
import dials.util.masking
from dials.util.options import ArgumentParser, flatten_experiments

Masks = List[Tuple[flex.bool, ...]]

log = logging.getLogger("dials.generate_mask")

phil_scope = phil.parse(
    """
    output {
        mask = pixels.mask
            .type = path
            .help = "Name of output mask file."
        experiments = None
            .type = path
            .help = "Name of output experiment list file.  If this is set, a copy of "
                    "the experiments, modified with the generated pixel masks, "
                    "will be saved to this location."
        log = 'dials.generate_masks.log'
            .type = str
            .help = "The log filename."
    }

    include scope dials.util.masking.phil_scope
    """,
    process_includes=True,
)


def generate_mask(
    experiments: ExperimentList,
    params: phil.scope_extract,
) -> Tuple[Masks, Optional[ExperimentList]]:
    """
    Generate a pixel mask for each imageset in an experiment list.

    Use the masking parameters :param:`params` and the experiments in the experiment
    list :param:`experiments` to define a pixel mask for each of the associated
    imagesets.  The masks are generated using :mod:`dials.util.masking`.

    The masks will be saved to disk at the location specified by
    :attr:`params.output.mask`.  If the experiment list contains more than one
    imageset, multiple mask files will be produced, with filenames differentiated by
    an appended number.  Optionally, if a path :attr:`params.output.experiments` is
    set, a modified copy of :param:`experiments` with the masks applied will be saved
    to that location.

    Args:
        experiments: An experiment list containing only one imageset.
        params: Masking parameters, having the structure defined in
            :data:`phil_scope`.

    Returns:
        A list of masks, one for each imageset.

        A copy of :param:`experiments` with the masks applied (optional,
        only returned if :attr:`params.output.experiments` is set).
    """
    imagesets = experiments.imagesets()
    masks = []

    # Create output mask filenames
    num_imagesets = len(imagesets)
    if num_imagesets == 1:
        filenames = [params.output.mask]
    else:
        # If there is more than one imageset, append a number to each output filename
        name, ext = os.path.splitext(params.output.mask)
        pad = len(str(num_imagesets))
        filenames = [
            "{name}_{num:0{pad}}{ext}".format(name=name, num=i + 1, pad=pad, ext=ext)
            for i in range(num_imagesets)
        ]

    for imageset, filename in zip(imagesets, filenames):
        mask = dials.util.masking.generate_mask(imageset, params)
        masks.append(mask)

        # Save the mask to file
        log.info("Writing mask to %s", filename)
        with open(filename, "wb") as fh:
            pickle.dump(mask, fh)

        if params.output.experiments:
            # Apply the mask to the imageset
            imageset.external_lookup.mask.data = ImageBool(mask)
            imageset.external_lookup.mask.filename = filename

    if params.output.experiments:
        # Save the experiment list
        log.info("Saving experiments to %s", params.output.experiments)
        experiments.as_file(params.output.experiments)
    else:
        experiments = None

    return masks, experiments


@dials.util.show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    """
    Parse command-line arguments, run the script.

    Uses the DIALS option parser to extract an experiment list and
    parameters, then passes these to :func:`generate_mask`.

    Args:
        phil: PHIL scope for option parser.
        args: Arguments to parse. If None, :data:`sys.argv[1:]` will be used.
    """
    # Create the parser
    usage = "usage: dials.generate_mask [options] models.expt"
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_experiments_from_images=True,
    )

    # Parse the command line arguments
    params, options = parser.parse_args(args=args, show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    # Configure logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Check number of args
    if len(experiments) == 0:
        parser.print_help()
        return

    # Run the script
    generate_mask(experiments, params)


if __name__ == "__main__":
    run()
