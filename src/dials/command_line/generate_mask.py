"""
Mask images to remove unwanted pixels.

This program is used to generate mask to specify which pixels should be
considered "invalid" during spot finding and integration. It provides a few
options to create simple masks using the detector trusted range, or from
simple shapes or by setting different resolution ranges.

Masks can also be combined by including them as arguments.

Examples::

  dials.generate_mask models.expt border=5

  dials.generate_mask models.expt \\
    untrusted.rectangle=50,100,50,100 \\
    untrusted.circle=200,200,100

  dials.generate_mask models.expt d_max=2.00

  dials.generate_mask models.expt d_max=2.00 existing.mask

  dials.generate_mask backstop.mask shadow.mask
"""

from __future__ import annotations

import logging
import os.path
import pickle

import libtbx.phil as phil
from dxtbx.format.image import ImageBool
from dxtbx.model.experiment_list import ExperimentList
from scitbx.array_family import flex

import dials.util
import dials.util.log
import dials.util.masking
from dials.util.options import ArgumentParser, flatten_experiments

Masks = list[tuple[flex.bool, ...]]

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
        log = 'dials.generate_mask.log'
            .type = str
            .help = "The log filename."
    }

    merge_imageset_mask = False
        .type = bool
        .help = "If True, merge pixel masks defined in the imageset, such as one specified "
                "during dials.import and one provided by the dxtbx class."

    include scope dials.util.masking.phil_scope
    """,
    process_includes=True,
)


def generate_mask(
    experiments: ExperimentList,
    params: phil.scope_extract,
    existing_masks: None | Masks = None,
) -> tuple[Masks, ExperimentList | None]:
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
        existing_masks: list of masks to combine with the mask being
            generated

    Returns:
        A list of masks, one for each imageset.

        A copy of :param:`experiments` with the masks applied (optional,
        only returned if :attr:`params.output.experiments` is set).
    """

    if existing_masks:
        existing_mask = list(existing_masks[0])
        for mask in existing_masks[1:]:
            for panel_idx in range(len(existing_mask)):
                existing_mask[panel_idx] &= mask[panel_idx]
        existing_mask = tuple(existing_mask)

    # Check if only combining masks
    if not experiments and existing_masks:
        # Save the mask to file
        log.info("Writing mask to %s", params.output.mask)
        with open(params.output.mask, "wb") as fh:
            pickle.dump(existing_mask, fh)
        return

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
        mask = list(mask)

        if params.merge_imageset_mask:
            imageset_mask = imageset.get_mask(0)
            for m1, m2 in zip(mask, imageset_mask):
                m1 &= m2
        if existing_masks:
            for m1, m2 in zip(mask, existing_mask):
                m1 &= m2

        mask = tuple(mask)
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
def run(args: list[str] = None, phil: phil.scope = phil_scope) -> None:
    """
    Parse command-line arguments, run the script.

    Uses the DIALS option parser to extract an experiment list and
    parameters, then passes these to :func:`generate_mask`.

    Args:
        phil: PHIL scope for option parser.
        args: Arguments to parse. If None, :data:`sys.argv[1:]` will be used.
    """
    # Create the parser
    usage = "usage: dials.generate_mask [options] [models.expt] [mask, mask, ...]"
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_experiments_from_images=True,
    )

    # Parse the command line arguments
    params, options, unhandled = parser.parse_args(
        args=args, show_diff_phil=True, return_unhandled=True
    )
    experiments = flatten_experiments(params.input.experiments)

    # Configure logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Read in any starting masks
    remains_unhandled = []
    existing_masks = []
    for arg in unhandled:
        if os.path.isfile(arg):
            try:
                with open(arg, "rb") as fh:
                    mask = pickle.load(fh, encoding="bytes")
            except Exception:
                remains_unhandled.append(arg)
            else:
                if (
                    isinstance(mask, tuple)
                    and mask
                    and all(type(m) is flex.bool for m in mask)
                ):
                    existing_masks.append(mask)
                else:
                    print("Invalid mask file:", arg)
                    remains_unhandled.append(arg)
    if remains_unhandled:
        print(
            "Couldn't recognize the following arguments:", ", ".join(remains_unhandled)
        )
        parser.print_help()
        return

    # Check number of args
    if len(experiments) == 0 and len(existing_masks) == 0:
        parser.print_help()
        return

    if not all(
        m[0].focus() == p.focus() for m in zip(*existing_masks) for p in list(m)
    ):
        print("Not all input masks are of the same shape")
        parser.print_help()
        return

    # Run the script
    generate_mask(experiments, params, existing_masks=existing_masks)


if __name__ == "__main__":
    run()
