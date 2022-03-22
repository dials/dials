# LIBTBX_SET_DISPATCHER_NAME dials.background
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1


from __future__ import annotations

import logging
import math

import iotbx.phil
from libtbx.phil import parse
from scitbx import matrix

import dials.util.masking
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.background")

help_message = """

Examples::

  dials.background image_*.cbf

  dials.background imported.expt
"""

phil_scope = iotbx.phil.parse(
    """\
n_bins = 100
  .type = int
images = None
  .type = ints
  .help = "Images on which to perform the analysis (otherwise use all images)"
corrected = False
  .type = bool
  .help = "Use corrected data (i.e after applying pedestal and gain) in analysis"

masking {
  include scope dials.util.masking.phil_scope
}

output {
    plot = None
      .type = path
      .help = "Save background plot to file"
    size_inches = None
      .type = floats(value_min=0, size=2)
    log = dials.background.log
      .type = str
}

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.background [options] image_*.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    # Ensure we have either a data block or an experiment list
    experiments = flatten_experiments(params.input.experiments)
    imagesets = experiments.imagesets()

    # Configure the logging
    log.config(logfile=params.output.log)
    logger.info(dials_version())

    if params.output.plot:
        import matplotlib

        matplotlib.use("agg")

        from matplotlib import pyplot

        from dials.util.matplotlib_utils import resolution_formatter

        fig = pyplot.figure(figsize=params.output.size_inches)
        ax = fig.add_subplot(111)

    for i_imgset, imageset in enumerate(imagesets):
        first, last = imageset.get_scan().get_image_range()
        images = range(first, last + 1)

        if params.images:
            if min(params.images) < first or max(params.images) > last:
                raise Sorry("image outside of scan range")
            images = params.images

        d_spacings = []
        intensities = []
        sigmas = []

        for indx in images:
            logger.info(f"For imageset {i_imgset} image {indx}:")
            d, I, sig = background(
                imageset,
                indx - first,  # indices passed to imageset.get_raw_data start from zero
                n_bins=params.n_bins,
                corrected=params.corrected,
                mask_params=params.masking,
                show_summary=True,
            )

            msg = [f"{'d':>8} {'I':>8} {'sig':>8}"]
            for j in range(len(I)):
                msg.append(f"{d[j]:8.3f} {I[j]:8.3f} {sig[j]:8.3f}")
            logger.info("\n".join(msg))

            d_spacings.append(d)
            intensities.append(I)
            sigmas.append(sig)

        if params.output.plot:
            ax.set_xlabel(r"resolution ($\AA$)")
            ax.set_ylabel(r"$\langle I_b \rangle$")
            for indx, d, I, sig in zip(images, d_spacings, intensities, sigmas):
                filenames = imageset.reader().paths()
                if len(imagesets) > 1:
                    label = (
                        f"{filenames[indx - first]}"
                        if len(filenames) > 1
                        else f"{filenames[0]} image {indx}"
                    )
                else:
                    label = f"image {indx}" if len(images) > 1 else ""
                ds2 = 1 / flex.pow2(d)
                ax.plot(ds2, I, label=label)
            ax.xaxis.set_major_formatter(resolution_formatter)

    if params.output.plot:
        try:
            if len(imagesets) > 1 or len(images) > 1:
                # Plot a legend if there are fewer lines than the number of colours
                # in the colour cycle
                if len(ax.lines) <= len(
                    pyplot.rcParams["axes.prop_cycle"].by_key()["color"]
                ):
                    pyplot.gca().legend()
            pyplot.savefig(params.output.plot)
        except ValueError:
            raise Sorry(f"Unable to save plot to {params.output.plot}")


def background(
    imageset, indx, n_bins, corrected=False, mask_params=None, show_summary=False
):
    if mask_params is None:
        # Default mask params for trusted range
        mask_params = phil_scope.fetch(parse("")).extract().masking

    detector = imageset.get_detector()
    beam = imageset.get_beam()

    # Only working with single panel detector for now
    assert len(detector) == 1
    panel = detector[0]
    imageset_mask = imageset.get_mask(indx)[0]
    mask = dials.util.masking.generate_mask(imageset, mask_params)[0]
    mask = imageset_mask & mask

    n = matrix.col(panel.get_normal()).normalize()
    b = matrix.col(beam.get_s0()).normalize()
    wavelength = beam.get_wavelength()

    if math.fabs(b.dot(n)) < 0.95:
        raise Sorry("Detector not perpendicular to beam")

    # Use corrected data to determine signal and background regions
    corrected_data = imageset.get_corrected_data(indx)
    assert len(corrected_data) == 1
    corrected_data = corrected_data[0].as_double()

    # Use choice of raw or corrected data to evaluate the background values
    if corrected:
        data = corrected_data
    else:
        data = imageset.get_raw_data(indx)[0].as_double()

    spot_params = spot_phil.fetch(source=parse("")).extract()
    threshold_function = SpotFinderFactory.configure_threshold(spot_params)
    peak_pixels = threshold_function.compute_threshold(corrected_data, mask)
    signal = data.select(peak_pixels.iselection())
    background_pixels = mask & ~peak_pixels
    background = data.select(background_pixels.iselection())

    # print some summary information
    if show_summary:
        logger.info(f"Mean background: {flex.sum(background) / background.size():.3f}")
        if len(signal) > 0:
            logger.info(
                f"Max/total signal pixels: {flex.max(signal):.0f} / {flex.sum(signal):.0f}"
            )
        else:
            logger.info("No signal pixels on this image")
        logger.info(
            "Peak/background/masked pixels: %d / %d / %d"
            % (peak_pixels.count(True), background.size(), mask.count(False))
        )

    # compute histogram of two-theta values, then same weighted
    # by pixel values, finally divide latter by former to get
    # the radial profile out, need to set the number of bins
    # sensibly; inspired by method in PyFAI

    two_theta_array = panel.get_two_theta_array(beam.get_s0())
    two_theta_array = two_theta_array.as_1d().select(background_pixels.iselection())

    # Use flex.weighted_histogram
    h0 = flex.weighted_histogram(two_theta_array, n_slots=n_bins)
    h1 = flex.weighted_histogram(two_theta_array, background, n_slots=n_bins)
    h2 = flex.weighted_histogram(
        two_theta_array, background * background, n_slots=n_bins
    )

    d0 = h0.slots()
    d1 = h1.slots()
    d2 = h2.slots()

    I = d1 / d0
    I2 = d2 / d0
    sig = flex.sqrt(I2 - flex.pow2(I))

    tt = h0.slot_centers()
    d_spacings = wavelength / (2.0 * flex.sin(0.5 * tt))

    return d_spacings, I, sig


if __name__ == "__main__":
    run()
