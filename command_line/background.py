# LIBTBX_SET_DISPATCHER_NAME dials.background
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import absolute_import, division, print_function

import math

import iotbx.phil
from libtbx.phil import parse
from scitbx import matrix

import dials.util.masking
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.array_family import flex
from dials.util import Sorry, show_mail_handle_errors
from dials.util.options import OptionParser, flatten_experiments

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
}

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.background [options] image_*.cbf"

    parser = OptionParser(
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

    if params.output.plot:
        import matplotlib

        matplotlib.use("agg")

        import matplotlib.ticker as mticker
        from matplotlib import pyplot

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
            print(f"For imageset {i_imgset} image {indx}:")
            d, I, sig = background(
                imageset,
                indx - first,  # indices passed to imageset.get_raw_data start from zero
                n_bins=params.n_bins,
                corrected=params.corrected,
                mask_params=params.masking,
            )

            print("%8s %8s %8s" % ("d", "I", "sig"))
            for j in range(len(I)):
                print("%8.3f %8.3f %8.3f" % (d[j], I[j], sig[j]))

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
                    label = f"image {indx}" if len(images) > 1 else f""
                ds2 = 1 / flex.pow2(d)
                ax.plot(ds2, I, label=label)
            xticks = ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(xticks))
            x_tick_labs = [
                "" if e <= 0.0 else "{:.2f}".format(math.sqrt(1.0 / e)) for e in xticks
            ]
            ax.set_xticklabels(x_tick_labs)

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


def background(imageset, indx, n_bins, corrected=False, mask_params=None):
    if mask_params is None:
        # Default mask params for trusted range
        mask_params = phil_scope.fetch(parse("")).extract().masking

    mask = dials.util.masking.generate_mask(imageset, mask_params)

    detector = imageset.get_detector()
    beam = imageset.get_beam()
    # Only working with single panel detector for now
    assert len(detector) == 1
    panel = detector[0]
    mask = mask[0]

    n = matrix.col(panel.get_normal()).normalize()
    b = matrix.col(beam.get_s0()).normalize()
    wavelength = beam.get_wavelength()

    if math.fabs(b.dot(n)) < 0.95:
        raise Sorry("Detector not perpendicular to beam")

    if corrected:
        data = imageset.get_corrected_data(indx)
    else:
        data = imageset.get_raw_data(indx)
    assert len(data) == 1
    data = data[0]

    data = data.as_double()

    spot_params = spot_phil.fetch(source=parse("")).extract()
    threshold_function = SpotFinderFactory.configure_threshold(spot_params)
    peak_pixels = threshold_function.compute_threshold(data, mask)
    signal = data.select(peak_pixels.iselection())
    background_pixels = mask & ~peak_pixels
    background = data.select(background_pixels.iselection())

    # print some summary information
    print("Mean background: %.3f" % (flex.sum(background) / background.size()))
    if len(signal) > 0:
        print(
            "Max/total signal pixels: %.0f / %.0f"
            % (flex.max(signal), flex.sum(signal))
        )
    else:
        print("No signal pixels on this image")
    print(
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
