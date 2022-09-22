from __future__ import annotations

import os
import pickle

import iotbx.phil
from dxtbx.imageset import ImageSet
from libtbx.math_utils import nearest_integer as nint
from scitbx.array_family import flex

from dials.algorithms.image.threshold import DispersionThresholdDebug
from dials.util import Sorry, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """

This program can be used to estimate the gain of the detector. For pixel array
detectors the gain is usually set to 1.00. This means that the pixels behave
according to Poisson statistics. However, for older CCD detectors the gain may
have a different value. This value is important because it can affect, amongst
other things, the ability of the spot finding algorithm which can result in
noise being identified as diffraction spots.

Examples::

  dials.estimate_gain models.expt
"""

phil_scope = iotbx.phil.parse(
    """\
  kernel_size = 10,10
    .type = ints(size=2, value_min=1)
  max_images = 1
    .type = int
    .help = "For multi-file images (NeXus for example), report a gain for each"
            "image, up to max_images, and then report an average gain"
  output {
    gain_map = None
      .type = str
      .help = "Name of output gain map file"
  }
""",
    process_includes=True,
)


def estimate_gain(
    imageset: ImageSet,
    kernel_size: tuple[int, int] = (10, 10),
    output_gain_map: str | bytes | os.PathLike | None = None,
    max_images: int = 1,
) -> float:
    detector = imageset.get_detector()

    gains = flex.double()

    for image_no in range(len(imageset)):
        raw_data = imageset.get_raw_data(image_no)

        gain_value = 1
        gain_map = [
            flex.double(raw_data[i].accessor(), gain_value)
            for i in range(len(detector))
        ]

        mask = imageset.get_mask(image_no)

        min_local = 0

        # dummy values, shouldn't affect results
        nsigma_b = 6
        nsigma_s = 3
        global_threshold = 0

        kabsch_debug_list = [
            DispersionThresholdDebug(
                raw_data[i_panel].as_double(),
                mask[i_panel],
                gain_map[i_panel],
                kernel_size,
                nsigma_b,
                nsigma_s,
                global_threshold,
                min_local,
            )
            for i_panel in range(len(detector))
        ]

        dispersion = flex.double()
        for kabsch in kabsch_debug_list:
            dispersion.extend(kabsch.index_of_dispersion().as_1d())

        sorted_dispersion = flex.sorted(dispersion)

        q1 = sorted_dispersion[nint(len(sorted_dispersion) / 4)]
        q2 = sorted_dispersion[nint(len(sorted_dispersion) / 2)]
        q3 = sorted_dispersion[nint(len(sorted_dispersion) * 3 / 4)]
        iqr = q3 - q1

        print(f"q1, q2, q3: {q1:.2f}, {q2:.2f}, {q3:.2f}")
        if iqr == 0.0:
            raise Sorry("Unable to robustly estimate the variation of pixel values.")

        inlier_sel = (sorted_dispersion > (q1 - 1.5 * iqr)) & (
            sorted_dispersion < (q3 + 1.5 * iqr)
        )
        sorted_dispersion = sorted_dispersion.select(inlier_sel)
        gain = sorted_dispersion[nint(len(sorted_dispersion) / 2)]
        print(f"Estimated gain: {gain:.2f}")
        gains.append(gain)

        if image_no == 0:
            gain0 = gain
        if image_no + 1 >= max_images:
            break

    if len(gains) > 1:
        stats = flex.mean_and_variance(gains)
        print(
            "Average gain: %.2f +/- %.2f"
            % (stats.mean(), stats.unweighted_sample_standard_deviation())
        )

    if output_gain_map:
        if len(gains) > 1:
            raw_data = imageset.get_raw_data(0)
        # write the gain map
        gain_map_flex = flex.double(flex.grid(raw_data[0].all()), gain0)
        with open(output_gain_map, "wb") as fh:
            pickle.dump(gain_map_flex, fh, protocol=pickle.HIGHEST_PROTOCOL)

    return gain0


@show_mail_handle_errors()
def run(args: list[str] | None = None) -> None:
    usage = "dials.estimate_gain [options] models.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, _ = parser.parse_args(args, show_diff_phil=False)

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        print("The following parameters have been modified:\n")
        print(diff_phil)

    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        return
    elif len(experiments.imagesets()) > 1:
        raise Sorry("Only one imageset can be processed at a time")
    else:
        imagesets = experiments.imagesets()

    assert len(imagesets) == 1
    imageset = imagesets[0]
    estimate_gain(
        imageset, params.kernel_size, params.output.gain_map, params.max_images
    )


if __name__ == "__main__":
    run()
