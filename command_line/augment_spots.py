from __future__ import annotations

import iotbx.phil

import dials.util
from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

help_message = """
Augment spot list with additional information - for example number of pixels
in peak region etc."""

phil_scope = iotbx.phil.parse(
    """
output {
  reflections = stronger.refl
    .type = path
}
find_max = False
  .type = bool
"""
)


def add_max_pixels_to_reflections(reflections):
    """Iterate through reflections, find max pixel in each shoebox and max
    valid, add columns thereof"""

    from dials.algorithms.shoebox import MaskCode

    good = MaskCode.Foreground | MaskCode.Valid

    max_pixel = flex.double(reflections.size(), 0.0)
    max_valid = flex.double(reflections.size(), 0.0)

    shoeboxes = reflections["shoebox"]

    for j, s in enumerate(shoeboxes):
        max_pixel[j] = flex.max(s.data)
        valid = s.mask.as_1d() == good
        max_valid[j] = flex.max(s.data.as_1d().select(valid))

    reflections["max_pixel"] = max_pixel
    reflections["max_valid"] = max_valid


def add_resolution_to_reflections(reflections, experiments):
    """Add column d to reflection list"""

    # will assume everything from the first detector at the moment - clearly this
    # could be incorrect, will have to do something a little smarter, later

    if "imageset_id" not in reflections:
        reflections["imageset_id"] = reflections["id"]

    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    reflections["d"] = 1 / reflections["rlp"].norms()


def augment_reflections(reflections, params, experiments=None):
    """Add extra columns of derived data."""

    from dials.algorithms.shoebox import MaskCode

    good = MaskCode.Foreground | MaskCode.Valid

    if params.find_max:
        add_max_pixels_to_reflections(reflections)

    x0, x1, y0, y1, z0, z1 = reflections["bbox"].parts()

    dx = x1 - x0
    dy = y1 - y0

    # compute signal pixels in each shoebox as an array
    n_signal = reflections["shoebox"].count_mask_values(good)

    reflections["dx"] = dx
    reflections["dy"] = dy
    reflections["n_signal"] = n_signal

    if experiments:
        add_resolution_to_reflections(reflections, experiments)

    return reflections


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util import Sorry

    usage = "dials.augment_spots [options] [models.expt] strong.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(reflections) != 1:
        parser.print_help()
        return
    if len(experiments) > 1:
        raise Sorry("0, 1 experiments required")

    stronger = augment_reflections(reflections[0], params, experiments=experiments)
    stronger.as_file(params.output.reflections)


if __name__ == "__main__":
    run()
