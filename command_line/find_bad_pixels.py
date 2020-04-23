# LIBTBX_SET_DISPATCHER_NAME dials.find_bad_pixels

from __future__ import absolute_import, division, print_function

import concurrent.futures
import math
import pickle
import sys

import iotbx.phil
from scitbx.array_family import flex
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments
from dxtbx.model.experiment_list import ExperimentList, Experiment

help_message = """

Examples::

  dials.find_bad_pixels data_master.h5 [nproc=8]

"""

phil_scope = iotbx.phil.parse(
    """
images = None
  .type = ints
  .help = "Images on which to perform the analysis (otherwise use all images)"
image_range = None
  .type = ints(value_min=0, size=2)
  .help = "Image range for analysis e.g. 1,1800"
nproc = 1
  .type = int(value_min=1)
    .help = "The number of processes to use."

output {
    mask = pixels.mask
        .type = path
        .help = "Output mask file name"
    png = pixels.png
        .type = path
        .help = "Bad pixel mask as image"
    print_values = False
        .type = bool
        .help = "Print bad pixel values"
}
"""
)


def find_constant_signal_pixels(imageset, images):
    """Find pixels which are constantly reporting as signal through the
    images in imageset: on every image the pixel dispersion index is computed,
    and signal pixels identified using the default settings. A map is then
    calculated of the number of times a pixel is identified as signal: if this
    is >= 50% of the images (say) that pixel is untrustworthy."""

    panels = imageset.get_detector()

    # only cope with monilithic detectors or the I23 Pilatus 12M
    assert len(panels) in (1, 24)

    # trusted range the same for all panels anyway
    detector = panels[0]
    trusted = detector.get_trusted_range()

    # construct an integer array same shape as image; accumulate number of
    # "signal" pixels in each pixel across data

    total = None

    for idx in images:
        pixels = imageset.get_raw_data(idx - 1)
        known_mask = imageset.get_mask(idx - 1)

        # apply known mask
        for _pixel, _panel, _mask in zip(pixels, panels, known_mask):
            _pixel.set_selected(~_mask, -1)
            for f0, s0, f1, s1 in _panel.get_mask():
                blank = flex.int(flex.grid(s1 - s0, f1 - f0), 0)
                _pixel.matrix_paste_block_in_place(blank, s0, f0)

        if len(pixels) == 1:
            data = pixels[0]
        else:
            ny, nx = pixels[0].focus()
            data = flex.int(flex.grid(24 * ny + 23 * 17, nx), -1)
            for j in range(24):
                data.matrix_paste_block_in_place(pixels[j], j * (ny + 17), 0)

        negative = data < int(round(trusted[0]))
        hot = data > int(round(trusted[1]))
        bad = negative | hot

        data = data.as_double()

        spot_params = spot_phil.fetch(
            source=iotbx.phil.parse("min_spot_size=1")
        ).extract()
        threshold_function = SpotFinderFactory.configure_threshold(
            spot_params,
            ExperimentList(
                [
                    Experiment(
                        beam=imageset.get_beam(),
                        detector=imageset.get_detector(),
                        goniometer=imageset.get_goniometer(),
                        scan=imageset.get_scan(),
                        imageset=imageset,
                    )
                ]
            ),
        )
        peak_pixels = threshold_function.compute_threshold(data, ~bad)

        if total is None:
            total = peak_pixels.as_1d().as_int()
        else:
            total += peak_pixels.as_1d().as_int()

    return total


def run(args):
    usage = "dev.dials.find_bad_pixels [options] data_master.h5"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
        parser.print_help()
        sys.exit("Please pass an experiment list\n")
        return

    imagesets = experiments.imagesets()

    if len(imagesets) != 1:
        sys.exit("Please pass an experiment list that contains one imageset")

    imageset = imagesets[0]
    panels = imageset.get_detector()
    detector = panels[0]
    nfast, nslow = detector.get_image_size()

    first, last = imageset.get_scan().get_image_range()
    images = range(first, last + 1)

    if params.images is None and params.image_range is not None:
        start, end = params.image_range
        params.images = list(range(start, end + 1))

    if params.images:
        if min(params.images) < first or max(params.images) > last:
            sys.exit("Image outside of scan range")
        images = params.images

    # work around issues with HDF5 and multiprocessing
    if hasattr(imageset.reader(), "nullify_format_instance"):
        imageset.reader().nullify_format_instance()

    n = int(math.ceil(len(images) / params.nproc))
    chunks = [images[i : i + n] for i in range(0, len(images), n)]

    assert len(images) == sum(len(chunk) for chunk in chunks)

    if len(chunks) < params.nproc:
        params.nproc = len(chunks)

    total = None

    with concurrent.futures.ProcessPoolExecutor(max_workers=params.nproc) as p:
        jobs = []
        for j in range(params.nproc):
            jobs.append(p.submit(find_constant_signal_pixels, imageset, chunks[j]))
        for job in concurrent.futures.as_completed(jobs):
            if total is None:
                total = job.result()
            else:
                total += job.result()

    hot_mask = total >= (len(images) // 2)
    hot_pixels = hot_mask.iselection()

    if params.output.mask:
        with open(params.output.mask, "wb") as fh:
            hot_mask.reshape(flex.grid(reversed(detector.get_image_size())))
            pickle.dump(~hot_mask, fh)

    if params.output.print_values:
        for h in hot_pixels:
            print("    mask[%d, %d] = 8" % (h % nfast, h // nfast))


if __name__ == "__main__":
    run(sys.argv[1:])
