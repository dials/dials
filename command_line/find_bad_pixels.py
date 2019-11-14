# LIBTBX_SET_DISPATCHER_NAME dev.dials.find_bad_pixels

from __future__ import absolute_import, division, print_function

import concurrent.futures
import math
import sys

import iotbx.phil
from scitbx.array_family import flex
from dials.util import Sorry
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments
from dxtbx.model.experiment_list import ExperimentList, Experiment
from libtbx import easy_pickle

help_message = """

Examples::

  dev.dials.find_bad_pixels data_master.h5

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
}
"""
)


def find_constant_signal_pixels(imageset, images):
    """Find pixels which are constantly reporting as signal through the
    images in imageset."""

    detectors = imageset.get_detector()
    assert len(detectors) == 1
    detector = detectors[0]
    trusted = detector.get_trusted_range()

    # construct an integer array same shape as image; accumulate number of
    # "signal" pixels in each pixel across data

    total = None

    for idx in images:
        pixels = imageset.get_raw_data(idx - 1)
        assert len(pixels) == 1
        data = pixels[0]

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

    params, options = parser.parse_args(show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
        parser.print_help()
        sys.exit("Please pass an experiment list\n")
        return

    imagesets = experiments.imagesets()

    if len(imagesets) != 1:
        sys.exit("Please pass an experiment list that contains one imageset")

    imageset = imagesets[0]
    detectors = imageset.get_detector()
    assert len(detectors) == 1
    detector = detectors[0]
    trusted = detector.get_trusted_range()

    first, last = imageset.get_scan().get_image_range()
    images = range(first, last + 1)

    if params.images is None and params.image_range is not None:
        start, end = params.image_range
        params.images = list(range(start, end + 1))

    if params.images:
        if min(params.images) < first or max(params.images) > last:
            raise Sorry("image outside of scan range")
        images = params.images

    # work around issues with HDF5 and multiprocessing
    imageset.reader().nullify_format_instance()

    n = int(math.ceil(len(images) / params.nproc))
    work = [images[i : i + n] for i in range(0, len(images), n)]
    assert len(images) == sum([len(chunk) for chunk in work])

    total = None

    with concurrent.futures.ProcessPoolExecutor(max_workers=params.nproc) as p:
        jobs = []
        for j in range(params.nproc):
            jobs.append(p.submit(find_constant_signal_pixels, imageset, work[j]))
        for job in concurrent.futures.as_completed(jobs):
            if total is None:
                total = job.result()
            else:
                total += job.result()

    hot_mask = total >= (len(images) // 2)
    hot_pixels = hot_mask.iselection()

    capricious_pixels = {}
    for h in hot_pixels:
        capricious_pixels[h] = []

    for idx in images:
        pixels = imageset.get_raw_data(idx - 1)
        data = pixels[0]

        for h in hot_pixels:
            capricious_pixels[h].append(data[h])

    nslow, nfast = data.focus()

    ffff = 0

    for h in hot_pixels:
        if total[h] == len(images) and data[h] >= trusted[1]:
            ffff += 1
            continue
        print("Pixel %d at %d %d" % (total[h], h // nfast, h % nfast))
        if len(set(capricious_pixels[h])) >= len(capricious_pixels[h]) // 2:
            print("  ... many possible values")
            continue
        values = set(capricious_pixels[h])
        result = [(capricious_pixels[h].count(value), value) for value in values]
        for count, value in reversed(sorted(result)):
            print("  %08x %d" % (value, count))

    print("Also found %d very hot pixels" % ffff)
    hot_mask.reshape(flex.grid(data.focus()))

    easy_pickle.dump(params.output.mask, (~hot_mask,))


if __name__ == "__main__":
    run(sys.argv[1:])
