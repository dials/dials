# LIBTBX_SET_DISPATCHER_NAME dev.dials.csv

from __future__ import annotations

import gzip
import io
import warnings

import iotbx.phil
from dxtbx.model import ExperimentList

import dials.util
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

phil_scope = iotbx.phil.parse(
    """
output {
  csv = rl.csv
    .type = path
    .help = 'Output filename for reciprocal mapped reflections'
  compress = False
    .type = bool
    .help = 'Compress file with gzip as written'
  dp = 0
    .type = int
    .help = 'Decimal places for output, 0 => %f'
}
"""
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dev.dials.csv [options] imported.expt strong.refl output.csv=rl.csv"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, options = parser.parse_args(args, show_diff_phil=False)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if not experiments or not reflections:
        parser.print_help()
        exit(0)

    spots = []

    for reflection in reflections:
        unique_ids = set(reflection["id"])
        for unique_id in sorted(unique_ids):
            spots.append(reflection.select(reflection["id"] == unique_id))
        if not reflection:  # If there are no reflections then export an empty list
            spots.append(reflection)

    assert len(experiments) == len(spots)

    if params.output.compress:
        fout = gzip.GzipFile(params.output.csv, "w")
        # GzipFile() always provides binary access only.
        # Replace the file object with one that allows writing text:
        fout = io.TextIOWrapper(fout)
        # Rely on garbage collection to close the underlying GzipFile.
    else:
        fout = open(params.output.csv, "w")

    fout.write("# x,y,z,experiment_id,imageset_id\n")

    dp = params.output.dp

    if dp <= 0:
        fmt = "%f,%f,%f,%d,%d\n"
    else:
        fmt = "%%.%df,%%.%df,%%.%df,%%d,%%d\n" % (dp, dp, dp)

    print("Using format:", fmt.strip())

    for k, (expt, refl) in enumerate(zip(experiments, spots)):
        if "imageset_id" not in refl:
            refl["imageset_id"] = refl["id"]

        refl.centroid_px_to_mm(ExperimentList([expt]))
        refl.map_centroids_to_reciprocal_space(ExperimentList([expt]))
        rlp = refl["rlp"]

        for _rlp in rlp:
            fout.write(fmt % (_rlp[0], _rlp[1], _rlp[2], k, k))

        print(f"Appended {len(rlp)} spots to {params.output.csv}")

    fout.close()


if __name__ == "__main__":
    warnings.warn(
        "dev.dials.csv is deprecated. Similar functionality is available"
        " with dials.export format=json",
        DeprecationWarning,
    )
    run()
