# DIALS_ENABLE_COMMAND_LINE_COMPLETION
"""
Split dials still-shot data files (expt and refl files) into separate files
based on image number or metadata. For example, splitting a dose series dataset.
Either the repeat unit can be specified with the series_repeat option, or
a yaml file can be provided.
If series_repeat=10, then data corresponding to every nth image will be split into
separate files:
group_0_0.expt
group_1_0.expt
....
group_9_0.expt
and the same for reflection files.
The second index relates to the input file index, so if more than one input file is input,
there will then be group_0_1.expt etc.

Example yaml syntax is
Point to a data array:

---
metadata:
  timepoint:
    '/path/to/example_master.h5' : '/path/to/meta.h5:/timepoint'
grouping:
  split_by:
    values:
      - timepoint

or for images specified by a template:
---
metadata:
  dose_point:
    '/path/to/example_#####.cbf' : 'repeat=10'
grouping:
  split_by:
    values:
      - dose_point

"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import iotbx.phil
from dxtbx.serialize import load
from libtbx import Auto

from dials.util import log
from dials.util.image_grouping import (
    FilePair,
    ParsedYAML,
    get_grouping_handler,
    series_repeat_to_groupings,
)
from dials.util.options import ArgumentParser
from dials.util.system import CPU_COUNT
from dials.util.version import dials_version

logger = logging.getLogger("dials")


phil_str = """
input {
  reflections = None
    .type = str
    .multiple = True
    .help = "Path to an integrated reflections file"
    .expert_level = 1
  experiments = None
    .type = str
    .multiple = True
    .help = "Path to an integrated experiments file"
    .expert_level = 1
}
grouping = None
  .type = str
  .help = "Path to a .yml file defining grouping structure during processing"
series_repeat = None
  .type = int(value_min=2)
  .expert_level = 2
  .help = "This option allows the user to specify that the data is a repeated series"
          "by providing the number of repeated measurements at each point. i.e. it"
          "is assumed that $series_repeat measurements are taken at each position"
          "and that these form consecutive images in the input image files."
nproc=Auto
  .type=int
output.log=dials.split_ssx_data.log
  .type=str
"""

phil_scope = iotbx.phil.parse(phil_str)


def run(args=sys.argv[1:]):
    parser = ArgumentParser(
        usage="dials.split_ssx_data integrated.expt integrated.refl series_repeat=10",
        read_experiments=False,
        read_reflections=False,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    params, options, unhandled = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )
    if unhandled:
        for item in unhandled:
            if item.endswith(".expt"):
                args[args.index(item)] = f"input.experiments = {item}"
            elif item.endswith(".refl"):
                args[args.index(item)] = f"input.reflections = {item}"
            else:
                raise ValueError(f"Unhandled argument: {item}")
        params, options = parser.parse_args(args=args, show_diff_phil=False)

    if not (params.input.reflections and params.input.experiments):
        raise ValueError("Reflections and experiments files must both be specified")
    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    reflection_files = sorted([Path(i).resolve() for i in params.input.reflections])
    experiment_files = sorted([Path(i).resolve() for i in params.input.experiments])

    datafiles = [FilePair(e, r) for e, r in zip(experiment_files, reflection_files)]
    expts = []
    for fp in datafiles:
        elist = load.experiment_list(fp.expt, check_format=False)
        for e in elist:
            if e.scan and e.scan.get_oscillation()[1] != 0.0:
                logger.info("dials.split_ssx_data can only be used on still-shot data")
                sys.exit()
        expts.append(elist)

    if params.nproc is Auto:
        params.nproc = CPU_COUNT
        logger.info(f"Using nproc={params.nproc}")

    if params.series_repeat:
        parsed = series_repeat_to_groupings(
            expts, params.series_repeat, groupname="split_by"
        )
        handler = get_grouping_handler(parsed, "split_by", nproc=params.nproc)
    elif params.grouping:
        try:
            parsed = ParsedYAML(params.grouping)
        except Exception as e:
            logger.info(f"Error: {e}\nPlease check input in yaml file.")
        groupings = parsed.groupings
        if len(groupings) != 1:
            logger.info("Only one grouping type can be specified in grouping yaml")
            sys.exit()
        handler = get_grouping_handler(
            parsed, list(parsed.groupings.keys())[0], nproc=params.nproc
        )
    try:
        _ = handler.split_files_to_groups(
            working_directory=Path.cwd(),
            data_file_pairs=datafiles,
        )
        logger.info("Finished splitting data.")
    except Exception as e:
        logger.info(f"Error: {e}\nPlease check input.")


if __name__ == "__main__":
    run()
