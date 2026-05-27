"""Convert a per-frame ImageSet .expt file to the stills-as-ImageSequence format.

Reads the classic format (one ImageSet per experiment, one monochromatic Beam per
experiment, no scans) and writes the new format produced by the
stills_process_imagesequence branch: one shared ImageSequence per source file, one
XFELBeam, detectors deduplicated, and a single consolidated scan object with
per-frame wavelengths.

Experiments are sorted by 0-based HDF5 frame index (ascending) so that the
``frame_numbers`` array in the consolidated scan is ordered. If a companion
reflection table is provided, the ``id`` column and ``experiment_identifiers``
map are remapped to match the new experiment ordering.

Examples::

  dials.expt_set_to_sequence integrated.expt
  dials.expt_set_to_sequence integrated.expt output.experiments_filename=seq.expt
  dials.expt_set_to_sequence integrated.expt integrated.refl \\
    output.experiments_filename=seq.expt output.reflections_filename=seq.refl
"""

from __future__ import annotations

import copy
import json
import logging
import sys
from collections import OrderedDict

from libtbx.phil import parse

import dials.util
from dials.util import log
from dials.util.options import ArgumentParser
from dials.util.stills_imageset_convert import (
    build_consolidated_scan,
)
from dials.util.stills_imageset_convert import (
    frame_index as _frame_index,
)
from dials.util.stills_imageset_convert import (
    source_file as _source_file,
)

logger = logging.getLogger(__name__)

help_message = __doc__

phil_scope = parse(
    """
  output {
    experiments_filename = sequence.expt
      .type = str
      .help = "Output filename for the ImageSequence-format experiments"

    reflections_filename = sequence.refl
      .type = str
      .help = "Output filename for the remapped reflection table"

    log = dials.expt_set_to_sequence.log
      .type = str
  }
"""
)


def convert_expt(input_path, output_path):
    """Convert JSON, return old_to_new permutation dict for refl remapping."""
    with open(input_path) as f:
        d = json.load(f)

    imagesets_in = d.get("imageset", [])
    beams_in = d.get("beam", [])
    detectors_in = d.get("detector", [])
    crystals_in = d.get("crystal", [])
    profiles_in = d.get("profile", [])
    scaling_models_in = d.get("scaling_model", [])
    experiments = d["experiment"]

    # Gather per-experiment metadata, then sort by (src, fi).
    exp_info = []
    for i, exp in enumerate(experiments):
        iset = imagesets_in[exp["imageset"]]
        exp_info.append(
            {
                "old_idx": i,
                "exp": exp,
                "src": _source_file(iset),
                "fi": _frame_index(iset),
                "wl": beams_in[exp["beam"]].get("wavelength", 0.0),
                "old_det": exp.get("detector", 0),
                "old_crystal": exp.get("crystal"),
                "old_profile": exp.get("profile"),
                "old_scaling": exp.get("scaling_model"),
            }
        )

    exp_info.sort(key=lambda x: (x["src"], x["fi"]))
    old_to_new = {info["old_idx"]: new_idx for new_idx, info in enumerate(exp_info)}

    # One ImageSequence per unique source file with a dense frame-index range.
    src_to_iset_idx = OrderedDict()
    new_imagesets = []
    src_frame_sets = {}
    for info in exp_info:
        src_frame_sets.setdefault(info["src"], set()).add(info["fi"])

    for src in src_frame_sets:
        fis = sorted(src_frame_sets[src])
        dense = list(range(fis[0], fis[-1] + 1))
        first_iset = next(
            imagesets_in[exp["imageset"]]
            for exp in experiments
            if _source_file(imagesets_in[exp["imageset"]]) == src
        )
        new_iset = {
            "__id__": "ImageSequence",
            "template": src,
            "single_file_indices": dense,
            "mask": first_iset.get("mask"),
            "gain": first_iset.get("gain"),
            "pedestal": first_iset.get("pedestal"),
            "dx": first_iset.get("dx"),
            "dy": first_iset.get("dy"),
            "params": first_iset.get("params", {}),
        }
        src_to_iset_idx[src] = len(new_imagesets)
        new_imagesets.append(new_iset)

    # One XFELBeam: strip wavelength from beam[0], change __id__.
    xfel_beam = copy.deepcopy(beams_in[0])
    xfel_beam["__id__"] = "xfel"
    xfel_beam.pop("wavelength", None)

    # Deduplicate detectors by JSON serialisation.
    det_key_to_new_idx = {}
    new_detectors = []
    old_to_new_det = {}
    for old_idx, det in enumerate(detectors_in):
        key = json.dumps(det, sort_keys=True)
        if key not in det_key_to_new_idx:
            det_key_to_new_idx[key] = len(new_detectors)
            new_detectors.append(copy.deepcopy(det))
        old_to_new_det[old_idx] = det_key_to_new_idx[key]

    # Build sorted experiments and per-experiment model lists.
    frame_numbers = []
    wavelengths = []
    new_crystals = []
    new_profiles = []
    new_scaling_models = []
    new_experiments = []

    for new_idx, info in enumerate(exp_info):
        frame_numbers.append(info["fi"] + 1)
        wavelengths.append(info["wl"])

        if info["old_crystal"] is not None and crystals_in:
            new_crystals.append(copy.deepcopy(crystals_in[info["old_crystal"]]))
        if info["old_profile"] is not None and profiles_in:
            new_profiles.append(copy.deepcopy(profiles_in[info["old_profile"]]))
        if info["old_scaling"] is not None and scaling_models_in:
            new_scaling_models.append(
                copy.deepcopy(scaling_models_in[info["old_scaling"]])
            )

        exp = copy.deepcopy(info["exp"])
        exp["imageset"] = src_to_iset_idx[info["src"]]
        exp["beam"] = 0
        exp["detector"] = old_to_new_det[info["old_det"]]
        exp["scan"] = 0
        exp["scan_point"] = new_idx
        if new_crystals:
            exp["crystal"] = new_idx
        if new_profiles:
            exp["profile"] = new_idx
        if new_scaling_models:
            exp["scaling_model"] = new_idx
        exp.pop("goniometer", None)
        new_experiments.append(exp)

    consolidated_scan = build_consolidated_scan(frame_numbers, wavelengths)

    d["experiment"] = new_experiments
    d["imageset"] = new_imagesets
    d["beam"] = [xfel_beam]
    d["detector"] = new_detectors
    d["scan"] = [consolidated_scan]
    d["goniometer"] = []
    if new_crystals:
        d["crystal"] = new_crystals
    if new_profiles:
        d["profile"] = new_profiles
    if new_scaling_models:
        d["scaling_model"] = new_scaling_models

    with open(output_path, "w") as f:
        json.dump(d, f, indent=2)

    m = len(new_imagesets)
    logger.info(
        f"Wrote {len(new_experiments)} experiments as ImageSequence format to "
        f"{output_path} ({m} imageset{'s' if m != 1 else ''} from "
        f"{m} source file{'s' if m != 1 else ''})"
    )
    return old_to_new, new_experiments


def remap_reflections(input_refl, output_refl, old_to_new, new_experiments):
    """Remap experiment_id column and identifier map to match sorted order."""
    from dials.array_family import flex

    refl = flex.reflection_table.from_file(input_refl)

    refl["id"] = flex.int([old_to_new[i] for i in refl["id"]])

    eid = refl.experiment_identifiers()
    old_uuid = {k: eid[k] for k in eid.keys()}
    for old_idx, new_idx in old_to_new.items():
        if old_idx in old_uuid:
            eid[new_idx] = old_uuid[old_idx]

    refl.as_file(output_refl)
    logger.info(f"Wrote remapped reflections to {output_refl}")


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = (
        "dials.expt_set_to_sequence [options] integrated.expt [integrated.refl] "
        "[output.experiments_filename=sequence.expt]"
    )

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )
    params, options = parser.parse_args(args, show_diff_phil=True)
    log.config(verbosity=options.verbose, logfile=params.output.log)

    if not params.input.experiments:
        logger.info("No experiments found in the input")
        parser.print_help()
        return

    if len(params.input.experiments) > 1:
        sys.exit("Error: provide a single .expt file as input")

    input_expt = params.input.experiments[0].filename
    output_expt = params.output.experiments_filename

    old_to_new, new_experiments = convert_expt(input_expt, output_expt)

    if params.input.reflections:
        if len(params.input.reflections) > 1:
            sys.exit("Error: provide at most one .refl file as input")
        remap_reflections(
            params.input.reflections[0].filename,
            params.output.reflections_filename,
            old_to_new,
            new_experiments,
        )


if __name__ == "__main__":
    run()
