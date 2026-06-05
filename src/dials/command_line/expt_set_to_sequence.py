"""Convert a per-frame ImageSet .expt file to the stills-as-ImageSequence format.

Reads the classic format (one ImageSet per experiment, one monochromatic Beam per
experiment, no scans) and writes the new format produced by the
stills_process_imagesequence branch: one shared ImageSequence per source file, one
XFELBeam, detectors deduplicated, and a single consolidated scan object with
per-frame wavelengths.

Experiments are sorted by 0-based HDF5 frame index (ascending) so that each
imageset's ``single_file_indices`` is sparse, ascending, and aligns with its
experiments in scan_point order (the consolidated reader's alignment invariant).
If a companion reflection table is provided, the ``id`` column and
``experiment_identifiers`` map are remapped to match the new experiment ordering.

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

    # A source file is "real" only if every imageset referencing it carries
    # single_file_indices. Legacy per-frame ImageSet lists (e.g. older merge
    # inputs) dropped single_file_indices and reused one dummy/empty source path
    # for every still, so frame identity was lost on disk. For such a source we
    # cannot recover the original HDF5 frame index -- but each experiment still
    # carries its own beam wavelength, so we synthesize a distinct 0-based frame
    # index per experiment (in original order) and treat each as its own unique
    # frame. Collapsing them to a single (src, fi) would fuse N distinct per-shot
    # wavelengths into one and corrupt every downstream wavelength read.
    src_has_sfi = {}
    for exp in experiments:
        iset = imagesets_in[exp["imageset"]]
        src = _source_file(iset)
        has = bool(iset.get("single_file_indices"))
        src_has_sfi[src] = src_has_sfi.get(src, True) and has

    # Gather per-experiment metadata, then sort by (src, fi).
    synthetic_fi = {}
    exp_info = []
    for i, exp in enumerate(experiments):
        iset = imagesets_in[exp["imageset"]]
        src = _source_file(iset)
        if src_has_sfi[src]:
            fi = _frame_index(iset)
        else:
            fi = synthetic_fi.get(src, 0)
            synthetic_fi[src] = fi + 1
        exp_info.append(
            {
                "old_idx": i,
                "exp": exp,
                "src": src,
                "fi": fi,
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
        # Sparse, ascending frame indices: one entry per unique frame this
        # source actually contributes. (Dense range() listed frames no
        # experiment references and broke positional alignment with scan_points.)
        fis = sorted(src_frame_sets[src])
        first_iset = next(
            imagesets_in[exp["imageset"]]
            for exp in experiments
            if _source_file(imagesets_in[exp["imageset"]]) == src
        )
        new_iset = {
            "__id__": "ImageSequence",
            "template": src,
            "single_file_indices": fis,
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

    # scan_point is a per-UNIQUE-FRAME global index: multi-lattice experiments on
    # the same (src, fi) share one scan_point (per-shot wavelength is identical
    # across lattices, so dedup is safe). It is distinct from new_idx, the
    # per-EXPERIMENT index for the crystal/profile/scaling_model/experiment
    # arrays (one new experiment per lattice). Assign scan_points per imageset in
    # ascending-frame order so that within each imageset ascending scan_point
    # corresponds to ascending frame -- the alignment the consolidated reader
    # relies on. The consolidated wavelength array is indexed by scan_point, so
    # it is built here in scan_point order (one entry per unique frame).
    wl_for_frame = {(info["src"], info["fi"]): info["wl"] for info in exp_info}
    scan_point_for_frame = {}
    consolidated_sfi = []
    wavelengths = []
    for src in src_frame_sets:
        for fi in sorted(src_frame_sets[src]):
            scan_point_for_frame[(src, fi)] = len(wavelengths)
            consolidated_sfi.append(fi)
            wavelengths.append(wl_for_frame[(src, fi)])

    # Build sorted experiments and per-experiment model lists.
    new_crystals = []
    new_profiles = []
    new_scaling_models = []
    new_experiments = []

    for new_idx, info in enumerate(exp_info):
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
        exp["scan_point"] = scan_point_for_frame[(info["src"], info["fi"])]
        if new_crystals:
            exp["crystal"] = new_idx
        if new_profiles:
            exp["profile"] = new_idx
        if new_scaling_models:
            exp["scaling_model"] = new_idx
        exp.pop("goniometer", None)
        new_experiments.append(exp)

    consolidated_scan = build_consolidated_scan(consolidated_sfi, wavelengths)

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
