"""Convert a stills-as-ImageSequence .expt file to the old per-frame ImageSet format.

Reads the new format produced by the stills_process_imagesequence branch (one shared
ImageSequence per source file, one XFELBeam, a consolidated scan with per-frame
wavelengths) and writes the classic format expected by legacy downstream tools: one
ImageSet per experiment, one monochromatic Beam per experiment, no scans.

Experiment order is preserved, so companion reflection tables (.refl files) do not
require any remapping — experiment_id values remain valid as-is.

Examples::

  dials.expt_sequence_to_set integrated.expt
  dials.expt_sequence_to_set integrated.expt output.experiments_filename=old.expt
  dials.expt_sequence_to_set integrated.expt integrated.refl \\
    output.experiments_filename=old.expt output.reflections_filename=old.refl
"""

from __future__ import annotations

import copy
import json
import logging
import shutil
import sys

import dials.util
from dials.util import log
from dials.util.options import ArgumentParser
from libtbx.phil import parse

logger = logging.getLogger(__name__)

help_message = __doc__

phil_scope = parse(
    """
  output {
    experiments_filename = imageset.expt
      .type = str
      .help = "Output filename for the ImageSet-format experiments"

    reflections_filename = imageset.refl
      .type = str
      .help = "Output filename for the (unchanged) reflection table"

    log = dials.expt_sequence_to_set.log
      .type = str
  }
"""
)


def _expand_consolidated_scan(scan_list):
    """Return (frame_numbers, wavelengths) from a __stills_consolidated scan list.

    Returns (None, None) if the scan list is absent or not consolidated.
    """
    if not scan_list or not scan_list[0].get("__stills_consolidated"):
        return None, None
    con = scan_list[0]
    fns = con["frame_numbers"]
    wls = con.get("properties", {}).get("wavelength", [None] * len(fns))
    return fns, wls


def _per_frame_scan_map(scan_list):
    """Build {scan_index: (fi, wl)} from non-consolidated per-frame scans.

    fi = 0-based HDF5 frame index = image_range[0] - 1.
    wl = wavelength from scan properties, or None.
    """
    result = {}
    for i, s in enumerate(scan_list):
        fi = s["image_range"][0] - 1
        wls = s.get("properties", {}).get("wavelength", [])
        result[i] = (fi, wls[0] if wls else None)
    return result


def convert_expt(input_path, output_path):
    with open(input_path) as f:
        d = json.load(f)

    frame_numbers, wavelengths = _expand_consolidated_scan(d.get("scan", []))
    consolidated = frame_numbers is not None
    if not consolidated:
        per_frame = _per_frame_scan_map(d.get("scan", []))

    imagesets_in = d.get("imageset", [])
    beams_in = d.get("beam", [])
    detectors_in = d.get("detector", [])

    new_imagesets = []
    new_beams = []
    new_detectors = []

    for exp in d["experiment"]:
        iset = imagesets_in[exp["imageset"]]
        beam_src = beams_in[exp["beam"]]

        if consolidated:
            sp = exp["scan_point"]
            fi = frame_numbers[sp] - 1
            wl = wavelengths[sp]
        else:
            scan_idx = exp.get("scan")
            if scan_idx is not None and scan_idx in per_frame:
                fi, wl = per_frame[scan_idx]
            else:
                fi = iset.get("single_file_indices", [0])[0]
                wl = None

        template = iset.get("template") or (iset.get("images") or [""])[0]
        new_iset = {
            "__id__": "ImageSet",
            "images": [template],
            "single_file_indices": [fi],
            "mask": iset.get("mask"),
            "gain": iset.get("gain"),
            "pedestal": iset.get("pedestal"),
            "dx": iset.get("dx"),
            "dy": iset.get("dy"),
            "params": iset.get("params", {}),
        }

        new_beam = copy.deepcopy(beam_src)
        new_beam["__id__"] = "monochromatic"
        if wl is not None:
            new_beam["wavelength"] = wl
        elif "wavelength" not in new_beam:
            new_beam["wavelength"] = 0.0

        det_idx = exp.get("detector", 0)
        new_det = copy.deepcopy(detectors_in[det_idx])

        new_idx = len(new_imagesets)
        new_imagesets.append(new_iset)
        new_beams.append(new_beam)
        new_detectors.append(new_det)

        exp["imageset"] = new_idx
        exp["beam"] = new_idx
        exp["detector"] = new_idx
        exp.pop("scan", None)
        exp.pop("scan_point", None)

    d["imageset"] = new_imagesets
    d["beam"] = new_beams
    d["detector"] = new_detectors
    d["scan"] = []

    with open(output_path, "w") as f:
        json.dump(d, f, indent=2)

    logger.info(
        f"Wrote {len(d['experiment'])} experiments as ImageSet format to {output_path}"
    )


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = (
        "dials.expt_sequence_to_set [options] integrated.expt [integrated.refl] "
        "[output.experiments_filename=imageset.expt]"
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

    convert_expt(input_expt, output_expt)

    if params.input.reflections:
        if len(params.input.reflections) > 1:
            sys.exit("Error: provide at most one .refl file as input")
        input_refl = params.input.reflections[0].filename
        output_refl = params.output.reflections_filename
        if input_refl != output_refl:
            shutil.copy(input_refl, output_refl)
            logger.info(
                f"Copied reflection table (no remapping needed) to {output_refl}"
            )


if __name__ == "__main__":
    run()
