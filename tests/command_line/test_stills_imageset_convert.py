"""Tests for the stills ImageSet <-> ImageSequence converters and their shared
consolidated-scan helpers (`dials.util.stills_imageset_convert`).

Guards P2: per-frame identity is stored only as the imageset's 0-based
`single_file_indices` (the redundant 1-based `frame_numbers` array was removed).
The per-frame `image_range`/`fi` reconstruction is implemented twice -- in dxtbx
`ExperimentListDict._expand_consolidated_scans` and in the dials helper
`consolidated_frame_map` (dxtbx cannot import dials) -- and the two must agree.
"""

from __future__ import annotations

import copy
import json

from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.FormatMultiImage import Reader as MultiImageReader
from dxtbx.imageset import ImageSequence, ImageSetData
from dxtbx.model import Beam, Detector, Experiment, ExperimentList, Scan
from dxtbx.model.experiment_list import (
    ExperimentListDict,
    ExperimentListFactory,
)

import dials.command_line.expt_sequence_to_set as seq_to_set
import dials.command_line.expt_set_to_sequence as set_to_seq
from dials.command_line.combine_experiments import (
    combine_experiments_no_reflections,
    phil_scope,
)
from dials.util.stills_imageset_convert import consolidated_frame_map


def _dxtbx_expand(obj):
    """Run dxtbx's `_expand_consolidated_scans` in isolation and return the
    per-experiment 0-based frame index recovered from the expanded scans."""
    eld = ExperimentListDict.__new__(ExperimentListDict)
    eld._obj = copy.deepcopy(obj)
    eld._check_format = False
    eld._directory = None
    eld._expand_consolidated_scans()
    out = []
    for e in eld._obj["experiment"]:
        scan = eld._obj["scan"][e["scan"]]
        out.append(scan["image_range"][0] - 1)  # 1-based image_range -> 0-based fi
    return out


def test_consolidated_reconstruction_dxtbx_dials_agree():
    """Multi-imageset + multi-lattice consolidated dict: the dxtbx and dials
    reconstructions must yield identical per-experiment frame indices.

    scan_point is global and per-unique-frame; single_file_indices is
    per-imageset, sparse, ascending. A frame with >1 lattice means
    len(experiments) > len(single_file_indices)."""
    obj = {
        "__id__": "ExperimentList",
        "imageset": [
            {
                "__id__": "ImageSequence",
                "template": "srcA.h5",
                "single_file_indices": [0, 2, 5],
            },
            {
                "__id__": "ImageSequence",
                "template": "srcB.h5",
                "single_file_indices": [1, 3],
            },
        ],
        "scan": [
            {
                "__stills_consolidated": True,
                "batch_offset": 0,
                # indexed by global scan_point: A0, A2, A5, B1, B3
                "properties": {"wavelength": [1.30, 1.32, 1.35, 1.31, 1.33]},
                "valid_image_ranges": {},
            },
        ],
        "experiment": [
            {"imageset": 0, "scan": 0, "scan_point": 0},  # srcA frame 0, lattice 1
            {"imageset": 0, "scan": 0, "scan_point": 0},  # srcA frame 0, lattice 2
            {"imageset": 0, "scan": 0, "scan_point": 1},  # srcA frame 2
            {"imageset": 0, "scan": 0, "scan_point": 2},  # srcA frame 5
            {"imageset": 1, "scan": 0, "scan_point": 3},  # srcB frame 1
            {"imageset": 1, "scan": 0, "scan_point": 4},  # srcB frame 3
        ],
    }
    expected_fi = [0, 0, 2, 5, 1, 3]

    dxtbx_fi = _dxtbx_expand(obj)

    fmap = consolidated_frame_map(copy.deepcopy(obj))
    dials_fi = [fmap[e["scan_point"]][0] for e in obj["experiment"]]
    dials_wl = [fmap[e["scan_point"]][1] for e in obj["experiment"]]

    assert dxtbx_fi == expected_fi
    assert dials_fi == expected_fi
    assert dxtbx_fi == dials_fi
    # wavelengths follow the per-experiment frame
    assert dials_wl == [1.30, 1.30, 1.32, 1.35, 1.31, 1.33]


def test_consolidated_frame_map_single_imageset():
    """Single-imageset collapse: frame_for_scan_point[sp] == single_file_indices[sp]."""
    obj = {
        "__id__": "ExperimentList",
        "imageset": [
            {
                "__id__": "ImageSequence",
                "template": "src.h5",
                "single_file_indices": [3, 7, 9],
            },
        ],
        "scan": [
            {
                "__stills_consolidated": True,
                "batch_offset": 0,
                "properties": {"wavelength": [1.1, 1.2, 1.3]},
                "valid_image_ranges": {},
            },
        ],
        "experiment": [
            {"imageset": 0, "scan": 0, "scan_point": 0},
            {"imageset": 0, "scan": 0, "scan_point": 1},
            {"imageset": 0, "scan": 0, "scan_point": 2},
        ],
    }
    fmap = consolidated_frame_map(copy.deepcopy(obj))
    assert [fmap[sp][0] for sp in (0, 1, 2)] == [3, 7, 9]
    assert _dxtbx_expand(obj) == [3, 7, 9]


def test_consolidated_frame_map_returns_none_for_nonconsolidated():
    obj = {
        "__id__": "ExperimentList",
        "scan": [{"image_range": [1, 1]}],
        "imageset": [],
        "experiment": [],
    }
    assert consolidated_frame_map(obj) is None
    assert (
        consolidated_frame_map({"scan": [], "imageset": [], "experiment": []}) is None
    )


def _classic_input():
    """A classic per-frame ImageSet .expt dict: multi-source, sparse frames, and
    one frame carrying two lattices. Each experiment is tagged so identity can be
    tracked across the converter round-trip."""
    experiments = [
        # (template, fi, wavelength)
        ("srcA.h5", 0, 1.30),  # frame A0 lattice 1
        ("srcA.h5", 0, 1.30),  # frame A0 lattice 2 (multi-lattice)
        ("srcA.h5", 2, 1.32),  # frame A2
        ("srcB.h5", 1, 1.31),  # frame B1
    ]
    d = {
        "__id__": "ExperimentList",
        "imageset": [],
        "beam": [],
        "detector": [{"panels": [], "_tag": "det"}],
        "crystal": [],
        "experiment": [],
    }
    for tag, (template, fi, wl) in enumerate(experiments):
        d["imageset"].append(
            {"__id__": "ImageSet", "images": [template], "single_file_indices": [fi]}
        )
        d["beam"].append(
            {"__id__": "monochromatic", "wavelength": wl, "direction": [0.0, 0.0, 1.0]}
        )
        d["crystal"].append({"_tag": tag})
        d["experiment"].append(
            {"imageset": tag, "beam": tag, "detector": 0, "crystal": tag, "_tag": tag}
        )
    return d, experiments


def test_converter_roundtrip_multisource_multilattice(tmp_path):
    """set->sequence->set reproduces the original per-frame (src, fi, wl) for each
    experiment, including the multi-lattice frame, and the intermediate sequence
    form is sparse with deduped scan_points."""
    classic, experiments = _classic_input()
    in_path = tmp_path / "classic.expt"
    seq_path = tmp_path / "sequence.expt"
    out_path = tmp_path / "roundtrip.expt"
    with open(in_path, "w") as f:
        json.dump(classic, f)

    # Forward: per-frame ImageSets -> shared ImageSequences + consolidated scan
    set_to_seq.convert_expt(str(in_path), str(seq_path))
    seq = json.loads(seq_path.read_text())

    # Sparse single_file_indices, one ImageSequence per source file
    sfi_by_template = {
        iset["template"]: iset["single_file_indices"] for iset in seq["imageset"]
    }
    assert sfi_by_template == {"srcA.h5": [0, 2], "srcB.h5": [1]}
    # frame_numbers gone; one wavelength per unique frame (3, not 4 experiments)
    con = seq["scan"][0]
    assert con.get("__stills_consolidated") is True
    assert "frame_numbers" not in con
    assert len(con["properties"]["wavelength"]) == 3
    # The two lattices on srcA frame 0 share one scan_point.
    sps = [e["scan_point"] for e in seq["experiment"]]
    assert len(set(sps)) == 3 and len(sps) == 4

    # Backward: -> per-frame ImageSets again
    seq_to_set.convert_expt(str(seq_path), str(out_path))
    back = json.loads(out_path.read_text())

    recovered = {}
    for e in back["experiment"]:
        iset = back["imageset"][e["imageset"]]
        beam = back["beam"][e["beam"]]
        recovered[e["_tag"]] = (
            iset["images"][0],
            iset["single_file_indices"][0],
            beam["wavelength"],
        )
    expected = {tag: experiments[tag] for tag in range(len(experiments))}
    assert recovered == expected


def _single_file_still_sequence(n):
    """A single-file-reader ImageSequence over n frames (no real data)."""
    reader = MultiImageReader(FormatMultiImage, ["dummy.h5"], num_images=n)
    return ImageSequence(
        ImageSetData(reader=reader, masker=None),
        scan=Scan(image_range=(1, n), oscillation=(0.0, 0.0)),
        beam=Beam(),
        detector=Detector(),
        goniometer=None,
    )


def test_split_combine_multilattice_roundtrip(tmp_path):
    """split_experiments -> combine_experiments on a multi-lattice composite.

    Regression for the post-combine case the handoff missed: combine preserves a
    separate Scan object per experiment (Scan identity is broken by serializing
    each split experiment to its own file), so multi-lattice experiments on the
    same frame would yield more scan_points than the imageset's
    single_file_indices. The dxtbx write side assigns scan_point per
    (imageset, unique frame), so the consolidated wavelength array stays
    one-entry-per-frame and the file decodes. The old frame_numbers array carried
    one entry per scan and never needed this; P2 does.
    """
    n = 3
    shared = _single_file_still_sequence(n)
    frames = [0, 0, 1, 2]  # frame 0 carries two lattices
    wls = {0: 1.30, 1: 1.31, 2: 1.32}

    composite = ExperimentList()
    for k, fi in enumerate(frames):
        scan = Scan((fi + 1, fi + 1), (0.0, 0.0))
        scan.set_property("wavelength", [wls[fi]])
        beam = Beam()
        beam.set_wavelength(wls[fi])
        expt = Experiment(
            beam=beam, detector=shared.get_detector(), scan=scan, imageset=shared
        )
        expt.identifier = "id%d" % k
        composite.append(expt)

    # split: serialize each experiment to its own file, then reload (this is what
    # breaks Scan identity for the two lattices on frame 0).
    reloaded = []
    for k, expt in enumerate(composite):
        el = ExperimentList()
        el.append(expt)
        p = tmp_path / ("split_%d.expt" % k)
        el.as_json(str(p))
        reloaded.append(
            ExperimentListFactory.from_json_file(str(p), check_format=False)
        )

    params = phil_scope.extract()
    params.reference_from_experiment.beam = 0
    params.reference_from_experiment.detector = 0
    params.reference_from_experiment.scan = 0
    combined = combine_experiments_no_reflections(params, reloaded)

    obj = combined.to_dict()
    con = next(s for s in obj["scan"] if s.get("__stills_consolidated"))
    assert "frame_numbers" not in con
    # one wavelength per unique frame (3), not per experiment (4)
    assert len(con["properties"]["wavelength"]) == 3
    # the two frame-0 lattices share one scan_point
    sps = [e["scan_point"] for e in obj["experiment"]]
    assert sps == [0, 0, 1, 2]

    restored = ExperimentListDict(obj, check_format=False).decode()
    assert [e.scan.get_image_range() for e in restored] == [
        (1, 1),
        (1, 1),
        (2, 2),
        (3, 3),
    ]
    assert [round(e.beam.get_wavelength(), 4) for e in restored] == [
        1.30,
        1.30,
        1.31,
        1.32,
    ]
    assert list(restored.identifiers()) == ["id0", "id1", "id2", "id3"]
