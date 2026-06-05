"""Shared helpers for the dials.expt_set_to_sequence / expt_sequence_to_set scripts.

These two command-line tools convert .expt files between the classic per-frame
ImageSet representation and the stills-as-ImageSequence representation produced
by the stills_process_imagesequence branch.  Both tools manipulate the .expt
JSON directly (rather than round-tripping through ExperimentList), so they share
the same low-level helpers for reading and writing the __stills_consolidated
scan dict, the imageset shape, and the per-frame frame-index lookup.

Kept in one place so the conversion logic does not drift between the two
scripts (a divergent definition of "what does a consolidated scan look like?"
would silently corrupt round-trips).
"""

from __future__ import annotations


def source_file(iset_dict):
    """Return the source-file path of an imageset dict, or "" if absent.

    ImageSequence dicts use ``template``; classic ImageSet dicts use the first
    entry of ``images``.
    """
    if iset_dict.get("template"):
        return iset_dict["template"]
    images = iset_dict.get("images", [])
    return images[0] if images else ""


def frame_index(iset_dict):
    """Return the 0-based HDF5 frame index of a per-frame imageset dict.

    Used to sort experiments by frame in the set→sequence direction.
    """
    sfi = iset_dict.get("single_file_indices", [])
    return sfi[0] if sfi else 0


def build_consolidated_scan(single_file_indices, wavelengths, batch_offset=0):
    """Return a __stills_consolidated scan dict.

    Per-frame identity is stored only as 0-based ``single_file_indices`` on the
    imageset (the redundant 1-based ``frame_numbers`` array was removed); the
    consolidated scan carries just the per-frame property arrays. This function
    builds the scan and enforces the alignment invariant that the wavelength
    array has exactly one entry per unique frame.

    Args:
      single_file_indices: the imageset's 0-based frame indices (sparse,
        ascending, one per unique frame). Used only to validate length.
      wavelengths:   list of per-unique-frame wavelengths in Angstrom.
      batch_offset:  integer batch offset (default 0).
    """
    assert len(single_file_indices) == len(wavelengths), (
        "consolidated stills: single_file_indices length %d != wavelength "
        "length %d" % (len(single_file_indices), len(wavelengths))
    )
    return {
        "__stills_consolidated": True,
        "batch_offset": batch_offset,
        "properties": {"wavelength": list(wavelengths)},
        "valid_image_ranges": {},
    }


def consolidated_frame_map(d):
    """Return ``{scan_point: (fi, wl)}`` for a consolidated-scan .expt dict.

    ``fi`` is the 0-based HDF5 frame index recovered from the referenced
    imageset's ``single_file_indices``; ``wl`` is the per-frame wavelength.
    Returns ``None`` if the scan is absent or not consolidated, so callers can
    fall back to :func:`per_frame_scan_map`.

    This is the dials-side mirror of dxtbx's
    ``ExperimentListDict._expand_consolidated_scans`` (dxtbx cannot import
    dials). The two must stay identical; a cross-repo unit test guards them.

    Invariant: within each imageset, ``single_file_indices`` is sparse and
    ascending and aligns positionally with that imageset's experiments taken in
    ascending (deduped) scan_point order. Multi-lattice frames share one
    scan_point, so dedup before zipping; group per imageset since scan_point is
    global while single_file_indices is per-imageset.
    """
    scan_list = d.get("scan", [])
    if not scan_list or not scan_list[0].get("__stills_consolidated"):
        return None
    wls = scan_list[0].get("properties", {}).get("wavelength", [])
    imagesets = d.get("imageset", [])

    sps_by_imageset = {}
    for exp in d["experiment"]:
        sp = exp.get("scan_point")
        if sp is None:
            continue
        sps_by_imageset.setdefault(exp["imageset"], set()).add(sp)

    frame_for_scan_point = {}
    for m, sp_set in sps_by_imageset.items():
        sps = sorted(sp_set)
        sfi_m = imagesets[m]["single_file_indices"]
        assert len(sfi_m) == len(sps), (
            "consolidated stills: single_file_indices length %d != %d unique "
            "scan_points for imageset %d" % (len(sfi_m), len(sps), m)
        )
        for local_pos, sp in enumerate(sps):
            frame_for_scan_point[sp] = sfi_m[local_pos]

    return {
        sp: (fi, wls[sp] if sp < len(wls) else None)
        for sp, fi in frame_for_scan_point.items()
    }


def per_frame_scan_map(scan_list):
    """Build ``{scan_index: (fi, wl)}`` from non-consolidated per-frame scan dicts.

    ``fi`` is the 0-based HDF5 frame index (``image_range[0] - 1``); ``wl`` is
    the first wavelength entry from the scan properties, or ``None`` if absent.
    """
    result = {}
    for i, s in enumerate(scan_list):
        fi = s["image_range"][0] - 1
        wls = s.get("properties", {}).get("wavelength", [])
        result[i] = (fi, wls[0] if wls else None)
    return result
