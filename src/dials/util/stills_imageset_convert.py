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


def build_consolidated_scan(frame_numbers, wavelengths, batch_offset=0):
    """Return a __stills_consolidated scan dict.

    Args:
      frame_numbers: list of 1-based image-range starts, length N (one per experiment).
      wavelengths:   list of per-frame wavelengths in Angstrom, length N.
      batch_offset:  integer batch offset (default 0).

    The returned dict is the only acceptable JSON form for a consolidated scan
    and matches what :func:`expand_consolidated_scan` knows how to read back.
    """
    return {
        "__stills_consolidated": True,
        "batch_offset": batch_offset,
        "frame_numbers": list(frame_numbers),
        "properties": {"wavelength": list(wavelengths)},
        "valid_image_ranges": {},
    }


def expand_consolidated_scan(scan_list):
    """Return ``(frame_numbers, wavelengths)`` from a __stills_consolidated scan list.

    Returns ``(None, None)`` if the scan list is absent or not consolidated, so
    callers can detect that case and fall back to per-frame scan dicts.
    """
    if not scan_list or not scan_list[0].get("__stills_consolidated"):
        return None, None
    con = scan_list[0]
    fns = con["frame_numbers"]
    wls = con.get("properties", {}).get("wavelength", [None] * len(fns))
    return fns, wls


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
