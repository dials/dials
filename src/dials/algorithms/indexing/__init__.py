from __future__ import annotations

import functools

import dials_algorithms_indexing_ext as ext

sampling_volume_map = ext.sampling_volume_map


@functools.cache
def rstbx_indexing_api_phil_scope():
    """Parse the rstbx ``indexing_api_defs`` PHIL defaults once and cache the
    parsed scope.

    Parsing this string is relatively expensive and was previously repeated
    every time a StillsIndexer or FFT1D strategy was constructed (i.e. once per
    still). Callers should ``.extract()`` the returned scope to obtain a fresh,
    independently mutable object.
    """
    import iotbx.phil
    from rstbx.phil.phil_preferences import indexing_api_defs

    return iotbx.phil.parse(input_string=indexing_api_defs)


class DialsIndexError(RuntimeError):
    pass


class DialsIndexRefineError(DialsIndexError):
    pass
