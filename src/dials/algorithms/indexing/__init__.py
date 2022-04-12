from __future__ import annotations

import dials_algorithms_indexing_ext as ext

sampling_volume_map = ext.sampling_volume_map


class DialsIndexError(RuntimeError):
    pass


class DialsIndexRefineError(DialsIndexError):
    pass
