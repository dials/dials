from __future__ import annotations

# The scaling model module loads all registered scaling models (including
# tof.model.TOFScalingModel) when it is imported, so it must be imported before
# tof.model, or importing tof.model first would be circular.
import dials.algorithms.scaling.model.model  # noqa: F401
