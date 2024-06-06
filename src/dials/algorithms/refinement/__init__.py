from __future__ import annotations

from dials.algorithms.refinement.refiner import Refiner, RefinerFactory


class DialsRefineConfigError(ValueError):
    pass


class DialsRefineRuntimeError(RuntimeError):
    pass


__all__ = [
    "DialsRefineConfigError",
    "DialsRefineRuntimeError",
    "Refiner",
    "RefinerFactory",
]
