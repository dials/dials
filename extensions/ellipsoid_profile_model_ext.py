from __future__ import annotations


class EllipsoidProfileModelExt:
    """An extension class implementing a reciprocal space multivariate normal profile model."""

    name = "ellipsoid"

    default = False

    @staticmethod
    def phil():
        from dials.algorithms.profile_model.ellipsoid.model import phil_scope

        return phil_scope

    @staticmethod
    def algorithm():
        from dials.algorithms.profile_model.ellipsoid.model import EllipsoidProfileModel

        return EllipsoidProfileModel

    @classmethod
    def from_dict(cls, d):
        return cls.algorithm().from_dict(d)
