from __future__ import annotations


class GaussianRSProfileModelExt:
    """An extension class implementing a reciprocal space gaussian profile model."""

    name = "gaussian_rs"

    default = True

    @staticmethod
    def phil():
        from dials.algorithms.profile_model.gaussian_rs import phil_scope

        return phil_scope

    @staticmethod
    def algorithm():
        from dials.algorithms.profile_model.gaussian_rs import Model

        return Model

    @classmethod
    def from_dict(cls, d):
        return cls.algorithm().from_dict(d)
