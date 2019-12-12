from __future__ import absolute_import, division, print_function


class GaussianRSProfileModelExt(object):
    """An extension class implementing a reciprocal space gaussian profile model."""

    name = "gaussian_rs"

    default = True

    @classmethod
    def phil(cls):
        from dials.algorithms.profile_model.gaussian_rs import phil_scope

        return phil_scope

    @classmethod
    def algorithm(cls):
        from dials.algorithms.profile_model.gaussian_rs import Model

        return Model

    @classmethod
    def from_dict(cls, d):
        return cls.algorithm().from_dict(d)
