#!/usr/bin/env python
#
# export_json.py
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division, print_function

import json

from dials.array_family import flex


class ReciprocalLatticeJson(object):
    def __init__(self, experiments, reflections):
        self.experiments = experiments
        self.reflections = self._map_centroids_to_reciprocal_space(
            experiments, reflections
        )

    def _map_centroids_to_reciprocal_space(self, experiments, reflections):
        refl = reflections
        reflections = flex.reflection_table()
        for i, expt in enumerate(self.experiments):
            refl_sel = refl.select(refl["imageset_id"] == i)
            refl_sel.centroid_px_to_mm(expt.detector, expt.scan)
            refl_sel.map_centroids_to_reciprocal_space(
                expt.detector, expt.beam, expt.goniometer
            )
            reflections.extend(refl_sel)
        return reflections

    def as_dict(self, n_digits=None):
        rlp = self.reflections["rlp"]
        if n_digits is not None:
            rlp = rlp.round(n_digits)
        flat_rlp = []
        for r in rlp:
            flat_rlp.extend(r)

        if "imageset_id" in self.reflections:
            imageset_id = list(self.reflections["imageset_id"])
            expt_id = list(self.reflections["id"])
        else:
            imageset_id = list(self.reflections["id"])
            expt_id = None

        d = {"rlp": flat_rlp, "imageset_id": imageset_id, "experiment_id": expt_id}
        return d

    def as_json(self, filename=None, compact=False, n_digits=None, experiments=None):
        d = self.as_dict(n_digits=n_digits)
        if experiments:
            d["experiments"] = experiments.to_dict()
        if compact:
            text = json.dumps(d, separators=(",", ":"), ensure_ascii=True)
        else:
            text = json.dumps(d, separators=(",", ": "), indent=1, ensure_ascii=True)
        if filename is not None:
            with open(filename, "wb") as f:
                f.write(text)
        else:
            return text
