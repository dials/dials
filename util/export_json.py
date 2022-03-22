from __future__ import annotations

import json


class ReciprocalLatticeJson:
    def __init__(self, experiments, reflections):
        self.experiments = experiments
        self.reflections = reflections
        self.reflections.centroid_px_to_mm(experiments)
        self.reflections.map_centroids_to_reciprocal_space(experiments)

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
            with open(filename, "w") as f:
                f.write(text)
        else:
            return text
