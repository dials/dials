from __future__ import annotations

from libtbx.str_utils import wordwrap


class ClusterInfo:
    def __init__(
        self, cluster_id, labels, multiplicity, completeness, unit_cell, height=None
    ):
        self.cluster_id = cluster_id
        self.labels = labels
        self.multiplicity = multiplicity
        self.completeness = completeness
        self.unit_cell = unit_cell
        self.height = height

    def __str__(self):
        lines = [
            "Cluster %i" % self.cluster_id,
            "  Number of datasets: %i" % len(self.labels),
            "  Completeness: %.1f %%" % (self.completeness * 100),
            "  Multiplicity: %.2f" % self.multiplicity,
            "  Datasets:" + ",".join("%s" % s for s in self.labels),
        ]
        if self.height is not None:
            lines.append("  height: %f" % self.height)
        return "\n".join(lines)

    @staticmethod
    def as_table(cluster_info):
        headers = [
            "Cluster",
            "No. datasets",
            "Datasets",
            "Height",
            "Multiplicity",
            "Completeness",
        ]
        rows = []
        for info in cluster_info:
            rows.append(
                [
                    "%i" % info.cluster_id,
                    "%i" % len(info.labels),
                    wordwrap(" ".join("%s" % l for l in info.labels)),
                    "%.2g" % info.height,
                    "%.1f" % info.multiplicity,
                    "%.2f" % info.completeness,
                ]
            )

        rows.insert(0, headers)
        return rows
