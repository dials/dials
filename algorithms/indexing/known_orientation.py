from __future__ import annotations

from dxtbx.model.experiment_list import Experiment, ExperimentList

from dials.algorithms.indexing.indexer import Indexer


class IndexerKnownOrientation(Indexer):
    def __init__(self, reflections, experiments, params, known_orientations):
        self.known_orientations = known_orientations
        super().__init__(reflections, experiments, params)

    def find_lattices(self):
        experiments = ExperimentList()
        for cm, expt in zip(self.known_orientations, self.experiments):
            # indexer expects crystals to be in primitive setting
            space_group = cm.get_space_group()
            cb_op_to_primitive = (
                space_group.info().change_of_basis_op_to_primitive_setting()
            )
            cm = cm.change_basis(cb_op_to_primitive)
            experiments.append(
                Experiment(
                    imageset=expt.imageset,
                    beam=expt.beam,
                    detector=expt.detector,
                    goniometer=expt.goniometer,
                    scan=expt.scan,
                    crystal=cm,
                    identifier=expt.identifier,
                )
            )
        return experiments
