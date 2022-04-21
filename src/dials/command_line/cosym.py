from __future__ import annotations

import logging
import random
import sys

import numpy as np

import iotbx.phil
from cctbx import sgtbx

from dials.algorithms.clustering.unit_cell import cluster_unit_cells
from dials.algorithms.symmetry.cosym import CosymAnalysis
from dials.algorithms.symmetry.cosym.observers import register_default_cosym_observers
from dials.array_family import flex
from dials.command_line.symmetry import (
    apply_change_of_basis_ops,
    change_of_basis_ops_to_minimum_cell,
    eliminate_sys_absent,
)
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.exclude_images import get_selection_for_valid_image_ranges
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
    select_datasets_on_identifiers,
    update_imageset_ids,
)
from dials.util.observer import Subject
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.cosym")

phil_scope = iotbx.phil.parse(
    """\
partiality_threshold = 0.4
  .type = float
  .help = "Use reflections with a partiality above the threshold."

unit_cell_clustering {
  threshold = 5000
    .type = float(value_min=0, allow_none=True)
    .help = 'Threshold value for the clustering'
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}

include scope dials.algorithms.symmetry.cosym.phil_scope

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

min_reflections = 10
  .type = int(value_min=1)
  .help = "The minimum number of reflections per experiment."

seed = 230
  .type = int(value_min=0)

output {
  suffix = "_reindexed"
    .type = str
  log = dials.cosym.log
    .type = str
  experiments = "symmetrized.expt"
    .type = path
  reflections = "symmetrized.refl"
    .type = path
  json = dials.cosym.json
    .type = path
  html = dials.cosym.html
    .type = path
}
""",
    process_includes=True,
)


class cosym(Subject):
    def __init__(self, experiments, reflections, params=None):
        super().__init__(events=["run_cosym", "performed_unit_cell_clustering"])
        if params is None:
            params = phil_scope.extract()
        self.params = params

        self._reflections = []
        for refl, expt in zip(reflections, experiments):
            sel = get_selection_for_valid_image_ranges(refl, expt)
            self._reflections.append(refl.select(sel))

        self._experiments, self._reflections = self._filter_min_reflections(
            experiments, self._reflections
        )
        self.ids_to_identifiers_map = {}
        for table in self._reflections:
            self.ids_to_identifiers_map.update(table.experiment_identifiers())
        self.identifiers_to_ids_map = {
            value: key for key, value in self.ids_to_identifiers_map.items()
        }

        if len(self._experiments) > 1 and self.params.unit_cell_clustering.threshold:
            # perform unit cell clustering
            identifiers = self._unit_cell_clustering(self._experiments)
            if len(identifiers) < len(self._experiments):
                logger.info(
                    "Selecting subset of %i datasets for cosym analysis: %s",
                    len(identifiers),
                    str(identifiers),
                )
                self._experiments, self._reflections = select_datasets_on_identifiers(
                    self._experiments, self._reflections, use_datasets=identifiers
                )

        # Map experiments and reflections to minimum cell
        cb_ops = change_of_basis_ops_to_minimum_cell(
            self._experiments,
            params.lattice_symmetry_max_delta,
            params.relative_length_tolerance,
            params.absolute_angle_tolerance,
        )
        exclude = [
            expt.identifier
            for expt, cb_op in zip(self._experiments, cb_ops)
            if not cb_op
        ]
        if len(exclude):
            logger.info(
                f"Rejecting {len(exclude)} datasets from cosym analysis "
                f"(couldn't determine consistent cb_op to minimum cell):\n"
                f"{exclude}",
            )
            self._experiments, self._reflections = select_datasets_on_identifiers(
                self._experiments, self._reflections, exclude_datasets=exclude
            )
            cb_ops = list(filter(None, cb_ops))

        # Eliminate reflections that are systematically absent due to centring
        # of the lattice, otherwise they would lead to non-integer miller indices
        # when reindexing to a primitive setting
        self._reflections = eliminate_sys_absent(self._experiments, self._reflections)

        self._experiments, self._reflections = apply_change_of_basis_ops(
            self._experiments, self._reflections, cb_ops
        )

        # transform models into miller arrays
        datasets = filtered_arrays_from_experiments_reflections(
            self.experiments,
            self.reflections,
            outlier_rejection_after_filter=False,
            partiality_threshold=params.partiality_threshold,
        )

        datasets = [
            ma.as_non_anomalous_array().merge_equivalents().array() for ma in datasets
        ]
        self.cosym_analysis = CosymAnalysis(datasets, self.params)

    @property
    def experiments(self):
        """Return the experiment list."""
        return self._experiments

    @property
    def reflections(self):
        """Return the list of reflection tables."""
        return self._reflections

    @Subject.notify_event(event="run_cosym")
    def run(self):
        self.cosym_analysis.run()
        reindexing_ops = self.cosym_analysis.reindexing_ops
        datasets_ = list(set(self.cosym_analysis.dataset_ids))

        # Log reindexing operators
        logger.info("Reindexing operators:")
        for cb_op in set(reindexing_ops):
            datasets = [d for d, o in zip(datasets_, reindexing_ops) if o == cb_op]
            logger.info(f"{cb_op}: {datasets}")

        self._apply_reindexing_operators(
            reindexing_ops, subgroup=self.cosym_analysis.best_subgroup
        )

    def export(self):
        """Output the datafiles for cosym.

        This includes the cosym.json, reflections and experiments files."""

        reindexed_reflections = flex.reflection_table()
        self._reflections = update_imageset_ids(self._experiments, self._reflections)
        for refl in self._reflections:
            reindexed_reflections.extend(refl)
        reindexed_reflections.reset_ids()

        logger.info(
            "Saving reindexed experiments to %s", self.params.output.experiments
        )
        self._experiments.as_file(self.params.output.experiments)
        logger.info(
            "Saving reindexed reflections to %s", self.params.output.reflections
        )
        reindexed_reflections.as_file(self.params.output.reflections)

    def _apply_reindexing_operators(self, reindexing_ops, subgroup=None):
        """Apply the reindexing operators to the reflections and experiments."""
        unique_ids = set(self.cosym_analysis.dataset_ids)
        for cb_op, dataset_id in zip(reindexing_ops, unique_ids):
            cb_op = sgtbx.change_of_basis_op(cb_op)
            logger.debug(
                "Applying reindexing op %s to dataset %i", cb_op.as_xyz(), dataset_id
            )
            expt = self._experiments[dataset_id]
            refl = self._reflections[dataset_id]
            expt.crystal = expt.crystal.change_basis(cb_op)
            if subgroup is not None:
                cb_op = subgroup["cb_op_inp_best"] * cb_op
                expt.crystal = expt.crystal.change_basis(cb_op)
                expt.crystal.set_space_group(
                    subgroup["best_subsym"].space_group().build_derived_acentric_group()
                )
            else:
                expt.crystal = expt.crystal.change_basis(cb_op)
            expt.crystal.set_unit_cell(
                expt.crystal.get_space_group().average_unit_cell(
                    expt.crystal.get_unit_cell()
                )
            )
            refl["miller_index"] = cb_op.apply(refl["miller_index"])
        # Allow for the case where some datasets are filtered out.
        if len(reindexing_ops) < len(self._experiments):
            to_delete = [
                i for i in range(len(self._experiments)) if i not in unique_ids
            ]
            for idx in sorted(to_delete, reverse=True):
                logger.info(
                    f"Removing dataset {idx} as unable to determine reindexing operator"
                )
                del self._experiments[idx]
                del self._reflections[idx]

    def _filter_min_reflections(self, experiments, reflections):
        identifiers = []

        for expt, refl in zip(experiments, reflections):
            if len(refl) >= self.params.min_reflections:
                identifiers.append(expt.identifier)

        return select_datasets_on_identifiers(
            experiments, reflections, use_datasets=identifiers
        )

    @Subject.notify_event("performed_unit_cell_clustering")
    def _unit_cell_clustering(self, experiments):
        crystal_symmetries = [
            expt.crystal.get_crystal_symmetry() for expt in experiments
        ]
        # lattice ids used to label plots, so want numerical ids
        lattice_ids = [
            self.identifiers_to_ids_map[i] for i in experiments.identifiers()
        ]

        clustering = cluster_unit_cells(
            crystal_symmetries,
            lattice_ids=lattice_ids,
        )
        self.unit_cell_clusters = clustering.clusters
        self.unit_cell_dendrogram = clustering.dendrogram

        logger.info(clustering)
        largest_cluster_lattice_ids = None
        for cluster in self.unit_cell_clusters:
            cluster_lattice_ids = cluster.lattice_ids
            if largest_cluster_lattice_ids is None:
                largest_cluster_lattice_ids = cluster_lattice_ids
            elif len(cluster_lattice_ids) > len(largest_cluster_lattice_ids):
                largest_cluster_lattice_ids = cluster_lattice_ids

        dataset_selection = largest_cluster_lattice_ids
        # now convert to actual identifiers for selection
        return [self.ids_to_identifiers_map[i] for i in dataset_selection]


help_message = """
This program implements the methods of `Gildea, R. J. & Winter, G. (2018).
Acta Cryst. D74, 405-410 <https://doi.org/10.1107/S2059798318002978>`_ for
determination of Patterson group symmetry from sparse multi-crystal data sets in
the presence of an indexing ambiguity.

The program takes as input a set of integrated experiments and reflections,
either in one file per experiment, or with all experiments combined in a single
models.expt and observations.refl file. It will perform analysis of the
symmetry elements present in the datasets and, if necessary, reindex experiments
and reflections as necessary to ensure that all output experiments and
reflections are indexed consistently.

Examples::

  dials.cosym models.expt observations.refl

  dials.cosym models.expt observations.refl space_group=I23

  dials.cosym models.expt observations.refl space_group=I23 lattice_group=I23
"""


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.cosym [options] models.expt observations.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        flex.set_random_seed(params.seed)
        np.random.seed(params.seed)
        random.seed(params.seed)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)
    if len(experiments) != len(reflections):
        raise Sorry(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        cosym_instance = cosym(
            experiments=experiments, reflections=reflections, params=params
        )
    except ValueError as e:
        raise Sorry(e)

    if params.output.html or params.output.json:
        register_default_cosym_observers(cosym_instance)
    cosym_instance.run()
    cosym_instance.export()


if __name__ == "__main__":
    run()
