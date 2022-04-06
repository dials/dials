from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def refine(params, reflections, experiments):
    if params.refinement.parameterisation.scan_varying:
        logger.warning(
            "scan_varying=True not supported in indexing: setting scan_varying=False"
        )
        params.refinement.parameterisation.scan_varying = False

    from dials.algorithms.refinement import RefinerFactory

    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments
    )

    outliers = None
    refined = refiner.run()
    return refiner, refined, outliers
