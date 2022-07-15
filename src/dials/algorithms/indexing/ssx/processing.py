from __future__ import annotations

import json
import logging
import math
import os
import pathlib
import sys
from dataclasses import dataclass, field
from multiprocessing import Pool
from typing import Any, List, Optional, Tuple, Union

from dxtbx.model import Crystal, Experiment, ExperimentList
from libtbx import Auto, phil

from dials.algorithms.indexing import DialsIndexError
from dials.algorithms.indexing.indexer import Indexer
from dials.algorithms.indexing.max_cell import find_max_cell
from dials.array_family import flex
from dials.command_line.combine_experiments import CombineWithReference

RAD2DEG = 180 / math.pi

logger = logging.getLogger("dials")


@dataclass
class InputToIndex:
    reflection_table: flex.reflection_table = None
    experiment: Experiment = None
    parameters: Any = None
    image_identifier: str = ""
    image_no: int = 0
    method_list: List[str] = field(default_factory=list)
    known_crystal_models: Optional[List[Crystal]] = None


@dataclass
class IndexingResult:
    image: str
    image_no: int
    reflection_table: flex.reflection_table = None
    experiments: ExperimentList = None
    n_strong: int = 0
    n_indexed: List[int] = field(default_factory=list)
    identifiers: List[str] = field(default_factory=list)
    rmsd_x: List[float] = field(default_factory=list)
    rmsd_y: List[float] = field(default_factory=list)
    rmsd_dpsi: List[float] = field(default_factory=list)


loggers_to_disable = [
    "dials.algorithms.refinement.reflection_processor",
    "dials.algorithms.refinement.refiner",
    "dials.algorithms.refinement.reflection_manager",
    "dials.algorithms.indexing.stills_indexer",
    "dials.algorithms.indexing.nave_parameters",
    "dials.algorithms.indexing.basis_vector_search.real_space_grid_search",
    "dials.algorithms.indexing.basis_vector_search.combinations",
    "dials.algorithms.indexing.indexer",
]
debug_loggers_to_disable = [
    "dials.algorithms.indexing.symmetry",
    "dials.algorithms.indexing.lattice_search",
    "dials.algorithms.refinement.engine",
]


def index_one(
    experiment: Experiment,
    reflection_table: flex.reflection_table,
    params: phil.scope_extract,
    method_list: List[str],
    image_no: int,
    known_crystal_models: List[Crystal] = None,
) -> Union[Tuple[ExperimentList, flex.reflection_table], Tuple[bool, bool]]:
    if not reflection_table:
        logger.info(
            f"Image {image_no+1}: Skipped indexing, no strong spots found/remaining."
        )
        return None, None
    # First suppress logging unless in verbose mode.
    if params.individual_log_verbosity < 2:
        for name in loggers_to_disable:
            logging.getLogger(name).disabled = True
    elif params.individual_log_verbosity == 2:
        for name in loggers_to_disable:
            logging.getLogger(name).setLevel(logging.INFO)
        for name in debug_loggers_to_disable:
            logging.getLogger(name).setLevel(logging.INFO)

    # Now try indexing with the chosen methods
    elist = ExperimentList([experiment])
    for method in method_list:
        params.indexing.method = method
        idxr = Indexer.from_parameters(
            reflection_table,
            elist,
            params=params,
            known_crystal_models=known_crystal_models,
        )
        try:
            idxr.index()
        except (DialsIndexError, AssertionError) as e:
            logger.info(
                f"Image {image_no+1}: Failed to index with {method} method, error: {e}"
            )
            if method == method_list[-1]:
                return None, None
        else:
            logger.info(
                f"Image {image_no+1}: Indexed {idxr.refined_reflections.size()}/{reflection_table.size()} spots with {method} method."
            )
            return idxr.refined_experiments, idxr.refined_reflections


def wrap_index_one(input_to_index: InputToIndex) -> IndexingResult:
    # First unpack the input and run the function
    expts, table = index_one(
        input_to_index.experiment,
        input_to_index.reflection_table,
        input_to_index.parameters,
        input_to_index.method_list,
        input_to_index.image_no,
        input_to_index.known_crystal_models,
    )

    # Now calculate some useful quantities for the result
    n_strong = input_to_index.reflection_table.get_flags(
        input_to_index.reflection_table.flags.strong
    ).count(True)
    if expts and table:
        result = IndexingResult(
            input_to_index.image_identifier,
            input_to_index.image_no,
            table,
            expts,
            n_strong=n_strong,
        )
        for id_, identifier in table.experiment_identifiers():
            selr = table.select(table["id"] == id_)
            calx, caly, _ = selr["xyzcal.px"].parts()
            obsx, obsy, _ = selr["xyzobs.px.value"].parts()
            delpsi = selr["delpsical.rad"]
            rmsd_x = flex.mean((calx - obsx) ** 2) ** 0.5
            rmsd_y = flex.mean((caly - obsy) ** 2) ** 0.5
            rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2) ** 0.5
            n_id_ = calx.size()
            result.n_indexed.append(n_id_)
            result.rmsd_x.append(rmsd_x)
            result.rmsd_y.append(rmsd_y)
            result.rmsd_dpsi.append(rmsd_z)
            result.identifiers.append(identifier)
    else:
        result = IndexingResult(
            input_to_index.image_identifier,
            input_to_index.image_no,
            n_strong=n_strong,
        )

    # If chosen, output a message to json to show live progress
    if input_to_index.parameters.output.nuggets:
        msg = {
            "image_no": input_to_index.image_no + 1,
            "n_indexed": (len(table) if table else 0),
            "n_strong": n_strong,
            "n_cryst": (len(expts) if expts else 0),
            "image": input_to_index.image_identifier,
        }
        with open(
            input_to_index.parameters.output.nuggets
            / f"nugget_index_{msg['image']}.json",
            "w",
        ) as f:
            f.write(json.dumps(msg))

    return result


def index_all_concurrent(
    experiments: ExperimentList,
    reflections: List[flex.reflection_table],
    params: phil.scope_extract,
    method_list: List[str],
) -> Tuple[ExperimentList, flex.reflection_table, dict]:

    input_iterable = []
    results_summary = {
        i: [] for i in range(len(experiments))
    }  # create to give results in order

    # Create a suitable iterable for passing to pool.map
    n = 0
    for iset in experiments.imagesets():
        for i in range(len(iset)):
            refl_index = i + n
            if reflections[refl_index]:
                expt = experiments[refl_index]
                idx_0 = i + iset.get_scan().get_batch_offset()  # slicing index
                new_iset = iset[idx_0 : idx_0 + 1]
                expt.imageset = new_iset
                input_iterable.append(
                    InputToIndex(
                        reflection_table=reflections[refl_index],
                        experiment=expt,
                        parameters=params,
                        image_identifier=iset.get_image_identifier(i).split("/")[-1],
                        image_no=refl_index,
                        method_list=method_list,
                    )
                )
            else:  # experiments that have already been filtered
                results_summary[refl_index].append(
                    {
                        "Image": iset.get_image_identifier(i).split("/")[-1],
                        "n_indexed": 0,
                        "n_strong": 0,
                    }
                )
        n += len(iset)

    with open(os.devnull, "w") as devnull:
        sys.stdout = devnull  # block printing from rstbx

        if params.indexing.nproc > 1:
            with Pool(params.indexing.nproc) as pool:
                results: List[IndexingResult] = pool.map(wrap_index_one, input_iterable)
        else:
            results = [wrap_index_one(i) for i in input_iterable]

    sys.stdout = sys.__stdout__
    # prepare tables for output
    indexed_experiments, indexed_reflections = _join_indexing_results(results)

    results_summary = _add_results_to_summary_dict(results_summary, results)

    return indexed_experiments, indexed_reflections, results_summary


def _join_indexing_results(
    results: List[IndexingResult],
) -> Tuple[ExperimentList, flex.reflection_table]:
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()
    n_tot = 0
    for res in results:
        if res.n_indexed:
            indexed_experiments.extend(res.experiments)
            table = res.reflection_table
            ids_map = dict(table.experiment_identifiers())
            for k in table.experiment_identifiers().keys():
                del table.experiment_identifiers()[k]
            table["id"] += n_tot
            for k, v in ids_map.items():
                table.experiment_identifiers()[k + n_tot] = v
            n_tot += len(ids_map.keys())
            indexed_reflections.extend(table)

    indexed_reflections.assert_experiment_identifiers_are_consistent(
        indexed_experiments
    )
    return indexed_experiments, indexed_reflections


def _add_results_to_summary_dict(
    results_summary: dict, results: List[IndexingResult]
) -> dict:
    # convert to results_summary dict current format
    for res in results:
        if res.n_indexed:
            for j in range(len(res.n_indexed)):
                results_summary[res.image_no].append(
                    {
                        "Image": res.image,
                        "identifier": res.identifiers[j],
                        "n_indexed": res.n_indexed[j],
                        "n_strong": res.n_strong,
                        "RMSD_X": res.rmsd_x[j],
                        "RMSD_Y": res.rmsd_y[j],
                        "RMSD_dPsi": res.rmsd_dpsi[j],
                    }
                )
        else:
            results_summary[res.image_no].append(
                {
                    "Image": res.image,
                    "n_indexed": 0,
                    "n_strong": res.n_strong,
                }
            )
    return results_summary


def preprocess(
    experiments: ExperimentList,
    observed: flex.reflection_table,
    params: phil.scope_extract,
) -> Tuple[List[flex.reflection_table], phil.scope_extract, List[str]]:
    reflections = observed.split_by_experiment_id()

    if len(reflections) != len(experiments):
        # spots may not have been found on every image. In this case, the length
        # of the list of reflection tables will be less than the length of experiments.
        # Add in empty items to the list, so that this can be reported on
        obs = set(observed["id"])
        no_refls = set(range(len(experiments))).difference(obs)
        # need to handle both cases where lots have no refls, or only a few do
        full_list = [None] * len(experiments)
        for i, o in enumerate(obs):
            full_list[o] = reflections[i]
        reflections = full_list
        logger.info(f"Filtered {len(no_refls)} images with no spots found.")
        if len(experiments) != len(reflections):
            raise ValueError(
                f"Unequal number of reflection tables {len(reflections)} and experiments {len(experiments)}"
            )

    # Calculate necessary quantities
    n_filtered_out = 0
    filtered_out_images = []
    for i, (refl, experiment) in enumerate(zip(reflections, experiments)):
        if refl:
            if refl.size() >= params.min_spots:
                elist = ExperimentList([experiment])
                refl["imageset_id"] = flex.int(
                    refl.size(), 0
                )  # needed for centroid_px_to_mm
                refl.centroid_px_to_mm(elist)
                refl.map_centroids_to_reciprocal_space(elist)
            else:
                n_filtered_out += 1
                reflections[i] = None
                filtered_out_images.append(i + 1)
        else:
            filtered_out_images.append(i + 1)
    if n_filtered_out:
        logger.info(
            f"Filtered {n_filtered_out} images with fewer than {params.min_spots} spots"
        )
        if params.output.nuggets:
            jstr = json.dumps({"filtered_images": filtered_out_images})
            with open(
                params.output.nuggets / "nugget_index_filtered_images.json", "w"
            ) as f:
                f.write(jstr)

    if n_filtered_out == len(reflections):
        logger.info("All images filtered out, none left to attempt indexing")
        return reflections, params, []

    # Determine the max cell if applicable
    if (params.indexing.max_cell is Auto) and (
        not params.indexing.known_symmetry.unit_cell
    ):
        if params.individual_log_verbosity <= 2:  # suppress the max cell debug log
            logging.getLogger("dials.algorithms.indexing.max_cell").setLevel(
                logging.INFO
            )
        max_cells = []
        for refl in reflections:
            if refl:
                try:
                    max_cells.append(find_max_cell(refl).max_cell)
                except (DialsIndexError, AssertionError):
                    pass
        if not max_cells:
            raise ValueError("Unable to find a max cell for any images")
        logger.info(f"Setting max cell to {max(max_cells):.1f} " + "\u212B")
        params.indexing.max_cell = max(max_cells)

    # Determine which methods to try
    method_list = params.method
    if "real_space_grid_search" in method_list:
        if not params.indexing.known_symmetry.unit_cell:
            logger.info("No unit cell given, real_space_grid_search will not be used")
            method_list.remove("real_space_grid_search")
    methods = ", ".join(method_list)
    pl = "s" if (len(method_list) > 1) else ""
    logger.info(f"Attempting indexing with {methods} method{pl}")

    return reflections, params, method_list


def index(
    experiments: ExperimentList,
    observed: flex.reflection_table,
    params: phil.scope_extract,
) -> Tuple[ExperimentList, flex.reflection_table, dict]:

    if params.output.nuggets:
        params.output.nuggets = pathlib.Path(
            params.output.nuggets
        )  # needed if not going via cli
        if not params.output.nuggets.is_dir():
            logger.warning(
                "output.nuggets not recognised as a valid directory path, no nuggets will be output"
            )
            params.output.nuggets = None

    reflection_tables, params, method_list = preprocess(experiments, observed, params)

    # Do the indexing and generate a summary of results
    indexed_experiments, indexed_reflections, results_summary = index_all_concurrent(
        experiments,
        reflection_tables,
        params,
        method_list,
    )

    # combine beam and detector models if not already
    if (len(indexed_experiments.detectors())) > 1 or (
        len(indexed_experiments.beams())
    ) > 1:
        combine = CombineWithReference(
            detector=indexed_experiments[0].detector, beam=indexed_experiments[0].beam
        )
        elist = ExperimentList()
        for expt in indexed_experiments:
            elist.append(combine(expt))
        indexed_experiments = elist

    return indexed_experiments, indexed_reflections, results_summary
