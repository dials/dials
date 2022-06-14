from __future__ import annotations

import concurrent.futures
import json
import logging
import math
import os
import pathlib
import sys
from collections import defaultdict
from typing import List, Tuple, Union

from dxtbx.model import Experiment, ExperimentList
from libtbx import Auto, phil

from dials.algorithms.indexing import DialsIndexError
from dials.algorithms.indexing.indexer import Indexer
from dials.algorithms.indexing.max_cell import find_max_cell
from dials.array_family import flex
from dials.command_line.combine_experiments import CombineWithReference

RAD2DEG = 180 / math.pi

logger = logging.getLogger("dials")

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
        idxr = Indexer.from_parameters(reflection_table, elist, params=params)
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


def index_all_concurrent(
    experiments: ExperimentList,
    reflections: List[flex.reflection_table],
    params: phil.scope_extract,
    method_list: List[str],
) -> Tuple[ExperimentList, flex.reflection_table, dict]:

    # first determine n_strong per image:
    n_strong_per_image = {}
    isets = experiments.imagesets()
    iset_to_expt_indices = {}
    for iset in isets:
        expt_indices = [
            i for i, expt in enumerate(experiments) if expt.imageset == iset
        ]
        iset_to_expt_indices[iset] = expt_indices
        for i, k in enumerate(expt_indices):
            table = reflections[k]
            img = iset.get_image_identifier(i).split("/")[-1]
            if table:
                n_strong = table.get_flags(table.flags.strong).count(True)
            else:
                n_strong = 0
            n_strong_per_image[img] = n_strong

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=params.indexing.nproc
    ) as pool, open(os.devnull, "w") as devnull:
        sys.stdout = devnull  # block printing from rstbx
        tables_list = [None] * len(reflections)
        expts_list = [None] * len(reflections)
        for isetno, (iset, expt_indices) in enumerate(iset_to_expt_indices.items()):
            if len(isets) > 1:
                logger.info(f"Indexing imageset {isetno}")
            futures = {}
            for i, k in enumerate(expt_indices):
                expt = experiments[k]
                table = reflections[k]
                if table:
                    futures[
                        pool.submit(index_one, expt, table, params, method_list, i)
                    ] = (i, k)
            for future in concurrent.futures.as_completed(futures):
                (j, l) = futures[future]
                if params.output.nuggets:
                    img = iset.get_image_identifier(j).split("/")[-1]
                    msg = {
                        "image_no": j + 1,
                        "n_indexed": 0,
                        "n_strong": n_strong_per_image[img],
                        "n_cryst": 0,
                        "image": img,
                    }
                try:
                    expts, refls = future.result()
                except Exception as e:
                    logger.info(e)
                else:
                    if refls and expts:
                        tables_list[l] = refls
                        elist = ExperimentList()
                        for jexpt in expts:
                            elist.append(
                                Experiment(
                                    identifier=jexpt.identifier,
                                    beam=jexpt.beam,
                                    detector=jexpt.detector,
                                    scan=jexpt.scan,
                                    goniometer=jexpt.goniometer,
                                    crystal=jexpt.crystal,
                                    imageset=jexpt.imageset[j : j + 1],
                                )
                            )
                        expts_list[l] = elist
                        if params.output.nuggets:
                            msg["n_indexed"] = len(refls)
                            msg["n_cryst"] = len(elist)
                if params.output.nuggets:
                    with open(
                        params.output.nuggets / f"nugget_index_{msg['image']}.json", "w"
                    ) as f:
                        f.write(json.dumps(msg))

    sys.stdout = sys.__stdout__  # restore printing

    # now postprocess - record generic information
    results_summary = defaultdict(list)
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()
    n_tot = 0
    for iset, expt_indices in iset_to_expt_indices.items():  # loop over imagesets
        for i, j in enumerate(expt_indices):
            elist = expts_list[j]
            table = tables_list[j]
            img = iset.get_image_identifier(i).split("/")[-1]
            n_strong = n_strong_per_image[img]
            if not (elist and table):
                results_summary[j].append(
                    {
                        "Image": img,
                        "n_indexed": 0,
                        "n_strong": n_strong,
                    }
                )
                continue
            indexed_experiments.extend(elist)
            ids_map = dict(table.experiment_identifiers())
            for k in table.experiment_identifiers().keys():
                del table.experiment_identifiers()[k]
            table["id"] += n_tot
            for k, v in ids_map.items():
                table.experiment_identifiers()[k + n_tot] = v
            n_tot += len(ids_map.keys())
            indexed_reflections.extend(table)

            # record some things for printing to output log/html
            for id_, identifier in table.experiment_identifiers():
                selr = table.select(table["id"] == id_)
                calx, caly, _ = selr["xyzcal.px"].parts()
                obsx, obsy, _ = selr["xyzobs.px.value"].parts()
                delpsi = selr["delpsical.rad"]
                rmsd_x = flex.mean((calx - obsx) ** 2) ** 0.5
                rmsd_y = flex.mean((caly - obsy) ** 2) ** 0.5
                rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2) ** 0.5
                n_id_ = calx.size()
                results_summary[j].append(
                    {
                        "Image": img,
                        "identifier": identifier,
                        "n_indexed": n_id_,
                        "n_strong": n_strong,
                        "RMSD_X": rmsd_x,
                        "RMSD_Y": rmsd_y,
                        "RMSD_dPsi": rmsd_z,
                    }
                )

    indexed_reflections.assert_experiment_identifiers_are_consistent(
        indexed_experiments
    )
    return indexed_experiments, indexed_reflections, results_summary


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
        no_refls = set(range(len(experiments))).difference(set(observed["id"]))
        for i in no_refls:
            reflections.insert(i, None)
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
