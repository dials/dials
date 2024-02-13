# Do not import this file directly. Use
#
#   from dials.array_family import flex
#

from __future__ import annotations

import collections
import copy
import functools
import itertools
import logging
import operator
import os
import pickle
from typing import List, Tuple

import numpy as np
import pandas as pd
from annlib_ext import AnnAdaptorSelfInclude

import boost_adaptbx.boost.python
import cctbx.array_family.flex
import cctbx.miller
import libtbx.smart_open
from scitbx import matrix

import dials.extensions.glm_background_ext
import dials.extensions.simple_centroid_ext
import dials.util.ext
import dials_array_family_flex_ext
from dials.algorithms.centroid import centroid_px_to_mm_panel
from dials.util.exclude_images import expand_exclude_multiples, set_invalid_images

__all__ = ["real", "reflection_table_selector"]

logger = logging.getLogger(__name__)

# Set the 'real' type to either float or double
if dials_array_family_flex_ext.get_real_type() == "float":
    real = cctbx.array_family.flex.float
elif dials_array_family_flex_ext.get_real_type() == "double":
    real = cctbx.array_family.flex.double
else:
    raise TypeError('unknown "real" type')


@boost_adaptbx.boost.python.inject_into(dials_array_family_flex_ext.reflection_table)
class _:
    """
    An injector class to add additional methods to the reflection table.
    """

    # Set the default algorithms. These are set as class variables so that if they
    # are changed in the class, all new instances of reflection table will have
    # the modified algorithms. If these are modified on the instance level, then
    # only the instance will have the modified algorithms and new instances will
    # have the defaults
    background_algorithm = functools.partial(
        dials.extensions.glm_background_ext.GLMBackgroundExt, None
    )
    centroid_algorithm = functools.partial(
        dials.extensions.simple_centroid_ext.SimpleCentroidExt, None
    )

    @staticmethod
    def from_predictions(
        experiment, dmin=None, dmax=None, margin=1, force_static=False, padding=0
    ):
        """
        Construct a reflection table from predictions.

        :param experiment: The experiment to predict from
        :param dmin: The maximum resolution
        :param dmax: The minimum resolution
        :param margin: The margin to predict around
        :param force_static: Do static prediction with a scan varying model
        :param padding: Padding in degrees
        :return: The reflection table of predictions
        """
        if experiment.profile is not None:
            return experiment.profile.predict_reflections(
                experiment.imageset,
                experiment.crystal,
                experiment.beam,
                experiment.detector,
                experiment.goniometer,
                experiment.scan,
                dmin=dmin,
                dmax=dmax,
                margin=margin,
                force_static=force_static,
                padding=padding,
            )
        from dials.algorithms.spot_prediction.reflection_predictor import (
            ReflectionPredictor,
        )

        predict = ReflectionPredictor(
            experiment,
            dmin=dmin,
            dmax=dmax,
            margin=margin,
            force_static=force_static,
            padding=padding,
        )
        return predict()

    @staticmethod
    def from_predictions_multi(
        experiments, dmin=None, dmax=None, margin=1, force_static=False, padding=0
    ):
        """
        Construct a reflection table from predictions.

        :param experiments: The experiment list to predict from
        :param dmin: The maximum resolution
        :param dmax: The minimum resolution
        :param margin: The margin to predict around
        :param force_static: Do static prediction with a scan varying model
        :param padding: Padding in degrees
        :return: The reflection table of predictions
        """
        result = dials_array_family_flex_ext.reflection_table()
        for i, e in enumerate(experiments):
            rlist = dials_array_family_flex_ext.reflection_table.from_predictions(
                e,
                dmin=dmin,
                dmax=dmax,
                margin=margin,
                force_static=force_static,
                padding=padding,
            )
            rlist["id"] = cctbx.array_family.flex.int(len(rlist), i)
            if e.identifier:
                rlist.experiment_identifiers()[i] = e.identifier
            result.extend(rlist)
        return result

    @staticmethod
    def from_observations(experiments, params=None, is_stills=False):
        """
        Construct a reflection table from observations.

        :param experiments: The experiments
        :param params: The input parameters
        :param is_stills:   [ADVANCED] Force still-handling of experiment
                            ID remapping for dials.stills_process. Do
                            not use for general processing unless you
                            know all the implications.
        :return: The reflection table of observations
        """
        from dials.algorithms.spot_finding.factory import SpotFinderFactory

        if params is None:
            from dials.command_line.find_spots import phil_scope
            from dials.util.phil import parse

            params = phil_scope.fetch(source=parse("")).extract()

        if params.spotfinder.filter.min_spot_size is libtbx.Auto:
            detector = experiments[0].imageset.get_detector()
            if detector[0].get_type() == "SENSOR_PAD":
                # smaller default value for pixel array detectors
                params.spotfinder.filter.min_spot_size = 3
            else:
                params.spotfinder.filter.min_spot_size = 6
            logger.info(
                "Setting spotfinder.filter.min_spot_size=%i",
                params.spotfinder.filter.min_spot_size,
            )

        # Set images to exclude in the imagesets
        if params.spotfinder.exclude_images_multiple:
            params.spotfinder.exclude_images = expand_exclude_multiples(
                experiments,
                params.spotfinder.exclude_images_multiple,
                params.spotfinder.exclude_images,
            )
        experiments = set_invalid_images(experiments, params.spotfinder.exclude_images)

        # Get the spot-finder from the input parameters
        logger.info("Configuring spot finder from input parameters")
        spotfinder = SpotFinderFactory.from_parameters(
            experiments=experiments, params=params, is_stills=is_stills
        )

        # Find the spots
        return spotfinder.find_spots(experiments)

    @staticmethod
    def from_pickle(filename):
        """
        Read the reflection table from pickle file.

        :param filename: The pickle filename
        :return: The reflection table
        """
        if filename and hasattr(filename, "__fspath__"):
            filename = filename.__fspath__()
        with libtbx.smart_open.for_reading(filename, "rb") as infile:
            result = pickle.load(infile, encoding="bytes")
            assert isinstance(result, dials_array_family_flex_ext.reflection_table)
            return result

    def as_msgpack_file(self, filename):
        """
        Write the reflection table to file in msgpack format
        """
        if filename and hasattr(filename, "__fspath__"):
            filename = filename.__fspath__()
        with libtbx.smart_open.for_writing(filename, "wb") as outfile:
            self.as_msgpack_to_file(dials.util.ext.streambuf(python_file_obj=outfile))

    @staticmethod
    def from_msgpack_file(filename):
        """
        Read the reflection table from file in msgpack format
        """
        if filename and hasattr(filename, "__fspath__"):
            filename = filename.__fspath__()
        with libtbx.smart_open.for_reading(filename, "rb") as infile:
            return dials_array_family_flex_ext.reflection_table.from_msgpack(
                infile.read()
            )

    def as_file(self, filename):
        """
        Write the reflection table to file in either msgpack or pickle format
        """
        if os.getenv("DIALS_USE_PICKLE"):
            self.as_pickle(filename)
        else:
            self.as_msgpack_file(filename)

    @staticmethod
    def from_file(filename):
        """
        Read the reflection table from either pickle or msgpack
        """
        try:
            return dials_array_family_flex_ext.reflection_table.from_msgpack_file(
                filename
            )
        except RuntimeError:
            return dials_array_family_flex_ext.reflection_table.from_pickle(filename)

    @staticmethod
    def empty_standard(nrows):
        """
        Create an empty table of specified number of rows with most of the standard
        keys

        :param nrows: The number of rows to create
        :return: The reflection table
        """

        assert nrows > 0
        table = dials_array_family_flex_ext.reflection_table(nrows)

        # General properties
        table["flags"] = cctbx.array_family.flex.size_t(nrows, 0)
        table["id"] = cctbx.array_family.flex.int(nrows, 0)
        table["panel"] = cctbx.array_family.flex.size_t(nrows, 0)

        # Predicted properties
        table["miller_index"] = cctbx.array_family.flex.miller_index(nrows)
        table["entering"] = cctbx.array_family.flex.bool(nrows)
        table["s1"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        table["xyzcal.mm"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        table["xyzcal.px"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        # table['ub_matrix'] = cctbx.array_family.flex.mat3_double(nrows, (0, 0, 0, 0, 0, 0, 0, 0, 0))

        # Observed properties
        table["xyzobs.px.value"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        table["xyzobs.px.variance"] = cctbx.array_family.flex.vec3_double(
            nrows, (0, 0, 0)
        )
        table["xyzobs.mm.value"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        table["xyzobs.mm.variance"] = cctbx.array_family.flex.vec3_double(
            nrows, (0, 0, 0)
        )
        table["rlp"] = cctbx.array_family.flex.vec3_double(nrows, (0, 0, 0))
        table["intensity.sum.value"] = cctbx.array_family.flex.double(nrows, 0)
        table["intensity.sum.variance"] = cctbx.array_family.flex.double(nrows, 0)
        table["intensity.prf.value"] = cctbx.array_family.flex.double(nrows, 0)
        table["intensity.prf.variance"] = cctbx.array_family.flex.double(nrows, 0)
        table["lp"] = cctbx.array_family.flex.double(nrows, 0)
        table["profile.correlation"] = cctbx.array_family.flex.double(nrows, 0)

        return table

    @staticmethod
    def plot(table, detector, key):
        """
        Plot a reflection table using matplotlib

        :param table: The reflection table
        :param detector: The detector model
        :param key: The key to plot
        """
        from matplotlib import pyplot as plt
        from matplotlib.patches import Polygon

        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")
        spots = table[key]
        if "px" in key:
            spots = [
                detector[table["panel"][i]].get_pixel_lab_coord(spots[i][0:2])
                for i in range(len(spots))
            ]
        else:
            assert "mm" in key
            spots = [
                detector[table["panel"][i]].get_lab_coord(spots[i][0:2])
                for i in range(len(spots))
            ]

        min_f = max_f = min_s = max_s = 0

        for i, panel in enumerate(detector):
            fs, ss = panel.get_image_size()
            p0 = panel.get_pixel_lab_coord((0, 0))
            p1 = panel.get_pixel_lab_coord((fs - 1, 0))
            p2 = panel.get_pixel_lab_coord((fs - 1, ss - 1))
            p3 = panel.get_pixel_lab_coord((0, ss - 1))
            p = Polygon(
                (p0[0:2], p1[0:2], p2[0:2], p3[0:2]),
                closed=True,
                color="green",
                fill=False,
                hatch="/",
            )

            if p.xy[:, 0].min() < min_f:
                min_f = p.xy[:, 0].min()
            if p.xy[:, 0].max() > max_f:
                max_f = p.xy[:, 0].max()
            if p.xy[:, 1].min() < min_s:
                min_s = p.xy[:, 1].min()
            if p.xy[:, 1].max() > max_s:
                max_s = p.xy[:, 1].max()

            ax.add_patch(p)

        ax.set_xlim((min_f - 10, max_f + 10))
        ax.set_ylim((min_s - 10, max_s + 10))
        plt.scatter([s[0] for s in spots], [s[1] for s in spots], c="blue", linewidth=0)
        plt.show()

    def as_pickle(self, filename):
        """
        Write the reflection table as a pickle file.

        :param filename: The output filename
        """
        # Clean up any removed experiments from the identifiers map
        self.clean_experiment_identifiers_map()

        with libtbx.smart_open.for_writing(filename, "wb") as outfile:
            pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    def as_h5(self, filename):
        """
        Write the reflection table as a HDF5 file.

        :param filename: The output filename
        """
        from dials.util.nexus_old import NexusFile

        handle = NexusFile(filename, "w")
        # Clean up any removed experiments from the identifiers map
        self.clean_experiment_identifiers_map()
        handle.set_reflections(self)
        handle.close()

    def as_miller_array(self, experiment, intensity="sum"):
        """Return a miller array with the chosen intensities.

        Use the provided experiment object and intensity choice to make a miller
        intensity array with sigmas (no scaling applied).

        Args:
            experiment (dxtbx.model.Experiment): An experiment object.
            intensity (str): The intensity type that will be used to make the
                miller array e.g 'prf', 'sum'.

        Returns:
            cctbx.miller.array: A miller array with intensities and sigmas.

        Raises:
            KeyError: If chosen intensity values cannot be found in the table.
        """

        try:
            intensities, variances = (
                self["intensity." + intensity + ".value"],
                self["intensity." + intensity + ".variance"],
            )
        except KeyError as e:
            logger.error(e, exc_info=True)
            raise KeyError(
                "Unable to find %s, %s in reflection table"
                % (
                    "intensity." + intensity + ".value",
                    "intensity." + intensity + ".variance",
                )
            )

        miller_set = cctbx.miller.set(
            crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
            indices=self["miller_index"],
            anomalous_flag=False,
        )
        i_obs = cctbx.miller.array(miller_set, data=intensities)
        i_obs.set_observation_type_xray_intensity()
        i_obs.set_sigmas(cctbx.array_family.flex.sqrt(variances))
        i_obs.set_info(
            cctbx.miller.array_info(source="DIALS", source_type="reflection_tables")
        )
        return i_obs

    def copy(self):
        """
        Copy everything.

        :return: A copy of the reflection table
        """
        return self.select(cctbx.array_family.flex.bool(len(self), True))

    def sort(self, name, reverse=False, order=None):
        """
        Sort the reflection table by a key.

        :param name: The name of the column
        :param reverse: Reverse the sort order
        :param order: For multi element items specify order
        """

        if type(self[name]) in (
            cctbx.array_family.flex.vec2_double,
            cctbx.array_family.flex.vec3_double,
            cctbx.array_family.flex.mat3_double,
            dials_array_family_flex_ext.int6,
            cctbx.array_family.flex.miller_index,
        ):
            data = self[name]
            if not order:
                perm = cctbx.array_family.flex.size_t(
                    sorted(range(len(self)), key=lambda x: data[x], reverse=reverse)
                )
            else:
                assert len(order) == len(data[0])
                perm = cctbx.array_family.flex.size_t(
                    sorted(
                        range(len(self)),
                        key=lambda x: tuple(data[x][i] for i in order),
                        reverse=reverse,
                    )
                )
        else:
            perm = cctbx.array_family.flex.sort_permutation(
                self[name], reverse=reverse, stable=True
            )
        self.reorder(perm)

    """
    Sorting the reflection table within an already sorted column
    """

    def subsort(self, key0, key1, reverse=False):
        """
        Sort the reflection based on key1 within a constant key0.

        :param key0: The name of the column values to sort within
        :param key1: The sorting key name within the selected column
        """
        uniq_values = self[key0]
        for ii in set(uniq_values):
            val = (uniq_values == ii).iselection()
            ref_tmp = copy.deepcopy(self[min(val) : (max(val) + 1)])
            ref_tmp.sort(key1, reverse)
            self[min(val) : (max(val) + 1)] = ref_tmp

    def match_by_hkle(
        self, other: dials_array_family_flex_ext.reflection_table
    ) -> Tuple[cctbx.array_family.flex.size_t, cctbx.array_family.flex.size_t]:
        """
        Match reflections with another set of reflections by the h, k, l
        and entering values. Uses pandas dataframe merge method to match
        the columns: assumes the key h, k, l, e is unique which is false
        if > 360 degree rotation.

        Args:
            other: reflection table to match against

        Returns:
            Indices in self, indices in other for matches
        """

        hkl = self["miller_index"].as_vec3_double().parts()
        hkl = (part.as_numpy_array().astype(int) for part in hkl)
        e = self["entering"].as_numpy_array().astype(int)
        n = np.arange(e.size)
        p0 = pd.DataFrame(dict(zip("hklen", (*hkl, e, n))), copy=False)

        hkl = other["miller_index"].as_vec3_double().parts()
        hkl = (part.as_numpy_array().astype(int) for part in hkl)
        e = other["entering"].as_numpy_array().astype(int)
        n = np.arange(e.size)
        p1 = pd.DataFrame(dict(zip("hklen", (*hkl, e, n))), copy=False)

        merged = pd.merge(p0, p1, on=["h", "k", "l", "e"], suffixes=[0, 1])

        n0 = cctbx.array_family.flex.size_t(merged.n0.values)
        n1 = cctbx.array_family.flex.size_t(merged.n1.values)

        return n0, n1

    @staticmethod
    def concat(
        tables: List[dials_array_family_flex_ext.reflection_table],
    ) -> dials_array_family_flex_ext.reflection_table:
        """
        Concatenate a list of reflection tables, taking care to correctly handle
        experiment identifiers and ids.
        :param tables: A list of reflection tables
        :return: A single combined reflection table
        """
        from dials.util.multi_dataset_handling import renumber_table_id_columns

        tables = renumber_table_id_columns(tables)
        first = tables[0]
        for table in tables[1:]:
            first.extend(table)
        return first

    def match_with_reference(self, other):
        """
        Match reflections with another set of reflections.

        :param other: The reflection table to match against
        :return: The matches
        """
        logger.info("Matching reference spots with predicted reflections")
        logger.info(" %d observed reflections input", len(other))
        logger.info(" %d reflections predicted", len(self))

        # Get the miller index, entering flag and turn number for
        # Both sets of reflections
        i1 = self["id"]
        h1 = self["miller_index"]
        e1 = self["entering"].as_int()
        x1, y1, z1 = self["xyzcal.px"].parts()
        p1 = self["panel"]

        i2 = other["id"]
        h2 = other["miller_index"]
        e2 = other["entering"].as_int()
        x2, y2, z2 = other["xyzcal.px"].parts()
        p2 = other["panel"]

        class Match:
            def __init__(self):
                self.a = []
                self.b = []

        # Create the match lookup
        lookup = collections.defaultdict(Match)
        for i in range(len(self)):
            item = h1[i] + (e1[i], i1[i], p1[i])
            lookup[item].a.append(i)

        # Add matches from input reflections
        for i in range(len(other)):
            item = h2[i] + (e2[i], i2[i], p2[i])
            if item in lookup:
                lookup[item].b.append(i)

        # Create the list of matches
        match1 = []
        match2 = []
        for item, value in lookup.items():
            if len(value.b) == 0:
                continue
            elif len(value.a) == 1 and len(value.b) == 1:
                match1.append(value.a[0])
                match2.append(value.b[0])
            else:
                matched = {}
                for i in value.a:
                    d = []
                    for j in value.b:
                        dx = x1[i] - x2[j]
                        dy = y1[i] - y2[j]
                        dz = z1[i] - z2[j]
                        d.append((i, j, dx**2 + dy**2 + dz**2))
                    i, j, d = min(d, key=lambda x: x[2])
                    if j not in matched:
                        matched[j] = (i, d)
                    elif d < matched[j][1]:
                        matched[j] = (i, d)
                for key1, value1 in matched.items():
                    match1.append(value1[0])
                    match2.append(key1)

        # Select everything which matches
        sind = cctbx.array_family.flex.size_t(match1)
        oind = cctbx.array_family.flex.size_t(match2)

        # Sort by self index
        sort_index = cctbx.array_family.flex.size_t(
            sorted(range(len(sind)), key=lambda x: sind[x])
        )
        sind = sind.select(sort_index)
        oind = oind.select(sort_index)

        s2 = self.select(sind)
        o2 = other.select(oind)
        h1 = s2["miller_index"]
        h2 = o2["miller_index"]
        e1 = s2["entering"]
        e2 = o2["entering"]
        assert (h1 == h2).all_eq(True)
        assert (e1 == e2).all_eq(True)
        x1, y1, z1 = s2["xyzcal.px"].parts()
        x2, y2, z2 = o2["xyzcal.px"].parts()
        distance = cctbx.array_family.flex.sqrt(
            cctbx.array_family.flex.pow2(x1 - x2)
            + cctbx.array_family.flex.pow2(y1 - y2)
            + cctbx.array_family.flex.pow2(z1 - z2)
        )
        mask = distance < 2
        logger.info(" %d reflections matched", len(o2))
        logger.info(" %d reflections accepted", mask.count(True))
        self.set_flags(sind.select(mask), self.flags.reference_spot)
        self.set_flags(sind.select(o2.get_flags(self.flags.strong)), self.flags.strong)
        self.set_flags(
            sind.select(o2.get_flags(self.flags.indexed)), self.flags.indexed
        )
        self.set_flags(
            sind.select(o2.get_flags(self.flags.used_in_refinement)),
            self.flags.used_in_refinement,
        )
        other_matched_indices = oind.select(mask)
        other_unmatched_mask = cctbx.array_family.flex.bool(len(other), True)
        other_unmatched_mask.set_selected(
            other_matched_indices,
            cctbx.array_family.flex.bool(len(other_matched_indices), False),
        )
        other_matched = other.select(other_matched_indices)
        other_unmatched = other.select(other_unmatched_mask)
        for key, column in self.select(sind.select(mask)).cols():
            other_matched[key] = column
        mask2 = cctbx.array_family.flex.bool(len(self), False)
        mask2.set_selected(sind.select(mask), True)
        return mask2, other_matched, other_unmatched

    def match(
        self,
        other: dials_array_family_flex_ext.reflection_table,
        *,
        max_separation: int = 2,
        key: str = "xyzobs.px.value",
        scale: Tuple[float, float, float] = (1, 1, 1),
    ) -> Tuple[
        cctbx.array_family.flex.int,
        cctbx.array_family.flex.int,
        cctbx.array_family.flex.double,
    ]:
        """
        Match reflections from this reflection list to another (reference) list

        The match is based on comparison of any chosen 3-vector column.

        Args:
            other: The reflection table to match against
            max_separations: Maximum distance in column-space to match on
            key: Name of 3-vector column
            scale:
                Relative scales to apply to vectors before matching in
                case e.g. very wide or very fine slicing.

        Returns:
            (self_index, other_index, distance) where:
                self_index: is the reindex in the this reflection table
                other_index: is the index in the other reflection table
                distance: The calculated distance of specified column-vector

            such that self[self_index][j] ~= other[other_index][j].
        """

        xyz = self[key]
        ref = other[key]

        if scale != (1, 1, 1):
            x, y, z = xyz.parts()
            x *= scale[0]
            y *= scale[1]
            z *= scale[2]
            xyz = cctbx.array_family.flex.vec3_double(x, y, z)

            x, y, z = ref.parts()
            x *= scale[0]
            y *= scale[1]
            z *= scale[2]
            ref = cctbx.array_family.flex.vec3_double(x, y, z)

        x = xyz.as_double().as_1d()
        r = ref.as_double().as_1d()
        ann = AnnAdaptorSelfInclude(r, 3)
        ann.query(x)

        mm = cctbx.array_family.flex.size_t_range(xyz.size())
        nn, distance = ann.nn.as_size_t(), cctbx.array_family.flex.sqrt(ann.distances)

        sel = distance <= max_separation

        mm = mm.select(sel)
        nn = nn.select(sel)
        distance = distance.select(sel)
        return mm, nn, distance

    def compute_zeta(self, experiment):
        """
        Compute zeta for each reflection.

        :param experiment: The experimental models
        :return: Zeta for each reflection
        """
        from dials.algorithms.profile_model.gaussian_rs import zeta_factor

        m2 = experiment.goniometer.get_rotation_axis()
        s0 = experiment.beam.get_s0()
        self["zeta"] = zeta_factor(m2, s0, self["s1"])
        return self["zeta"]

    def compute_zeta_multi(self, experiments):
        """
        Compute zeta for each reflection.

        :param experiments: The list of experiments
        :return: Zeta for each reflection
        """
        from dials.algorithms.profile_model.gaussian_rs import zeta_factor

        m2 = cctbx.array_family.flex.vec3_double(len(experiments))
        s0 = cctbx.array_family.flex.vec3_double(len(experiments))
        for i, e in enumerate(experiments):
            m2[i] = e.goniometer.get_rotation_axis()
            s0[i] = e.beam.get_s0()
        self["zeta"] = zeta_factor(m2, s0, self["s1"], self["id"])
        return self["zeta"]

    def compute_d_single(self, experiment):
        """
        Compute the resolution for each reflection.

        :param experiment: The experimental models
        :return: The resolution for each reflection
        """
        uc = dials_array_family_flex_ext.unit_cell(1)
        uc[0] = experiment.crystal.get_unit_cell()
        self["d"] = uc.d(
            self["miller_index"], cctbx.array_family.flex.size_t(len(self), 0)
        )
        return self["d"]

    def compute_d(self, experiments):
        """
        Compute the resolution for each reflection.

        :param experiments: The experiment list
        :return: The resolution for each reflection
        """
        uc = dials_array_family_flex_ext.unit_cell(len(experiments))
        for i, e in enumerate(experiments):
            uc[i] = e.crystal.get_unit_cell()
        assert self["id"].all_ge(0)
        self["d"] = uc.d(
            self["miller_index"], cctbx.array_family.flex.size_t(list(self["id"]))
        )
        return self["d"]

    def compute_bbox(self, experiments, sigma_b_multiplier=2.0):
        """
        Compute the bounding boxes.

        :param experiments: The list of experiments
        :param profile_model: The profile models
        :param sigma_b_multiplier: Multiplier to cover extra background
        :return: The bounding box for each reflection
        """
        self["bbox"] = dials_array_family_flex_ext.int6(len(self))
        for expr, indices in self.iterate_experiments_and_indices(experiments):
            self["bbox"].set_selected(
                indices,
                expr.profile.compute_bbox(
                    self.select(indices),
                    expr.crystal,
                    expr.beam,
                    expr.detector,
                    expr.goniometer,
                    expr.scan,
                    sigma_b_multiplier=sigma_b_multiplier,
                ),
            )
        return self["bbox"]

    def compute_partiality(self, experiments):
        """
        Compute the reflection partiality.

        :param experiments: The experiment list
        :param profile_model: The profile models
        :return: The partiality for each reflection
        """
        self["partiality"] = cctbx.array_family.flex.double(len(self))
        for expr, indices in self.iterate_experiments_and_indices(experiments):
            self["partiality"].set_selected(
                indices,
                expr.profile.compute_partiality(
                    self.select(indices),
                    expr.crystal,
                    expr.beam,
                    expr.detector,
                    expr.goniometer,
                    expr.scan,
                ),
            )
        return self["partiality"]

    def compute_mask(self, experiments, image_volume=None, overlaps=None):
        """
        Apply a mask to the shoeboxes.

        :param experiments: The list of experiments
        :param profile_model: The profile model
        """
        for expr, indices in self.iterate_experiments_and_indices(experiments):
            result = expr.profile.compute_mask(
                self.select(indices),
                expr.crystal,
                expr.beam,
                expr.detector,
                expr.goniometer,
                expr.scan,
                image_volume=image_volume,
            )
            if result is not None:
                if "fraction" not in self:
                    self["fraction"] = cctbx.array_family.flex.double(len(self))
                self["fraction"].set_selected(indices, result)

    def iterate_experiments_and_indices(self, experiments):
        """
        A helper function to iterate through experiments and indices of reflections
        for each experiment
        """
        assert len(experiments) > 0
        index_list = self.split_indices_by_experiment_id(len(experiments))
        assert len(experiments) == len(index_list)
        tot = 0
        for l in index_list:
            tot += len(l)
        assert tot == len(self)
        yield from zip(experiments, index_list)

    def random_split(self, n=2):
        """Randomly split table into n tables.

        Not all tables will be the same length."""
        n = int(n)
        if n <= 1:
            return [self]
        list_of_reflection_tables = []

        perm = cctbx.array_family.flex.random_permutation(self.size())
        divisions = [int(i * self.size() / n) for i in range(0, n + 1)]
        for i, v in enumerate(divisions[:-1]):
            isel = perm[v : divisions[i + 1]]
            if isel:
                list_of_reflection_tables.append(self.select(isel))
        return list_of_reflection_tables

    def compute_background(self, experiments, image_volume=None):
        """
        Helper function to compute the background.

        :param experiments: The list of experiments
        """
        success = self.background_algorithm(experiments).compute_background(
            self, image_volume
        )
        self.set_flags(~success, self.flags.failed_during_background_modelling)

    def compute_centroid(self, experiments, image_volume=None):
        """
        Helper function to compute the centroid.

        :param experiments: The list of experiments
        """
        self.centroid_algorithm(experiments).compute_centroid(
            self, image_volume=image_volume
        )

    def compute_summed_intensity(
        self, image_volume: dials.model.data.MultiPanelImageVolume = None
    ) -> None:
        """
        Compute intensity via summation integration.
        """
        from dials.algorithms.integration.sum import sum_integrate_and_update_table

        success = sum_integrate_and_update_table(self, image_volume=image_volume)
        self.set_flags(~success, self.flags.failed_during_summation)

    def compute_fitted_intensity(self, fitter):
        """
        Helper function to compute the intensity.

        :param experiments: The list of experiments
        :param profile_model: The profile model
        """
        success = fitter.fit(self)
        self.set_flags(~success, self.flags.failed_during_profile_fitting)

    def compute_corrections(self, experiments):
        """
        Helper function to correct the intensity.

        :param experiments: The list of experiments
        :return: The LP correction for each reflection
        """
        from dials.algorithms.integration import Corrections, CorrectionsMulti

        compute = CorrectionsMulti()
        for experiment in experiments:
            if (
                experiment.goniometer is not None
                and experiment.scan is not None
                and (experiment.scan.get_oscillation()[1] != 0.0)
            ):
                compute.append(
                    Corrections(
                        experiment.beam, experiment.goniometer, experiment.detector
                    )
                )
            else:
                compute.append(Corrections(experiment.beam, experiment.detector))
        lp = compute.lp(self["id"], self["s1"])
        self["lp"] = lp
        if experiment.detector[0].get_mu() > 0:
            qe = compute.qe(self["id"], self["s1"], self["panel"])
            self["qe"] = qe
        return lp

    def extract_shoeboxes(self, imageset, mask=None, nthreads=1):
        """
        Helper function to read a load of shoebox data.

        :param imageset: The imageset
        :param mask: The mask to apply
        :param nthreads: The number of threads to use
        :return: A tuple containing read time and extract time
        """
        from time import time

        from dials.model.data import make_image

        assert "shoebox" in self
        detector = imageset.get_detector()
        try:
            frame0, frame1 = imageset.get_array_range()
        except Exception:
            frame0, frame1 = (0, len(imageset))
        extractor = dials_array_family_flex_ext.ShoeboxExtractor(
            self, len(detector), frame0, frame1
        )
        logger.info(" Beginning to read images")
        read_time = 0
        extract_time = 0
        for i in range(len(imageset)):
            logger.debug("  reading image %d", i)
            st = time()
            image = imageset.get_corrected_data(i)
            mask2 = imageset.get_mask(i)
            if mask is not None:
                assert len(mask) == len(mask2)
                mask2 = tuple(m1 & m2 for m1, m2 in zip(mask, mask2))
            read_time += time() - st
            st = time()
            extractor.next(make_image(image, mask2))
            extract_time += time() - st
            del image
        assert extractor.finished()
        logger.info("  successfully read %d images", frame1 - frame0)
        logger.info("  read time: %.1f seconds", read_time)
        logger.info("  extract time: %.1f seconds", extract_time)
        return read_time, extract_time

    def is_overloaded(self, experiments):
        """
        Check if the shoebox contains overloaded pixels.

        :param experiments: The experiment list
        :return: True/False overloaded for each reflection
        """
        from dials.algorithms.shoebox import OverloadChecker

        assert "shoebox" in self
        assert "id" in self
        detectors = [expr.detector for expr in experiments]
        checker = OverloadChecker()
        for detector in detectors:
            checker.add(
                cctbx.array_family.flex.double(
                    p.get_trusted_range()[1] for p in detector
                )
            )
        result = checker(self["id"], self["shoebox"])
        self.set_flags(result, self.flags.overloaded)
        return result

    def contains_invalid_pixels(self):
        """
        Check if the shoebox contains invalid pixels.

        :return: True/False invalid for each reflection
        """
        from dials.algorithms.shoebox import MaskCode

        assert "shoebox" in self
        x0, x1, y0, y1, z0, z1 = self["bbox"].parts()
        ntotal = (x1 - x0) * (y1 - y0) * (z1 - z0)
        assert ntotal.all_gt(0)
        sbox = self["shoebox"]
        nvalid = sbox.count_mask_values(MaskCode.Valid)
        nbackg = sbox.count_mask_values(MaskCode.Background)
        nforeg = sbox.count_mask_values(MaskCode.Foreground)
        nvalbg = sbox.count_mask_values(MaskCode.Background | MaskCode.Valid)
        nvalfg = sbox.count_mask_values(MaskCode.Foreground | MaskCode.Valid)
        ninvbg = nbackg - nvalbg
        ninvfg = nforeg - nvalfg
        assert ninvbg.all_ge(0)
        assert ninvfg.all_ge(0)
        self.set_flags(ninvbg > 0, self.flags.background_includes_bad_pixels)
        self.set_flags(ninvfg > 0, self.flags.foreground_includes_bad_pixels)
        return (ntotal - nvalid) > 0

    def find_overlaps(self, experiments=None, border=0):
        """
        Check for overlapping reflections.

        :param experiments: The experiment list
        :param tolerance: A positive integer specifying border around shoebox
        :return: The overlap list
        """
        from dials.algorithms.shoebox import OverlapFinder

        # Expand the bbox if necessary
        if border > 0:
            x0, x1, y0, y1, z0, z1 = self["bbox"].parts()
            x0 -= border
            x1 += border
            y0 -= border
            y1 += border
            z0 -= border
            z1 += border
            bbox = dials_array_family_flex_ext.int6(x0, x1, y0, y1, z0, z1)
        else:
            bbox = self["bbox"]

        # Get the panel and id
        panel = self["panel"]

        # Group according to imageset
        if experiments is not None:
            groups = itertools.groupby(
                range(len(experiments)), lambda x: experiments[x].imageset
            )

            # Get the experiment ids we're to treat together
            lookup = {}
            for j, (key, indices) in enumerate(groups):
                for i in indices:
                    lookup[i] = j
            group_id = cctbx.array_family.flex.size_t([lookup[i] for i in self["id"]])
        elif "imageset_id" in self:
            imageset_id = self["imageset_id"]
            assert imageset_id.all_ge(0)
            group_id = cctbx.array_family.flex.size_t(list(imageset_id))
        else:
            raise RuntimeError("Either need to supply experiments or have imageset_id")

        # Create the overlap finder
        find_overlapping = OverlapFinder()

        # Find the overlaps
        overlaps = find_overlapping(group_id, panel, bbox)
        assert overlaps.num_vertices() == len(self)

        # Return the overlaps
        return overlaps

    def assert_experiment_identifiers_are_consistent(self, experiments=None):
        """
        Check the experiment identifiers
        """
        identifiers = self.experiment_identifiers()
        if len(identifiers) > 0:
            values = list(identifiers.values())
            assert len(set(values)) == len(values), (len(set(values)), len(values))
            if "id" in self:
                index = set(self["id"]).difference({-1})
                for i in index:
                    assert i in identifiers, (i, list(identifiers))
        if experiments is not None:
            if len(identifiers) > 0:
                assert len(identifiers) == len(experiments), (
                    len(identifiers),
                    len(experiments),
                )
                assert len(identifiers) == len(set(experiments.identifiers()))
                for experiment in experiments:
                    assert (
                        experiment.identifier in identifiers.values()
                    ), experiment.identifier

    def are_experiment_identifiers_consistent(self, experiments=None):
        """
        Check the experiment identifiers
        """
        try:
            self.assert_experiment_identifiers_are_consistent(experiments)
        except AssertionError:
            return False
        return True

    def compute_miller_indices_in_asu(self, experiments):
        """
        Compute miller indices in the asu
        """
        self["miller_index_asu"] = cctbx.array_family.flex.miller_index(len(self))
        for idx, experiment in enumerate(experiments):

            # Create the crystal symmetry object
            uc = experiment.crystal.get_unit_cell()
            sg = experiment.crystal.get_space_group()
            cs = cctbx.crystal.symmetry(uc, space_group=sg)

            # Get the selection and compute the miller indices
            selection = self["id"] == idx
            h = self["miller_index"].select(selection)
            ms = cctbx.miller.set(cs, h)
            ms_asu = ms.map_to_asu()
            h_asu = ms_asu.indices()

            # Set the miller indices
            self["miller_index_asu"].set_selected(selection, h_asu)

        return self["miller_index_asu"]

    def select_on_experiment_identifiers(self, list_of_identifiers):
        """
        Given a list of experiment identifiers (strings), perform a selection
        and return a reflection table with properly configured experiment_identifiers
        map.
        """
        # First get the reverse of the map i.e. ids for a given exp_identifier
        id_values = []
        for k, v in zip(
            self.experiment_identifiers().keys(), self.experiment_identifiers().values()
        ):
            if v in list_of_identifiers:
                id_values.append(k)
        if len(id_values) != len(list_of_identifiers):
            logger.warning(
                """Not all requested identifiers
found in the table's map, has the experiment_identifiers() map been created?
Requested %s:
Found %s"""
                % (list_of_identifiers, id_values)
            )
        # Build up a selection and use this
        sel = cctbx.array_family.flex.bool(self.size(), False)
        for id_val in id_values:
            id_sel = self["id"] == id_val
            sel.set_selected(id_sel, True)
        self = self.select(sel)
        # Remove entries from the experiment_identifiers map
        for k in self.experiment_identifiers().keys():
            if k not in id_values:
                del self.experiment_identifiers()[k]
        return self

    def remove_on_experiment_identifiers(self, list_of_identifiers):
        """
        Remove datasets from the table, given a list of experiment
        identifiers (strings).
        """
        # First get the reverse of the map i.e. ids for a given exp_identifier
        if len(self):
            assert "id" in self
        id_values = []
        for k, v in zip(
            self.experiment_identifiers().keys(), self.experiment_identifiers().values()
        ):
            if v in list_of_identifiers:
                id_values.append(k)
        if len(id_values) != len(list_of_identifiers):
            logger.warning(
                """Not all requested identifiers
found in the table's map, has the experiment_identifiers() map been created?
Requested %s:
Found %s"""
                % (list_of_identifiers, id_values)
            )
        # Now delete the selections, also removing the entry from the map
        for id_val in id_values:
            sel = self["id"] == id_val
            self.del_selected(sel)
            del self.experiment_identifiers()[id_val]
        return self

    def clean_experiment_identifiers_map(self):
        """
        Remove any entries from the identifier map that do not have any
        data in the table. Primarily to call as saving data to give a
        consistent table and map.
        """
        dataset_ids_in_table = set(self["id"]).difference({-1})
        dataset_ids_in_map = set(self.experiment_identifiers().keys())
        ids_to_remove = dataset_ids_in_map.difference(dataset_ids_in_table)
        for i in ids_to_remove:
            del self.experiment_identifiers()[i]

    def reset_ids(self):
        """
        Reset the 'id' column such that the experiment identifiers are
        numbered 0 .. n-1.
        """
        reverse_map = {v: k for k, v in self.experiment_identifiers()}
        for k in self.experiment_identifiers().keys():
            del self.experiment_identifiers()[k]
        if not len(self):
            return
        orig_id = self["id"].deep_copy()
        for i_exp, exp_id in enumerate(reverse_map.keys()):
            sel_exp = orig_id == reverse_map[exp_id]
            self["id"].set_selected(sel_exp, i_exp)
            self.experiment_identifiers()[i_exp] = exp_id

    def centroid_px_to_mm(self, experiments):
        """
        Map spot centroids from pixel/image number to mm/radian.

        Used to convert spot centroids coming from e.g. dials.find_spots which are
        in pixel/image number units to mm/radian units as required for indexing and
        refinement.

        Args:
          experiments (dxtbx.model.ExperimentList): A list of experiments.
        """

        self["xyzobs.mm.value"] = cctbx.array_family.flex.vec3_double(len(self))
        self["xyzobs.mm.variance"] = cctbx.array_family.flex.vec3_double(len(self))
        # e.g. data imported from XDS; no variance known then; since is used
        # only for weights assign as 1 => uniform weights
        if "xyzobs.px.variance" not in self:
            self["xyzobs.px.variance"] = cctbx.array_family.flex.vec3_double(
                len(self), (1, 1, 1)
            )
        panel_numbers = cctbx.array_family.flex.size_t(self["panel"])

        for i, expt in enumerate(experiments):
            if "imageset_id" in self:
                sel_expt = self["imageset_id"] == i
            else:
                sel_expt = self["id"] == i
            for i_panel in range(len(expt.detector)):
                sel = sel_expt & (panel_numbers == i_panel)
                centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
                    expt.detector[i_panel],
                    expt.scan,
                    self["xyzobs.px.value"].select(sel),
                    self["xyzobs.px.variance"].select(sel),
                    cctbx.array_family.flex.vec3_double(sel.count(True), (1, 1, 1)),
                )
                self["xyzobs.mm.value"].set_selected(sel, centroid_position)
                self["xyzobs.mm.variance"].set_selected(sel, centroid_variance)

    def map_centroids_to_reciprocal_space(
        self, experiments, calculated=False, crystal_frame=False
    ):
        """Map mm/radian spot centroids to reciprocal space.

        Used to convert spot centroids provided in mm/radian units to reciprocal space
        as required for indexing. Adds the column 'rlp' to the reflection table, which
        contains a :py:class:`.flex.vec3_double` array of the reciprocal lattice vectors.

        Args:
          experiments (dxtbx.model.ExperimentList): A list of experiments.
          calculated (Bool): use calculated positions.
          crystal_frame (Bool): return x, y, z positions as divided by U matrix
        """

        self["s1"] = cctbx.array_family.flex.vec3_double(len(self))
        self["rlp"] = cctbx.array_family.flex.vec3_double(len(self))
        panel_numbers = cctbx.array_family.flex.size_t(self["panel"])

        # crystal_frame is not used in indexing, but is a feature of the
        # reciprocal lattice viewer -> but if we are looking at the crystal
        # coordinate frame we need to look at the experiments independently

        for i, expt in enumerate(experiments):
            if not crystal_frame and "imageset_id" in self:
                sel_expt = self["imageset_id"] == i
            else:
                sel_expt = self["id"] == i

            for i_panel in range(len(expt.detector)):
                sel = sel_expt & (panel_numbers == i_panel)
                if calculated:
                    x, y, rot_angle = self["xyzcal.mm"].select(sel).parts()
                else:
                    x, y, rot_angle = self["xyzobs.mm.value"].select(sel).parts()
                s1 = expt.detector[i_panel].get_lab_coord(
                    cctbx.array_family.flex.vec2_double(x, y)
                )
                s1 = s1 / s1.norms() * (1 / expt.beam.get_wavelength())
                self["s1"].set_selected(sel, s1)
                S = s1 - expt.beam.get_s0()
                if expt.goniometer is not None:
                    setting_rotation = matrix.sqr(
                        expt.goniometer.get_setting_rotation()
                    )
                    rotation_axis = expt.goniometer.get_rotation_axis_datum()
                    sample_rotation = matrix.sqr(expt.goniometer.get_fixed_rotation())
                    if expt.crystal and crystal_frame:
                        sample_rotation *= matrix.sqr(expt.crystal.get_U())

                    self["rlp"].set_selected(sel, tuple(setting_rotation.inverse()) * S)
                    self["rlp"].set_selected(
                        sel,
                        self["rlp"]
                        .select(sel)
                        .rotate_around_origin(rotation_axis, -rot_angle),
                    )
                    self["rlp"].set_selected(
                        sel, tuple(sample_rotation.inverse()) * self["rlp"].select(sel)
                    )
                else:
                    self["rlp"].set_selected(sel, S)

    def calculate_entering_flags(self, experiments):
        """Calculate the entering flags for the reflections.

        Calculate a unit vector normal to the spindle-beam plane for this experiment,
        such that the vector placed at the centre of the Ewald sphere points to the
        hemispere in which reflections cross from inside to outside of the sphere
        (reflections are exiting). Adds the array of boolean entering flags to self
        as the "entering" column.

        Note:
            NB this vector is in +ve Y direction when using imgCIF coordinate frame.

        Args:

            experiments: The experiment list to use in calculating the entering flags.
        """

        assert "s1" in self

        # Init entering flags. These are always False for experiments that have no
        # rotation axis.
        enterings = cctbx.array_family.flex.bool(len(self), False)

        for iexp, experiment in enumerate(experiments):
            if not experiment.goniometer:
                continue
            axis = matrix.col(experiment.goniometer.get_rotation_axis())
            s0 = matrix.col(experiment.beam.get_s0())
            vec = s0.cross(axis)
            sel = self["id"] == iexp
            enterings.set_selected(sel, self["s1"].dot(vec) < 0.0)

        self["entering"] = enterings

    def get(self, key, default=None):
        """
        Get item from object for given key (ex: reflection_table column).

        Returns default value if not found.
        """
        if key in self:
            return self[key]
        return default


class reflection_table_selector:
    """
    A class to select columns from reflection table.

    This is mainly useful for specifying selections from phil parameters
    """

    def __init__(self, column, op, value):
        """
        Initialise the selector

        :param col: The column name
        :param op: The operator
        :param value: The value
        """
        # Set the column and value
        self.column = column
        self.value = value

        # Set the operator
        if isinstance(op, str):
            if op == "<":
                self.op = operator.lt
            elif op == "<=":
                self.op = operator.le
            elif op == "==":
                self.op = operator.eq
            elif op == "!=":
                self.op = operator.ne
            elif op == ">=":
                self.op = operator.ge
            elif op == ">":
                self.op = operator.gt
            elif op == "&":
                self.op = operator.and_
            else:
                raise RuntimeError("Unknown operator")
        else:
            self.op = op

    @property
    def op_string(self):
        """
        Return the operator as a string
        """

        if self.op == operator.lt:
            string = "<"
        elif self.op == operator.le:
            string = "<="
        elif self.op == operator.eq:
            string = "=="
        elif self.op == operator.ne:
            string = "!="
        elif self.op == operator.ge:
            string = ">="
        elif self.op == operator.gt:
            string = ">"
        elif self.op == operator.and_:
            string = "&"
        else:
            raise RuntimeError("Unknown operator")
        return string

    def __call__(self, reflections):
        """
        Select the reflections

        :param reflections: The reflections

        :return: The selection as a mask
        """
        if self.column == "intensity.sum.i_over_sigma":
            I = reflections["intensity.sum.value"]
            V = reflections["intensity.sum.variance"]
            mask1 = V > 0
            I = I.select(mask1)
            V = V.select(mask1)
            data = I / cctbx.array_family.flex.sqrt(V)
        elif self.column == "intensity.prf.i_over_sigma":
            I = reflections["intensity.prf.value"]
            V = reflections["intensity.prf.variance"]
            mask1 = V > 0
            I = I.select(mask1)
            V = V.select(mask1)
            data = I / cctbx.array_family.flex.sqrt(V)
        else:
            mask1 = None
            data = reflections[self.column]
        if isinstance(data, cctbx.array_family.flex.double):
            value = float(self.value)
        elif isinstance(data, cctbx.array_family.flex.int):
            value = int(self.value)
        elif isinstance(data, cctbx.array_family.flex.size_t):
            value = int(self.value)
        elif isinstance(data, cctbx.array_family.flex.std_string):
            value = self.value
        elif isinstance(data, cctbx.array_family.flex.vec3_double):
            raise RuntimeError("Comparison not implemented")
        elif isinstance(data, cctbx.array_family.flex.vec2_double):
            raise RuntimeError("Comparison not implemented")
        elif isinstance(data, cctbx.array_family.flex.mat3_double):
            raise RuntimeError("Comparison not implemented")
        elif isinstance(data, dials_array_family_flex_ext.int6):
            raise RuntimeError("Comparison not implemented")
        elif isinstance(data, dials_array_family_flex_ext.shoebox):
            raise RuntimeError("Comparison not implemented")
        else:
            raise RuntimeError("Unknown column type")
        mask2 = self.op(data, value)
        if mask1 is not None:
            mask1.set_selected(
                cctbx.array_family.flex.size_t(range(len(mask1))).select(mask1), mask2
            )
        else:
            mask1 = mask2
        return mask1
