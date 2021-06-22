#!/usr/bin/env python
#
# dials.potato.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

import json
import logging
from collections import OrderedDict
from math import floor, pi, sqrt

import numpy as np
from jinja2 import ChoiceLoader, Environment, PackageLoader

from libtbx.phil import parse
from scitbx import matrix

from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator, MaskCalculator
from dials.algorithms.profile_model.gaussian_rs.calculator import (
    ComputeEsdBeamDivergence,
)
from dials.algorithms.profile_model.potato import chisq_pdf, chisq_quantile
from dials.algorithms.profile_model.potato.model import ProfileModelFactory
from dials.algorithms.profile_model.potato.parameterisation import ModelState
from dials.algorithms.profile_model.potato.refiner import Refiner as ProfileRefiner
from dials.algorithms.profile_model.potato.refiner import RefinerData
from dials.algorithms.refinement.corrgram import corrgram
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.statistics.fast_mcd import FastMCD, maha_dist_sq
from dials.array_family import flex

logger = logging.getLogger("dials." + __name__)

# Parameters
phil_scope = parse(
    """

  profile
  {

    rlp_mosaicity {

      model = simple1 simple6 *angular2 angular4 angular6
        .type = choice

    }

    wavelength_spread {

      model = *delta gaussian
        .type = choice

    }

    unit_cell {

      fixed = False
        .type = bool

    }

    orientation {

      fixed = False
        .type = bool

    }

  }

  indexing {

    fail_on_bad_index = False
      .type = bool

  }

  refinement {

    max_separation = 2
      .type = float

    outlier_probability = 0.975
      .type = float

    n_macro_cycles = 3
      .type = int

    n_cycles = 3
      .type = int

    min_n_reflections=10
      .type = int

  }

  prediction {
    d_min = None
      .type = float

    probability = 0.9973
      .type = float
  }

  integration {

    use_crude_shoebox_mask = False
      .type = bool

    shoebox {

      probability = 0.9973
        .type = float

    }

    corrections {
      lp = True
        .type = bool
      dqe = True
        .type = bool
      partiality = True
        .type = bool
    }

  }

  debug {
    output {

      strong_spots = False
        .type = bool

      shoeboxes = True
        .type = bool

      profile_model = True
        .type = bool

      history = True
        .type = bool

      plots = False
        .type = bool

      print_shoeboxes = False
        .type = bool
    }
  }


"""
)


class Indexer(object):
    """
    A class to reindex the strong spot list

    """

    def __init__(self, params, experiments, reflections):
        """
        Do the indexing

        """

        # Save some state
        self.params = params
        self.experiments = experiments
        self.reflections = reflections

        # Do the processing
        self._index()
        self._predict()
        self._filter_reflections_based_on_centroid_distance()

    def _index(self):
        """
        Index the strong spots

        """

        # Get some stuff from experiment
        A = matrix.sqr(self.experiments[0].crystal.get_A())
        s0 = matrix.col(self.experiments[0].beam.get_s0())
        detector = self.experiments[0].detector

        # Create array if necessary
        if "miller_index" not in self.reflections:
            self.reflections["miller_index"] = flex.miller_index(len(self.reflections))

        # Index all the reflections
        xyz_list = self.reflections["xyzobs.px.value"]
        miller_index = self.reflections["miller_index"]
        selection = flex.size_t()
        num_reindexed = 0
        for i in range(len(self.reflections)):

            # Get the observed pixel coordinate
            x, y, _ = xyz_list[i]

            # Get the lab coord
            s1 = (
                matrix.col(detector[0].get_pixel_lab_coord((x, y))).normalize()
                * s0.length()
            )

            # Get the reciprocal lattice vector
            r = s1 - s0

            # Compute the fractional miller index
            hf = A.inverse() * r

            # Compute the integer miller index
            h = matrix.col(
                (
                    int(floor(hf[0] + 0.5)),
                    int(floor(hf[1] + 0.5)),
                    int(floor(hf[2] + 0.5)),
                )
            )

            # Print warning if reindexing
            if tuple(h) != miller_index[i]:
                logger.warn(
                    "Reindexing (% 3d, % 3d, % 3d) -> (% 3d, % 3d, % 3d)"
                    % (miller_index[i] + tuple(h))
                )
                num_reindexed += 1
                miller_index[i] = h
                if self.params.indexing.fail_on_bad_index:
                    raise RuntimeError("Bad index")

            # If its not indexed as 0, 0, 0 then append
            if h != matrix.col((0, 0, 0)) and (h - hf).length() < 0.3:
                selection.append(i)

        # Print some info
        logger.info(
            "Reindexed %d/%d input reflections" % (num_reindexed, len(self.reflections))
        )
        logger.info(
            "Selected %d/%d input reflections" % (len(selection), len(self.reflections))
        )

        # Select all the indexed reflections
        self.reflections.set_flags(selection, self.reflections.flags.indexed)
        self.reflections = self.reflections.select(selection)

    def _predict(self):
        """
        Predict the position of the spots

        """

        # Get some stuff from experiment
        A = matrix.sqr(self.experiments[0].crystal.get_A())
        s0 = matrix.col(self.experiments[0].beam.get_s0())

        # Compute the vector to the reciprocal lattice point
        # since this is not on the ewald sphere, lets call it s2
        h = self.reflections["miller_index"]
        s1 = flex.vec3_double(len(h))
        s2 = flex.vec3_double(len(h))
        for i in range(len(self.reflections)):
            r = A * matrix.col(h[i])
            s2[i] = s0 + r
            s1[i] = matrix.col(s2[i]).normalize() * s0.length()
        self.reflections["s1"] = s1
        self.reflections["s2"] = s2
        self.reflections["entering"] = flex.bool(len(h), False)

        # Compute the ray intersections
        xyzpx = flex.vec3_double()
        xyzmm = flex.vec3_double()
        for i in range(len(s2)):
            ss = s1[i]
            mm = self.experiments[0].detector[0].get_ray_intersection(ss)
            px = self.experiments[0].detector[0].millimeter_to_pixel(mm)
            xyzpx.append(px + (0,))
            xyzmm.append(mm + (0,))
        self.reflections["xyzcal.mm"] = xyzmm
        self.reflections["xyzcal.px"] = xyzpx
        logger.info("Do prediction for %d reflections" % len(self.reflections))

    def _filter_reflections_based_on_centroid_distance(self):
        """
        Filter reflections too far from predicted position

        """

        # Compute the x and y residuals
        Xobs, Yobs, _ = self.reflections["xyzobs.px.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.px"].parts()
        Xres = Xobs - Xcal
        Yres = Yobs - Ycal

        # Compute the epsilon residual
        s0_length = 1.0 / self.experiments[0].beam.get_wavelength()
        s1x, s1y, s1z = self.reflections["s2"].parts()
        s1_length = flex.sqrt(s1x ** 2 + s1y ** 2 + s1z ** 2)
        Eres = s1_length - s0_length

        # Initialise the fast_mcd outlier algorithm
        # fast_mcd = FastMCD((Xres, Yres, Eres))
        fast_mcd = FastMCD((Xres, Yres))

        # get location and MCD scatter estimate
        T, S = fast_mcd.get_corrected_T_and_S()

        # get squared Mahalanobis distances
        # d2s = maha_dist_sq((Xres, Yres, Eres), T, S)
        d2s = maha_dist_sq((Xres, Yres), T, S)

        # Compute the cutoff
        mahasq_cutoff = chisq_quantile(2, self.params.refinement.outlier_probability)

        # compare to the threshold and select reflections
        selection1 = d2s < mahasq_cutoff
        selection2 = (
            flex.sqrt(Xres ** 2 + Yres ** 2) < self.params.refinement.max_separation
        )
        selection = selection1 & selection2
        self.reflections = self.reflections.select(selection)

        # Print some stuff
        logger.info("-" * 80)
        logger.info("Centroid outlier rejection")
        logger.info(
            " Using MCD algorithm with probability = %f"
            % self.params.refinement.outlier_probability
        )
        logger.info(" Max X residual: %f" % flex.max(flex.abs(Xres)))
        logger.info(" Max Y residual: %f" % flex.max(flex.abs(Yres)))
        logger.info(" Max E residual: %f" % flex.max(flex.abs(Eres)))
        logger.info(" Mean X RMSD: %f" % (sqrt(flex.sum(Xres ** 2) / len(Xres))))
        logger.info(" Mean Y RMSD: %f" % (sqrt(flex.sum(Yres ** 2) / len(Yres))))
        logger.info(" Mean E RMSD: %f" % (sqrt(flex.sum(Eres ** 2) / len(Eres))))
        logger.info(" MCD location estimate: %.4f, %.4f" % tuple(T))
        logger.info(
            """ MCD scatter estimate:
      %.7f, %.7f,
      %.7f, %.7f"""
            % tuple(list(S))
        )
        # logger.info(" MCD location estimate: %.4f, %.4f, %.4f" % tuple(T))
        # logger.info(''' MCD scatter estimate:
        #   %.7f, %.7f, %.7f,
        #   %.7f, %.7f, %.7f,
        #   %.7f, %.7f, %.7f''' % tuple(list(S)))
        logger.info(" Number of outliers: %d" % selection1.count(False))
        logger.info(
            " Number of reflections with residual > %0.2f pixels: %d"
            % (self.params.refinement.max_separation, selection2.count(False))
        )
        logger.info(
            " Number of reflections selection for refinement: %d"
            % len(self.reflections)
        )
        logger.info("-" * 80)

        # Throw exception
        if len(self.reflections) < self.params.refinement.min_n_reflections:
            raise RuntimeError(
                "Too few reflections to perform refinement: got %d, expected %d"
                % (len(self.reflections), self.params.refinement.min_n_reflections)
            )


class InitialIntegrator(object):
    """
    A class to do an initial integration of strong spots

    """

    def __init__(self, params, experiments, reflections):
        """
        Do the initial integration

        """

        # Save the experiments and reflections
        self.experiments = experiments
        self.reflections = reflections
        sel = self.reflections.get_flags(self.reflections.flags.strong)
        self.reflections = self.reflections.select(sel)

        # Save the old shoeboxes
        self.shoeboxes = self.reflections["shoebox"]

        # Do the processing
        self._compute_sigma_d()
        self._compute_beam_vector()
        self._compute_bbox()
        self._allocate_shoebox()
        self._extract_shoebox()
        self._compute_mask()
        self._compute_background()
        self._compute_intensity()
        self._compute_centroid()

        self.reflections.compute_d(self.experiments)

        # Output some strong spots
        if params.debug.output.strong_spots:
            self.reflections.as_pickle("debug.strong.pickle")

        # Print shoeboxes
        if params.debug.output.print_shoeboxes:
            self._print_shoeboxes()

    def _compute_sigma_d(self):
        """
        Compute and initial spot size estimate

        """
        logger.info(
            "Computing initial sigma d estimate for %d reflections"
            % len(self.reflections)
        )
        compute_sigma_d = ComputeEsdBeamDivergence(
            self.experiments[0].detector, self.reflections
        )
        self.sigma_d = compute_sigma_d.sigma()
        logger.info("Sigma D: %.5f degrees" % (self.sigma_d * 180 / pi))
        logger.info("")

    def _compute_beam_vector(self):
        """
        Compute the obseved beam vector

        """
        panel = self.experiments[0].detector[0]
        xyz = self.reflections["xyzobs.px.value"]
        s1_obs = flex.vec3_double(len(self.reflections))
        for i in range(len(s1_obs)):
            x, y, z = xyz[i]
            s1_obs[i] = panel.get_pixel_lab_coord((x, y))
        self.reflections["s1_obs"] = s1_obs

    def _compute_bbox(self):
        """
        Compute the bounding box

        """

        logger.info(
            "Computing the bounding box for %d reflections" % len(self.reflections)
        )

        # Initialise the bounding box calculator
        compute_bbox = BBoxCalculator(
            self.experiments[0].crystal,
            self.experiments[0].beam,
            self.experiments[0].detector,
            self.experiments[0].goniometer,
            self.experiments[0].scan,
            self.sigma_d * 6,
            0,
        )

        # Compute the bounding box
        bbox = compute_bbox(
            self.reflections["s1_obs"],
            self.reflections["xyzcal.px"].parts()[2],
            self.reflections["panel"],
        )

        # Set in the reflection table
        self.reflections["bbox"] = bbox

    def _allocate_shoebox(self):
        """
        Allocate the shoebox

        """
        self.reflections["shoebox"] = flex.shoebox(
            self.reflections["panel"], self.reflections["bbox"], allocate=True
        )

    def _compute_mask(self):
        """
        Compute the spot mask

        """
        logger.info(
            "Creating the foreground mask for %d reflections" % len(self.reflections)
        )

        # Initialise the mask calculator
        mask_foreground = MaskCalculator(
            self.experiments[0].crystal,
            self.experiments[0].beam,
            self.experiments[0].detector,
            self.experiments[0].goniometer,
            self.experiments[0].scan,
            self.sigma_d * 3,
            0,
        )

        # Compute the reflection mask
        mask_foreground(
            self.reflections["shoebox"],
            self.reflections["s1_obs"],
            self.reflections["xyzcal.px"].parts()[2],
            self.reflections["panel"],
        )

        # Apply strong spot mask
        assert len(self.reflections) == len(self.shoeboxes)
        new_shoeboxes = self.reflections["shoebox"]
        old_shoeboxes = self.shoeboxes
        for s in range(len(new_shoeboxes)):
            bbox_old = old_shoeboxes[s].bbox
            mask_old = old_shoeboxes[s].mask
            bbox_new = new_shoeboxes[s].bbox
            mask_new = new_shoeboxes[s].mask
            for j in range(mask_old.all()[1]):
                for i in range(mask_old.all()[2]):
                    ii = bbox_old[0] + i - bbox_new[0]
                    jj = bbox_old[2] + j - bbox_new[2]
                    if mask_old[0, j, i] == 5:
                        if (
                            ii >= 0
                            and jj >= 0
                            and jj < mask_new.all()[1]
                            and ii < mask_new.all()[2]
                        ):
                            assert mask_new[0, jj, ii] & (1 << 0)
                            mask_new[0, jj, ii] |= (1 << 2) | (1 << 3)

    def _extract_shoebox(self):
        """
        Extract the shoebox

        """
        logger.info(
            "Extracting shoebox from image for %d reflections" % len(self.reflections)
        )
        self.reflections.extract_shoeboxes(self.experiments[0].imageset)

    def _compute_background(self):
        """
        Compute the reflection background

        """
        logger.info("Computing background for %d reflections" % len(self.reflections))
        self.reflections.compute_background(self.experiments)

    def _compute_intensity(self):
        """
        Compute the reflection intensity

        """
        logger.info("Computing intensity for %d reflections" % len(self.reflections))
        self.reflections.compute_summed_intensity()
        logger.info(
            "%d reflections integrated"
            % self.reflections.get_flags(self.reflections.flags.integrated_sum).count(
                True
            )
        )

    def _compute_centroid(self):
        """
        Compute the reflection centroid

        """
        logger.info("Computing centroid for %d reflections" % len(self.reflections))
        self.reflections.compute_centroid(self.experiments)

    def _print_shoeboxes(self):
        """
        Print the shoeboxes

        """
        sbox = self.reflections["shoebox"]
        for r in range(len(sbox)):
            data = sbox[r].data
            mask = sbox[r].mask
            logger.info(mask.as_numpy_array())
            logger.info(data.as_numpy_array())


class Refiner(object):
    """
    A class to do the refinement

    """

    def __init__(self, experiments, reflections, sigma_d, params):
        """
        Initialise the refiner

        """

        # Save the experiments and reflections
        self.params = params
        self.experiments = experiments
        self.reflections = reflections
        self.sigma_d = sigma_d
        self.plots_data = OrderedDict()

        # Set the M params
        if not hasattr(self.experiments[0].crystal, "mosaicity"):
            self.profile = ProfileModelFactory.from_sigma_d(self.params, sigma_d)
        else:
            self.profile = self.experiments[0].crystal.mosaicity

        # Preprocess the reflections
        self._preprocess()

        # Do the refinement
        for i in range(self.params.refinement.n_cycles):
            self._refine_profile()
            self._refine_crystal()

        # Post process the reflections
        self._postprocess()

    def _preprocess(self):
        """
        Preprocess the reflections

        """

        # Make some plots
        if self.params.output.html:
            self._plot_distance_from_ewald_sphere("Initial")

        # Construct the profile refiner data
        self._refiner_data = RefinerData.from_reflections(
            self.experiments[0], self.reflections
        )

    def _postprocess(self):
        """
        Postprocess the reflections

        """

        # Update predictions
        self._predict()

        # Compute prob
        self._compute_prediction_probability()

        # Save the profile model
        if self.params.debug.output.profile_model:
            self._save_profile_model()

        # Save the history
        if self.params.debug.output.history:
            self._save_history()

        # Make some plots
        if self.params.output.html:
            self._plot_distance_from_ewald_sphere("Final")

    def _refine_profile(self):
        """
        Do the profile refinement

        """
        logger.info("\n" + "=" * 80 + "\nRefining profile parmameters")

        # Create the parameterisation
        state = ModelState(
            self.experiments[0],
            self.profile.parameterisation(),
            fix_orientation=True,
            fix_unit_cell=True,
            fix_wavelength_spread=self.params.profile.wavelength_spread.model
            == "delta",
        )

        # Create the refiner and refine
        refiner = ProfileRefiner(state, self._refiner_data)
        refiner.refine()

        # Save the history
        self.history = refiner.history

        # Set the profile parameters
        self.profile.update_model(state)

        # Set the mosaicity
        self.experiments[0].crystal.mosaicity = self.profile

        # Plot the corrgram
        if self.params.debug.output.plots:
            self._plot_corrgram(refiner.correlation(), refiner.labels())

    def _refine_crystal(self):
        """
        Do the crystal refinement

        """
        if (
            self.params.profile.unit_cell.fixed is True
            and self.params.profile.orientation.fixed is True
        ):
            return

        logger.info("\n" + "=" * 80 + "\nRefining crystal parmameters")

        # Create the parameterisation
        state = ModelState(
            self.experiments[0],
            self.profile.parameterisation(),
            fix_mosaic_spread=True,
            fix_unit_cell=self.params.profile.unit_cell.fixed,
            fix_orientation=self.params.profile.orientation.fixed,
            fix_wavelength_spread=self.params.profile.wavelength_spread.model
            == "delta",
        )

        # Create the refiner and refine
        refiner = ProfileRefiner(state, self._refiner_data)
        refiner.refine()

    def _predict(self):
        """
        Predict the position of the spots

        """

        # Get some stuff from experiment
        A = matrix.sqr(self.experiments[0].crystal.get_A())
        s0 = matrix.col(self.experiments[0].beam.get_s0())

        # Compute the vector to the reciprocal lattice point
        # since this is not on the ewald sphere, lets call it s2
        h = self.reflections["miller_index"]
        s1 = flex.vec3_double(len(h))
        s2 = flex.vec3_double(len(h))
        for i in range(len(self.reflections)):
            r = A * matrix.col(h[i])
            s2[i] = s0 + r
            s1[i] = matrix.col(s2[i]).normalize() * s0.length()
        self.reflections["s1"] = s1
        self.reflections["s2"] = s2
        self.reflections["entering"] = flex.bool(len(h), False)

        # Compute the ray intersections
        xyzpx = flex.vec3_double()
        xyzmm = flex.vec3_double()
        for i in range(len(s2)):
            ss = s1[i]
            mm = self.experiments[0].detector[0].get_ray_intersection(ss)
            px = self.experiments[0].detector[0].millimeter_to_pixel(mm)
            xyzpx.append(px + (0,))
            xyzmm.append(mm + (0,))
        self.reflections["xyzcal.mm"] = xyzmm
        self.reflections["xyzcal.px"] = xyzpx
        logger.info("Do prediction for %d reflections" % len(self.reflections))

    def _compute_prediction_probability(self):

        # Get stuff from experiment
        s0 = matrix.col(self.experiments[0].beam.get_s0())
        profile = self.experiments[0].crystal.mosaicity

        # Loop through reflections
        min_p = None
        for i in range(len(self.reflections)):
            s2 = matrix.col(self.reflections[i]["s2"])
            s3 = s2.normalize() * s0.length()
            r = s2 - s0
            epsilon = s3 - s2
            sigma = profile.sigma_for_reflection(s0, r)
            sigma_inv = sigma.inverse()
            d = (epsilon.transpose() * sigma_inv * epsilon)[0]
            p = chisq_pdf(3, d)
            if min_p is None or p < min_p:
                min_p = p

        # Print some stuff
        logger.info(
            "Quantile required to predicted all observed reflections = %.5f"
            % (1 - min_p)
        )

    def _plot_distance_from_ewald_sphere(self, prefix):
        """
        Plot distance from Ewald sphere

        """

        s0 = matrix.col(self.experiments[0].beam.get_s0())
        s2 = self.reflections["s2"]
        D = flex.double(s0.length() - matrix.col(s).length() for s in s2)
        Dmean = flex.sum(D) / len(D)
        Dvar = flex.sum(flex.double([(d - Dmean) ** 2 for d in D])) / len(D)
        hist, bin_edges = np.histogram(
            D,
            bins=max(5, min(int(0.2 * len(s2)), 20)),
        )
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        self.plots_data.update(
            {
                f"{prefix}_epsilon_distribution": {
                    "data": [
                        (
                            {
                                "x": bin_centers.tolist(),
                                "y": hist.tolist(),
                                "type": "bar",
                            }
                        )
                    ],
                    "layout": {
                        "title": f"{prefix} epsilon distribution. <br>Mean(epsilon) = {Dmean:.2e}, Variance(epsilon) = {Dvar:.2e}",
                        "xaxis": {"title": "Distance from Ewald sphere (epsilon)"},
                        "yaxis": {"title": "Frequency"},
                        "bargap": 0,
                    },
                }
            }
        )

    def _plot_corrgram(self, corrmat, labels):
        """
        Plot a corrgram of correlations between parameters

        """
        plt = corrgram(corrmat, labels)
        plt.savefig("corrgram.png", dpi=300)
        plt.clf()

    def _save_profile_model(self):
        """
        Save the profile model to file

        """
        with open("profile_model.json", "w") as outfile:
            data = {
                "rlp_mosaicity": tuple(self.experiments[0].crystal.mosaicity.sigma())
            }
            json.dump(data, outfile, indent=2)

    def _save_history(self):
        """
        Save the history

        """
        with open("history.json", "w") as outfile:
            json.dump(self.history, outfile, indent=2)


class FinalIntegrator(object):
    """
    Do the final refinement

    """

    def __init__(self, params, experiments, reflections, sigma_d):
        """
        Initialise the refiner

        """

        # Save some stuff
        self.params = params
        self.experiments = experiments
        self.reflections = reflections
        self.sigma_d = sigma_d
        self.plots_data = OrderedDict()

        logger.info("\n" + "=" * 80 + "\nIntegrating reflections")

        # Do the processing
        self._compute_bbox()
        self._allocate_shoebox()
        self._extract_shoebox()
        self._compute_mask()
        self._compute_background()
        self._compute_intensity()
        self._compute_centroid()
        self._compute_partiality()

        self._print_report()
        # Plot the partialities
        if params.output.html:
            self._plot_partiality()

    def _print_report(self):

        from dials.algorithms.integration.report import IntegrationReport

        report = IntegrationReport(self.experiments, self.reflections)
        logger.info(report.tables[2].as_str())

    def _compute_bbox(self):
        """
        Do crude bbox calculation from sigma_b or from model

        """
        if self.params.integration.use_crude_shoebox_mask:
            self._compute_bbox_from_sigma_d()
        else:
            self._compute_bbox_from_model()

    def _compute_bbox_from_sigma_d(self):
        """
        Compute the bounding box

        """

        logger.info(
            "Computing the bounding box for %d reflections" % len(self.reflections)
        )

        # Initialise the bounding box calculator
        compute_bbox = BBoxCalculator(
            self.experiments[0].crystal,
            self.experiments[0].beam,
            self.experiments[0].detector,
            self.experiments[0].goniometer,
            self.experiments[0].scan,
            self.sigma_d * 6,
            0,
        )

        # Compute the bounding box
        bbox = compute_bbox(
            self.reflections["s1"],
            self.reflections["xyzcal.px"].parts()[2],
            self.reflections["panel"],
        )

        # Set in the reflection table
        self.reflections["bbox"] = bbox

        # Select reflections within detector
        x0, x1, y0, y1, _, _ = self.reflections["bbox"].parts()
        xsize, ysize = self.experiments[0].detector[0].get_image_size()
        selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)
        self.reflections = self.reflections.select(selection)
        logger.info("Filtered reflections with bbox outside image range")
        logger.info("Kept %d reflections" % len(self.reflections))

    def _compute_bbox_from_model(self):
        """
        Compute the bounding box

        """
        # Compute the bounding box
        profile = self.experiments[0].crystal.mosaicity
        profile.compute_bbox(
            self.experiments,
            self.reflections,
            self.params.integration.shoebox.probability,
        )

        # Select reflections within detector
        x0, x1, y0, y1, _, _ = self.reflections["bbox"].parts()
        xsize, ysize = self.experiments[0].detector[0].get_image_size()
        selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)
        self.reflections = self.reflections.select(selection)
        logger.info("Filtered reflections with bbox outside image range")
        logger.info("Kept %d reflections" % len(self.reflections))

    def _allocate_shoebox(self):
        """
        Allocate the shoebox

        """
        self.reflections["shoebox"] = flex.shoebox(
            self.reflections["panel"], self.reflections["bbox"], allocate=True
        )

    def _compute_mask(self):
        """
        Compute the reflection mask

        """
        if self.params.integration.use_crude_shoebox_mask:
            self._compute_mask_from_sigma_d()
        else:
            self._compute_mask_from_model()

    def _compute_mask_from_sigma_d(self):
        """
        Compute the spot mask

        """
        logger.info(
            "Creating the foreground mask for %d reflections" % len(self.reflections)
        )

        # Initialise the mask calculator
        mask_foreground = MaskCalculator(
            self.experiments[0].crystal,
            self.experiments[0].beam,
            self.experiments[0].detector,
            self.experiments[0].goniometer,
            self.experiments[0].scan,
            self.sigma_d * 3,
            0,
        )

        # Compute the reflection mask
        mask_foreground(
            self.reflections["shoebox"],
            self.reflections["s1"],
            self.reflections["xyzcal.px"].parts()[2],
            self.reflections["panel"],
        )

    def _compute_mask_from_model(self):
        """
        Compute the reflection mask

        """
        profile = self.experiments[0].crystal.mosaicity
        profile.compute_mask(
            self.experiments,
            self.reflections,
            self.params.integration.shoebox.probability,
        )

    def _extract_shoebox(self):
        """
        Extract the shoebox

        """
        logger.info(
            "Extracting shoebox from image for %d reflections" % len(self.reflections)
        )
        self.reflections.extract_shoeboxes(self.experiments[0].imageset)

    def _compute_background(self):
        """
        Compute the reflection background

        """
        logger.info("Computing background for %d reflections" % len(self.reflections))
        self.reflections.compute_background(self.experiments)

    def _compute_intensity(self):
        """
        Compute the reflection intensity

        """
        logger.info("Computing intensity for %d reflections" % len(self.reflections))
        self.reflections.compute_summed_intensity()
        self.reflections.compute_corrections(self.experiments)
        logger.info(
            "%d reflections integrated"
            % self.reflections.get_flags(self.reflections.flags.integrated_sum).count(
                True
            )
        )

    def _compute_centroid(self):
        """
        Compute the reflection centroid

        """
        logger.info("Computing centroid for %d reflections" % len(self.reflections))
        self.reflections.compute_centroid(self.experiments)

    def _compute_partiality(self):
        """
        Compute the partiality

        """
        profile = self.experiments[0].crystal.mosaicity
        profile.compute_partiality(self.experiments, self.reflections)

    def _plot_partiality(self):
        """
        Plot the partiality

        """

        hist, bin_edges = np.histogram(
            self.reflections["partiality"],
            bins=max(5, min(int(0.2 * self.reflections.size()), 20)),
        )
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        self.plots_data.update(
            {
                "partiality_distribution": {
                    "data": [
                        (
                            {
                                "x": bin_centers.tolist(),
                                "y": hist.tolist(),
                                "type": "bar",
                            }
                        )
                    ],
                    "layout": {
                        "title": "Partiality distribution",
                        "xaxis": {"title": "Partiality"},
                        "yaxis": {"title": "Frequency"},
                        "bargap": 0,
                    },
                }
            }
        )


class Integrator(object):
    """
    Class to perform integration of stills in the following way:

    1. Do an initial integration of strong reflections
    2. Refine profile and crystal parameters
    3. Do a final integration of all reflections

    """

    def __init__(self, experiments, reflections, params=None):
        """
        Initialise the integrator

        """

        # Only use single experiment at the moment
        if len(experiments) > 1:
            raise RuntimeError("Only 1 experiment can be processed")

        # Set the parameters
        if params is not None:
            self.params = params
        else:
            self.params = phil_scope.extract(parse(""))

        # Save some stuff
        self.experiments = experiments
        self.strong = reflections
        self.reference = None
        self.reflections = None
        self.sigma_d = None
        self.plots_data = OrderedDict()

    def reindex_strong_spots(self):
        """
        Reindex the strong spot list

        """
        indexer = Indexer(self.params, self.experiments, self.strong)
        self.reference = indexer.reflections

    def integrate_strong_spots(self):
        """
        Do an initial integration of the strong spots

        """
        integrator = InitialIntegrator(self.params, self.experiments, self.reference)
        self.reference = integrator.reflections
        self.sigma_d = integrator.sigma_d

    def refine(self):
        """
        Do the refinement of profile and crystal parameters

        """
        refiner = Refiner(self.experiments, self.reference, self.sigma_d, self.params)
        self.experiments = refiner.experiments
        self.reference = refiner.reflections
        self.plots_data.update(refiner.plots_data)

    def predict(self):
        """
        Predict the reflections

        """
        logger.info("\n" + "=" * 80 + "\nPredicting reflections")

        # Set a resolution range
        if self.params.prediction.d_min is None:
            s0 = self.experiments[0].beam.get_s0()
            d_min = self.experiments[0].detector.get_max_resolution(s0)
        else:
            d_min = self.params.predictions.d_min

        # Create the index generator
        index_generator = IndexGenerator(
            self.experiments[0].crystal.get_unit_cell(),
            self.experiments[0].crystal.get_space_group().type(),
            d_min,
        )

        # Get an array of miller indices
        miller_indices_to_test = index_generator.to_array()
        logger.info("Generated %d miller indices" % len(miller_indices_to_test))

        # Get the covariance matrix
        profile = self.experiments[0].crystal.mosaicity
        id_map = dict(self.strong.experiment_identifiers())
        self.reflections = profile.predict_reflections(
            self.experiments, miller_indices_to_test, self.params.prediction.probability
        )

        # Do the prediction
        self.reference = self.reference
        self.reflections.compute_d(self.experiments)
        logger.info("Predicted %d reflections" % len(self.reflections))

        _, _, unmatched = self.reflections.match_with_reference(self.reference)

        # now set the identifiers
        ids_ = set(self.reflections["id"])
        assert ids_ == set(id_map.keys()), f"{ids_}, {id_map.keys()}"
        for id_ in ids_:
            self.reflections.experiment_identifiers()[id_] = id_map[id_]

        # Add unmatched
        # columns = flex.std_string()
        # for col in unmatched.keys():
        #   if col in self.reflections:
        #     columns.append(col)
        # unmatched = unmatched.select(columns)
        # unmatched['id'] = flex.size_t(list(unmatched['id']))
        # self.reflections.extend(unmatched)

    def integrate(self):
        """
        Do an final integration of the reflections

        """
        integrator = FinalIntegrator(
            self.params, self.experiments, self.reflections, self.sigma_d
        )
        self.reflections = integrator.reflections
        self.plots_data.update(integrator.plots_data)

        # Delete shoeboxes if necessary
        if not self.params.debug.output.shoeboxes:
            del self.reflections["shoebox"]

        # Delete corrections if specified
        if not self.params.integration.corrections.lp:
            del self.reflections["lp"]
        if not self.params.integration.corrections.dqe:
            del self.reflections["dqe"]
        if not self.params.integration.corrections.partiality:
            del self.reflections["partiality"]
            del self.reflections["partiality.inv.variance"]


def generate_html_report(integrator, filename):
    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("simple_report.html")
    html = template.render(
        page_title="DIALS SSX integration report",
        panel_title="Integration plots",
        panel_id="ewald",
        graphs=integrator.plots_data,
    )
    logger.info(f"Writing html report to {filename}")
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))
