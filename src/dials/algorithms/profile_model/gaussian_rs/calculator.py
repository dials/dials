# FIXME Mosaicity seems to be overestimated
# FIXME Don't know how XDS REFLECTING_RANGE is calculated
# FIXME Don't know what XDS REFLECTION_RANGE is used for
# FIXME Don't know what XDS BEAM_DIVERGENCE is used for
# FIXME Should maybe be scan varying
# FIXME Don't know how XDS calculated the n_sigma


from __future__ import annotations

import collections
import copy
import logging
import math

import numpy as np

import scitbx.math
from dxtbx import flumpy

from dials.array_family import flex

logger = logging.getLogger(__name__)


class ComputeEsdBeamDivergence:
    """Calculate the E.s.d of the beam divergence."""

    def __init__(self, detector, reflections, centroid_definition="s1"):
        """Calculate the E.s.d of the beam divergence.

        Params:
            detector The detector class
            reflections The reflections
            centroid_definition ENUM com or s1
        """
        # Calculate the beam direction variances
        variance = self._beam_direction_variance_list(
            detector, reflections, centroid_definition
        )

        # Calculate and return the e.s.d of the beam divergence
        self._sigma = math.sqrt(np.sum(variance) / variance.size)

    def sigma(self):
        """Return the E.S.D of the beam divergence."""
        return self._sigma

    def _beam_direction_variance_list(
        self, detector, reflections, centroid_definition="s1"
    ):
        """Calculate the variance in beam direction for each spot.

        Params:
            detector The detector model
            reflections The list of reflections
            centroid_definition ENUM com or s1

        Returns:
            The list of variances
        """
        # Get the reflection columns
        shoebox = reflections["shoebox"]
        xyz = reflections["xyzobs.px.value"]

        # Loop through all the reflections
        variance = np.array([], dtype=np.float64)

        if centroid_definition == "com":
            # Calculate the beam vector at the centroid
            s1_centroid = []
            for r in range(len(reflections)):
                panel = shoebox[r].panel
                s1_centroid.append(detector[panel].get_pixel_lab_coord(xyz[r][0:2]))
        else:
            s1_centroid = reflections["s1"]

        for r in range(len(reflections)):
            # Get the coordinates and values of valid shoebox pixels
            # FIXME maybe I note in Kabsch (2010) s3.1 step (v) is
            # background subtraction, appears to be missing here.
            mask = shoebox[r].mask != 0
            values = flumpy.to_numpy(shoebox[r].values(mask))
            s1 = shoebox[r].beam_vectors(detector, mask)

            angles = flumpy.to_numpy(s1.angle(s1_centroid[r], deg=False))

            if np.sum(values) > 1:
                var = np.sum(values * np.square(angles)) / (np.sum(values) - 1)
                variance = np.append(variance, var)

        # Return a list of variances
        return variance


class FractionOfObservedIntensity:
    """Calculate the fraction of observed intensity for different sigma_m."""

    def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
        """Initialise the algorithm. Calculate the list of tau and zetas.

        Params:
            reflections The list of reflections
            experiment The experiment object
        """
        # Get the oscillation width
        dphi2 = scan.get_oscillation(deg=False)[1] / 2.0

        # Calculate a list of angles and zeta's
        tau, zeta = self._calculate_tau_and_zeta(
            crystal, beam, detector, goniometer, scan, reflections
        )

        # Calculate zeta * (tau +- dphi / 2) / math.sqrt(2)
        self.e1 = (tau + dphi2) * flex.abs(zeta) / math.sqrt(2.0)
        self.e2 = (tau - dphi2) * flex.abs(zeta) / math.sqrt(2.0)

    def _calculate_tau_and_zeta(
        self, crystal, beam, detector, goniometer, scan, reflections
    ):
        """Calculate the list of tau and zeta needed for the calculation.

        Params:
            reflections The list of reflections
            experiment The experiment object.

        Returns:
            (list of tau, list of zeta)
        """
        from dials.algorithms.shoebox import MaskCode

        mask_code = MaskCode.Valid | MaskCode.Foreground

        # Calculate the list of frames and z coords
        sbox = reflections["shoebox"]
        bbox = reflections["bbox"]
        phi = reflections["xyzcal.mm"].parts()[2]

        # Calculate the zeta list
        zeta = reflections["zeta"]

        # Calculate the list of tau values
        tau = []
        zeta2 = []
        for s, b, p, z in zip(sbox, bbox, phi, zeta):
            for z0, f in enumerate(range(b[4], b[5])):
                phi0 = scan.get_angle_from_array_index(int(f), deg=False)
                phi1 = scan.get_angle_from_array_index(int(f) + 1, deg=False)
                m = s.mask[z0 : z0 + 1, :, :]
                if m.count(mask_code) > 0:
                    tau.append((phi1 + phi0) / 2.0 - p)
                    zeta2.append(z)

        # Return the list of tau and zeta
        return flex.double(tau), flex.double(zeta2)

    def __call__(self, sigma_m):
        """Calculate the fraction of observed intensity for each observation.

        Params:
            sigma_m The mosaicity

        Returns:
            A list of log intensity fractions
        """
        # Tiny value
        TINY = 1e-10
        assert sigma_m > TINY

        # Calculate the two components to the fraction
        a = scitbx.math.erf(self.e1 / sigma_m)
        b = scitbx.math.erf(self.e2 / sigma_m)

        # Calculate the fraction of observed reflection intensity
        R = (a - b) / 2.0

        # Set any points <= 0 to 1e-10 (otherwise will get a floating
        # point error in log calculation below).
        assert R.all_ge(0)
        mask = R < TINY
        assert mask.count(True) < len(mask)
        R.set_selected(mask, TINY)

        # Return the logarithm of r
        return flex.log(R)


class ComputeEsdReflectingRange:
    """calculate the e.s.d of the reflecting range (mosaicity)."""

    class Estimator:
        """Estimate E.s.d reflecting range by maximum likelihood estimation."""

        def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
            """Initialise the optimization."""
            from scitbx import simplex

            # FIXME in here this code is very unstable or actually broken if
            # we pass in a few lone images i.e. screening shots - propose need
            # for alternative algorithm, in meantime can avoid code death if
            # something like this
            #
            # if scan.get_num_images() == 1:
            #   self.sigma = 0.5 * scan.get_oscillation_range()[1] * math.pi / 180.0
            #   return
            #
            # is used... @JMP please could we discuss? assert best method to get
            # sigma_m in that case is actually to look at present / absent
            # reflections as per how Mosflm does it...
            # Initialise the function used in likelihood estimation.
            self._R = FractionOfObservedIntensity(
                crystal, beam, detector, goniometer, scan, reflections
            )

            # Set the starting values to try 1, 3 degrees seems sensible for
            # crystal mosaic spread
            start = 1 * math.pi / 180
            stop = 3 * math.pi / 180
            starting_simplex = [flex.double([start]), flex.double([stop])]

            # Initialise the optimizer
            optimizer = simplex.simplex_opt(
                1, matrix=starting_simplex, evaluator=self, tolerance=1e-7
            )

            # Get the solution
            self.sigma = math.exp(optimizer.get_solution()[0])

        def target(self, log_sigma):
            """The target for minimization."""
            return -flex.sum(self._R(math.exp(log_sigma[0])))

    class CrudeEstimator:
        """If the main estimator failed make a crude estimate"""

        def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
            # Calculate a list of angles and zeta's
            tau, zeta = self._calculate_tau_and_zeta(
                crystal, beam, detector, goniometer, scan, reflections
            )

            # Calculate zeta * (tau +- dphi / 2) / math.sqrt(2)
            X = tau * zeta
            mv = flex.mean_and_variance(X)
            self.sigma = math.sqrt(mv.unweighted_sample_variance())

        def _calculate_tau_and_zeta(
            self, crystal, beam, detector, goniometer, scan, reflections
        ):
            """Calculate the list of tau and zeta needed for the calculation.

            Params:
                reflections The list of reflections
                experiment The experiment object.

            Returns:
                (list of tau, list of zeta)
            """
            from dials.algorithms.shoebox import MaskCode

            mask_code = MaskCode.Valid | MaskCode.Foreground

            # Calculate the list of frames and z coords
            sbox = reflections["shoebox"]
            bbox = reflections["bbox"]
            phi = reflections["xyzcal.mm"].parts()[2]

            # Calculate the zeta list
            zeta = reflections["zeta"]

            # Calculate the list of tau values
            tau = []
            zeta2 = []
            for s, b, p, z in zip(sbox, bbox, phi, zeta):
                for z0, f in enumerate(range(b[4], b[5])):
                    phi0 = scan.get_angle_from_array_index(int(f), deg=False)
                    phi1 = scan.get_angle_from_array_index(int(f) + 1, deg=False)
                    m = s.mask[z0 : z0 + 1, :, :]
                    if m.count(mask_code) > 0:
                        tau.append((phi1 + phi0) / 2.0 - p)
                        zeta2.append(z)

            # Return the list of tau and zeta
            return flex.double(tau), flex.double(zeta2)

    class ExtendedEstimator:
        """Try to estimate using knowledge of intensities"""

        def __init__(
            self,
            crystal,
            beam,
            detector,
            goniometer,
            scan,
            reflections,
            n_macro_cycles=10,
        ):
            from scitbx import simplex

            # Get the oscillation width
            dphi2 = scan.get_oscillation(deg=False)[1] / 2.0

            # Calculate a list of angles and zeta's
            tau, zeta, n, indices = self._calculate_tau_and_zeta(
                crystal, beam, detector, goniometer, scan, reflections
            )

            # Calculate zeta * (tau +- dphi / 2) / math.sqrt(2)
            self.e1 = (tau + dphi2) * flex.abs(zeta) / math.sqrt(2.0)
            self.e2 = (tau - dphi2) * flex.abs(zeta) / math.sqrt(2.0)
            self.indices = indices
            if len(self.e1) == 0:
                raise RuntimeError(
                    "Something went wrong. Zero pixels selected for estimation of profile parameters."
                )

            # Compute intensity
            self.K = flex.double()
            self.nj = []
            for i0, i1 in zip(self.indices[:-1], self.indices[1:]):
                nj = n[i0:i1]
                self.K.append(flex.sum(nj))
                self.nj.append(nj)

            # Set the starting values to try 1, 3 degrees seems sensible for
            # crystal mosaic spread
            start = math.log(0.1 * math.pi / 180)
            stop = math.log(1 * math.pi / 180)
            starting_simplex = [flex.double([start]), flex.double([stop])]

            # Initialise the optimizer
            optimizer = simplex.simplex_opt(
                1, matrix=starting_simplex, evaluator=self, tolerance=1e-3
            )

            # Get the solution
            sigma = math.exp(optimizer.get_solution()[0])

            # Save the result
            self.sigma = sigma

        def target(self, log_sigma):
            """The target for minimization."""
            sigma_m = math.exp(log_sigma[0])

            # Tiny value
            TINY = 1e-10
            assert sigma_m > TINY

            # Calculate the two components to the fraction
            a = scitbx.math.erf(self.e1 / sigma_m)
            b = scitbx.math.erf(self.e2 / sigma_m)

            # Calculate the fraction of observed reflection intensity
            zi = (a - b) / 2.0

            # Set any points <= 0 to 1e-10 (otherwise will get a floating
            # point error in log calculation below).
            assert zi.all_ge(0)
            mask = zi < TINY
            assert mask.count(True) < len(mask)
            zi.set_selected(mask, TINY)

            # Compute the likelihood
            #
            # The likelihood here is a result of the sum of two log likelihood
            # functions:
            #
            # The first is the same as the one in Kabsch2010 as applied to the
            # reflection as a whole. This results in the term log(Z)
            #
            # The second is the likelihood for each reflection modelling as a Poisson
            # distributions with shape given by sigma M. This gives sum(ci log(zi)) -
            # sum(ci)*log(sum(zi))
            #
            # If the reflection is recorded on 1 frame, the second component is zero
            # and so the likelihood is dominated by the first term which can be seen
            # as a prior for sigma, which accounts for which reflections were actually
            # recorded.
            #
            L = 0
            for (kj, nj, i0, i1) in zip(
                self.K, self.nj, self.indices[:-1], self.indices[1:]
            ):
                zj = zi[i0:i1]
                logZ = math.log(flex.sum(zj))
                L += flex.sum(nj * flex.log(zj)) - kj * logZ + logZ
            logger.debug("Sigma M: %f, log(L): %f", sigma_m * 180 / math.pi, L)

            # Return the logarithm of r
            return -L

        def _calculate_tau_and_zeta(
            self, crystal, beam, detector, goniometer, scan, reflections
        ):
            """Calculate the list of tau and zeta needed for the calculation.

            Params:
                reflections The list of reflections
                experiment The experiment object.

            Returns:
                (list of tau, list of zeta)
            """
            from dials.algorithms.shoebox import MaskCode

            mask_code = MaskCode.Valid | MaskCode.Foreground

            # Calculate the list of frames and z coords
            sbox = reflections["shoebox"]
            phi = reflections["xyzcal.mm"].parts()[2]

            # Calculate the zeta list
            zeta = reflections["zeta"]

            # Calculate the list of tau values
            tau = []
            zeta2 = []
            num = []
            indices = [0]
            for s, p, z in zip(sbox, phi, zeta):
                b = s.bbox
                for z0, f in enumerate(range(b[4], b[5])):
                    phi0 = scan.get_angle_from_array_index(int(f), deg=False)
                    phi1 = scan.get_angle_from_array_index(int(f) + 1, deg=False)
                    d = s.data[z0 : z0 + 1, :, :]
                    m = s.mask[z0 : z0 + 1, :, :]
                    d = flex.sum(d.as_1d().select(m.as_1d() == mask_code))
                    if d > 0:
                        tau.append((phi1 + phi0) / 2.0 - p)
                        zeta2.append(z)
                        num.append(d)
                if len(zeta2) > indices[-1]:
                    indices.append(len(zeta2))

            # Return the list of tau and zeta
            return (
                flex.double(tau),
                flex.double(zeta2),
                flex.double(num),
                flex.size_t(indices),
            )

    def __init__(
        self, crystal, beam, detector, goniometer, scan, reflections, algorithm="basic"
    ):
        """initialise the algorithm with the scan.

        params:
            scan the scan object
        """

        if algorithm == "basic":

            # Calculate sigma_m
            try:
                estimator = ComputeEsdReflectingRange.Estimator(
                    crystal, beam, detector, goniometer, scan, reflections
                )
            except Exception:
                logger.info("Using Crude Mosaicity estimator")
                estimator = ComputeEsdReflectingRange.CrudeEstimator(
                    crystal, beam, detector, goniometer, scan, reflections
                )

        elif algorithm == "extended":
            estimator = ComputeEsdReflectingRange.ExtendedEstimator(
                crystal, beam, detector, goniometer, scan, reflections
            )

        # Save the solution
        self._sigma = estimator.sigma

    def sigma(self):
        """Return the E.S.D reflecting rang."""
        return self._sigma


class ProfileModelCalculator:
    """Class to help calculate the profile model."""

    def __init__(
        self,
        reflections,
        crystal,
        beam,
        detector,
        goniometer,
        scan,
        min_zeta=0.05,
        algorithm="basic",
        centroid_definition="s1",
    ):
        """Calculate the profile model."""
        from dxtbx.model.experiment_list import Experiment

        # Check input has what we want
        assert reflections is not None
        assert "miller_index" in reflections
        assert "s1" in reflections
        assert "shoebox" in reflections
        assert "xyzobs.px.value" in reflections
        assert "xyzcal.mm" in reflections

        assert centroid_definition in ("s1", "com")

        n_all = reflections.size()

        # stills images behave differently in here
        if goniometer is None or scan is None or scan.is_still():
            logger.info("Using %d reflections for sigma calculation", n_all)
            logger.info("Calculating E.S.D Beam Divergence.")
            beam_divergence = ComputeEsdBeamDivergence(
                detector, reflections, centroid_definition
            )
            self._sigma_b = beam_divergence.sigma()
            # FIXME calculate properly
            self._sigma_m = 0.0
        else:
            # filter reflections before determining the reflection profile
            # parameters - first by used_in_refinement then by zeta

            reflections = _select_reflections_for_sigma_calc(reflections)

            zeta = reflections.compute_zeta(
                Experiment(
                    crystal=crystal,
                    beam=beam,
                    detector=detector,
                    goniometer=goniometer,
                    scan=scan,
                )
            )
            reflections = reflections.select(flex.abs(zeta) >= min_zeta)
            n_use = reflections.size()

            logger.info("Using %d / %d reflections for sigma calculation", n_use, n_all)
            logger.info("Calculating E.S.D Beam Divergence.")
            beam_divergence = ComputeEsdBeamDivergence(
                detector, reflections, centroid_definition
            )

            self._sigma_b = beam_divergence.sigma()

            logger.info("Calculating E.S.D Reflecting Range (mosaicity).")
            reflecting_range = ComputeEsdReflectingRange(
                crystal,
                beam,
                detector,
                goniometer,
                scan,
                reflections,
                algorithm=algorithm,
            )

            self._sigma_m = reflecting_range.sigma()

        logger.info(" sigma b: %f degrees", self._sigma_b * 180 / math.pi)
        logger.info(" sigma m: %f degrees", self._sigma_m * 180 / math.pi)

    def sigma_b(self):
        """Return the E.S.D beam divergence."""
        return self._sigma_b

    def sigma_m(self):
        """Return the E.S.D reflecting range."""
        return self._sigma_m


class ScanVaryingProfileModelCalculator:
    """Class to help calculate the profile model."""

    def __init__(
        self,
        reflections,
        crystal,
        beam,
        detector,
        goniometer,
        scan,
        min_zeta=0.05,
        algorithm="basic",
        centroid_definition="s1",
    ):
        """Calculate the profile model."""
        from dxtbx.model.experiment_list import Experiment

        # Check input has what we want
        assert reflections is not None
        assert "miller_index" in reflections
        assert "s1" in reflections
        assert "shoebox" in reflections
        assert "xyzobs.px.value" in reflections
        assert "xyzcal.mm" in reflections

        assert centroid_definition in ("s1", "com")

        # Select by zeta
        zeta = reflections.compute_zeta(
            Experiment(
                crystal=crystal,
                beam=beam,
                detector=detector,
                goniometer=goniometer,
                scan=scan,
            )
        )
        mask = flex.abs(zeta) >= min_zeta
        reflections = reflections.select(mask)

        # Split the reflections into partials
        reflections = copy.deepcopy(reflections)
        reflections.split_partials_with_shoebox()

        # Get a list of reflections for each frame
        bbox = reflections["bbox"]
        index_list = collections.defaultdict(list)
        for i, (x0, x1, y0, y1, z0, z1) in enumerate(bbox):
            assert z1 == z0 + 1
            index_list[z0].append(i)
        reflection_list = {
            key: reflections.select(flex.size_t(value))
            for key, value in index_list.items()
        }

        # The range of frames
        z0, z1 = scan.get_array_range()
        min_z = min(index_list)
        max_z = max(index_list)
        assert z0 == min_z
        assert z1 == max_z + 1

        # Compute for all frames
        self._num = []
        sigma_b = flex.double()
        sigma_m = flex.double()
        for i in range(z0, z1):

            # Get reflections at the index
            reflections = reflection_list[i]

            self._num.append(len(reflections))

            # Calculate the E.S.D of the beam divergence
            beam_divergence = ComputeEsdBeamDivergence(detector, reflections)

            # Set the sigma b
            sigma_b.append(beam_divergence.sigma())

            # Calculate the E.S.D of the reflecting range
            reflecting_range = ComputeEsdReflectingRange(
                crystal, beam, detector, goniometer, scan, reflections
            )

            logger.info(
                "Computing profile model for frame %d: sigma_b = %.4f degrees, sigma_m = %.4f degrees",
                i,
                beam_divergence.sigma() * 180 / math.pi,
                reflecting_range.sigma() * 180 / math.pi,
            )

            # Set the sigmas
            sigma_m.append(reflecting_range.sigma())

        def convolve(data, kernel):
            assert len(kernel) & 1
            result = []
            mid = len(kernel) // 2
            for i in range(len(data)):
                r = 0
                for j in range(len(kernel)):
                    k = i - mid + j
                    if k < 0:
                        r += kernel[j] * data[0]
                    elif k >= len(data):
                        r += kernel[j] * data[-1]
                    else:
                        r += kernel[j] * data[k]
                result.append(r)
            return result

        def gaussian_kernel(n):
            assert n & 1
            mid = n // 2
            sigma = mid / 3.0
            kernel = [math.exp(-((i - mid) ** 2) / (2 * sigma**2)) for i in range(n)]
            kernel = [k / sum(kernel) for k in kernel]
            return kernel

        # Smooth the parameters
        kernel = gaussian_kernel(51)
        sigma_b_sq_new = convolve(flex.pow2(sigma_b), kernel)
        sigma_m_sq_new = convolve(flex.pow2(sigma_m), kernel)

        # Print the output - mean as is scan varying
        mean_sigma_b = math.sqrt(sum(flex.pow2(sigma_b)) / len(sigma_b))
        mean_sigma_m = math.sqrt(sum(flex.pow2(sigma_m)) / len(sigma_m))

        # Save the smoothed parameters
        self._sigma_b = flex.sqrt(flex.double(sigma_b_sq_new))
        self._sigma_m = flex.sqrt(flex.double(sigma_m_sq_new))
        assert len(self._sigma_b) == len(self._sigma_m)

        # Print out smoothed profile parameters
        for i in range(len(sigma_b)):
            logger.info(
                "Smoothed profile model for frame %d: sigma_b = %.4f degrees, sigma_m = %.4f degrees",
                i,
                self._sigma_b[i] * 180 / math.pi,
                self._sigma_m[i] * 180 / math.pi,
            )

        # Print the mean parameters
        logger.info(" sigma b: %f degrees", mean_sigma_b * 180 / math.pi)
        logger.info(" sigma m: %f degrees", mean_sigma_m * 180 / math.pi)

    def num(self):
        """The number of reflections used."""
        return self._num

    def sigma_b(self):
        """Return the E.S.D beam divergence."""
        return self._sigma_b

    def sigma_m(self):
        """Return the E.S.D reflecting range."""
        return self._sigma_m


def _select_reflections_for_sigma_calc(reflections, min_number_of_refl=10000):
    """Determine a subset of reflections to use for sigma_m calculation."""
    n_ref = reflections.size()
    if n_ref >= min_number_of_refl:
        # ideally use well-sampled selection from refinement
        used_in_ref = reflections.get_flags(reflections.flags.used_in_refinement)
        n_used_in_ref = used_in_ref.count(True)
        if n_used_in_ref >= min_number_of_refl:
            selected_reflections = reflections.select(used_in_ref)
            logger.debug(
                "Using %s reflections with used_in_refinement flag for sigma calculation",
                n_used_in_ref,
            )

        # handle case where used_in_refl not set/data not refined
        else:
            # more than min_no of refls, but refinement didn't use more than
            # min_no: still use the ones from refinement and 'top up'
            if n_used_in_ref:
                selected_reflections = reflections.select(used_in_ref)
                other_refls = reflections.select(~used_in_ref)
            else:
                selected_reflections = flex.reflection_table()
                other_refls = reflections
            # top up with every nth reflection needed to make up the number needed
            n_sel = selected_reflections.size()
            n_other = other_refls.size()
            step_size = int(math.floor(n_other / (min_number_of_refl - n_sel)))
            sel = flex.size_t(i for i in range(0, n_other, step_size))
            selected_reflections.extend(other_refls.select(sel))
            logger.debug(
                "Using %s reflections for sigma calculation",
                selected_reflections.size(),
            )
        return selected_reflections

    logger.debug("Using all suitable reflections for sigma calculation")
    return reflections
