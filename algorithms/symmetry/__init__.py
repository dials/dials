"""Methods for symmetry determination.

This module provides a base class for symmetry determination algorithms.
"""
from __future__ import division, absolute_import, print_function

import logging

logger = logging.getLogger(__name__)

from six.moves import cStringIO as StringIO

import libtbx
from scitbx.array_family import flex
from cctbx import adptbx
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews


class symmetry_base(object):
    """Base class for symmetry analysis."""

    def __init__(
        self,
        intensities,
        normalisation="ml_aniso",
        lattice_symmetry_max_delta=2.0,
        d_min=libtbx.Auto,
        min_i_mean_over_sigma_mean=4,
        min_cc_half=0.6,
        relative_length_tolerance=None,
        absolute_angle_tolerance=None,
    ):
        """Initialise a symmetry_base object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            symmetry anaylsis.
          normalisation (str): The normalisation method to use. Possible choices are
            'kernel', 'quasi', 'ml_iso' and 'ml_aniso'. Set to None to switch off
            normalisation altogether.
          lattice_symmetry_max_delta (float): The maximum value of delta for
            determining the lattice symmetry using the algorithm of Le Page (1982).
          d_min (float): Optional resolution cutoff to be applied to the input
            intensities. If set to :data:`libtbx.Auto` then d_min will be
            automatically determined according to the parameters
            ``min_i_mean_over_sigma_mean`` and ``min_cc_half``.
          min_i_mean_over_sigma_mean (float): minimum value of |I|/|sigma(I)| for
            automatic determination of resolution cutoff.
          min_cc_half (float): minimum value of CC1/2 for automatic determination of
            resolution cutoff.
          relative_length_tolerance (float): Relative length tolerance in checking
            consistency of input unit cells against the median unit cell.
          absolute_angle_tolerance (float): Absolute angle tolerance in checking
            consistency of input unit cells against the median unit cell.

        """
        self.input_intensities = intensities

        uc_params = [flex.double() for i in range(6)]
        for d in self.input_intensities:
            for i, p in enumerate(d.unit_cell().parameters()):
                uc_params[i].append(p)
        self.median_unit_cell = uctbx.unit_cell(
            parameters=[flex.median(p) for p in uc_params]
        )
        self._check_unit_cell_consistency(
            relative_length_tolerance, absolute_angle_tolerance
        )

        self.intensities = self.input_intensities[0]
        self.dataset_ids = flex.double(self.intensities.size(), 0)
        for i, d in enumerate(self.input_intensities[1:]):
            self.intensities = self.intensities.concatenate(
                d, assert_is_similar_symmetry=False
            )
            self.dataset_ids.extend(flex.double(d.size(), i + 1))
        self.intensities = self.intensities.customized_copy(
            unit_cell=self.median_unit_cell
        )
        self.intensities.set_observation_type_xray_intensity()
        sys_absent_flags = self.intensities.sys_absent_flags(integral_only=True).data()
        self.intensities = self.intensities.select(~sys_absent_flags)
        self.dataset_ids = self.dataset_ids.select(~sys_absent_flags)

        self.lattice_symmetry_max_delta = lattice_symmetry_max_delta
        self.subgroups = metric_subgroups(
            self.intensities.crystal_symmetry(),
            max_delta=self.lattice_symmetry_max_delta,
            bravais_types_only=False,
        )

        self.cb_op_inp_min = self.subgroups.cb_op_inp_minimum
        self.intensities = (
            self.intensities.change_basis(self.cb_op_inp_min)
            .customized_copy(space_group_info=sgtbx.space_group_info("P1"))
            .map_to_asu()
            .set_info(self.intensities.info())
        )

        self.lattice_group = (
            self.subgroups.result_groups[0]["subsym"].space_group().make_tidy()
        )
        self.patterson_group = (
            self.lattice_group.build_derived_patterson_group().make_tidy()
        )
        logger.info("Patterson group: %s" % self.patterson_group.info())

        sel = self.patterson_group.epsilon(self.intensities.indices()) == 1
        self.intensities = self.intensities.select(sel)
        self.dataset_ids = self.dataset_ids.select(sel)

        # Correct SDs by "typical" SD factors
        self._correct_sigmas(sd_fac=2.0, sd_b=0.0, sd_add=0.03)

        self._resolution_filter(d_min, min_i_mean_over_sigma_mean, min_cc_half)

    def _check_unit_cell_consistency(
        self, relative_length_tolerance, absolute_angle_tolerance
    ):
        for d in self.input_intensities:
            if (
                relative_length_tolerance is not None
                and absolute_angle_tolerance is not None
            ):
                assert d.unit_cell().is_similar_to(
                    self.median_unit_cell,
                    relative_length_tolerance,
                    absolute_angle_tolerance,
                ), (str(d.unit_cell()), str(self.median_unit_cell))

    def _normalise(self):
        if normalisation is None:
            return
        elif normalisation == "kernel":
            normalise = self.kernel_normalisation
        elif normalisation == "quasi":
            normalise = self.quasi_normalisation
        elif normalisation == "ml_iso":
            normalise = self.ml_iso_normalisation
        elif normalisation == "ml_aniso":
            normalise = self.ml_aniso_normalisation

        for i in range(int(flex.max(self.dataset_ids) + 1)):
            logger.info("Normalising intensities for dataset %i" % (i + 1))
            intensities = self.intensities.select(self.dataset_ids == i)
            if i == 0:
                normalised_intensities = normalise(intensities)
            else:
                normalised_intensities = normalised_intensities.concatenate(
                    normalise(intensities)
                )
        self.intensities = normalised_intensities.set_info(
            self.intensities.info()
        ).set_observation_type_xray_intensity()

    def _correct_sigmas(self, sd_fac, sd_b, sd_add):
        # sd' = SDfac * Sqrt(sd^2 + SdB * I + (SDadd * I)^2)
        variance = flex.pow2(self.intensities.sigmas())
        si2 = flex.pow2(sd_add * self.intensities.data())
        ssc = variance + sd_b * self.intensities.data() + si2
        MINVARINFRAC = 0.1
        ssc.set_selected(ssc < MINVARINFRAC * variance, MINVARINFRAC * variance)
        sd = sd_fac * flex.sqrt(ssc)
        self.intensities = self.intensities.customized_copy(sigmas=sd).set_info(
            self.intensities.info()
        )

    @staticmethod
    def kernel_normalisation(intensities):
        """Kernel normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.

        """
        normalisation = absolute_scaling.kernel_normalisation(
            intensities, auto_kernel=True
        )
        return normalisation.normalised_miller.deep_copy().set_info(intensities.info())

    @staticmethod
    def quasi_normalisation(intensities):
        """Quasi-normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.

        """
        # handle negative reflections to minimise effect on mean I values.
        intensities.data().set_selected(intensities.data() < 0.0, 0.0)

        # set up binning objects
        if intensities.size() > 20000:
            n_refl_shells = 20
        elif intensities.size() > 15000:
            n_refl_shells = 15
        else:
            n_refl_shells = 10
        d_star_sq = intensities.d_star_sq().data()
        step = (flex.max(d_star_sq) - flex.min(d_star_sq) + 1e-8) / n_refl_shells
        intensities.setup_binner_d_star_sq_step(d_star_sq_step=step)

        normalisations = intensities.intensity_quasi_normalisations()
        return intensities.customized_copy(
            data=(intensities.data() / normalisations.data()),
            sigmas=(intensities.sigmas() / normalisations.data()),
        )

    @staticmethod
    def ml_aniso_normalisation(intensities):
        """Anisotropic maximum-likelihood normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.

        """
        return symmetry_base._ml_normalisation(intensities, aniso=True)

    @staticmethod
    def ml_iso_normalisation(intensities):
        """Isotropic maximum-likelihood normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.

        """
        return symmetry_base._ml_normalisation(intensities, aniso=False)

    @staticmethod
    def _ml_normalisation(intensities, aniso):
        # estimate number of residues per unit cell
        mr = matthews.matthews_rupp(intensities.crystal_symmetry())
        n_residues = mr.n_residues

        # estimate B-factor and scale factors for normalisation
        if aniso:
            normalisation = absolute_scaling.ml_aniso_absolute_scaling(
                intensities, n_residues=n_residues
            )
            u_star = normalisation.u_star
        else:
            normalisation = absolute_scaling.ml_iso_absolute_scaling(
                intensities, n_residues=n_residues
            )
            u_star = adptbx.b_as_u(
                adptbx.u_iso_as_u_star(intensities.unit_cell(), normalisation.b_wilson)
            )

        # record output in log file
        if aniso:
            b_cart = normalisation.b_cart
            logger.info("ML estimate of overall B_cart value:")
            logger.info(
                """\
  %5.2f, %5.2f, %5.2f
  %12.2f, %5.2f
  %19.2f"""
                % (b_cart[0], b_cart[3], b_cart[4], b_cart[1], b_cart[5], b_cart[2])
            )
        else:
            logger.info("ML estimate of overall B value:")
            logger.info("   %5.2f A**2" % normalisation.b_wilson)
        logger.info("ML estimate of  -log of scale factor:")
        logger.info("  %5.2f" % (normalisation.p_scale))

        s = StringIO()
        mr.show(out=s)
        normalisation.show(out=s)
        logger.debug(s.getvalue())

        # apply scales
        return intensities.customized_copy(
            data=scaling.ml_normalise_aniso(
                intensities.indices(),
                intensities.data(),
                normalisation.p_scale,
                intensities.unit_cell(),
                u_star,
            ),
            sigmas=scaling.ml_normalise_aniso(
                intensities.indices(),
                intensities.sigmas(),
                normalisation.p_scale,
                intensities.unit_cell(),
                u_star,
            ),
        )

    def _resolution_filter(self, d_min, min_i_mean_over_sigma_mean, min_cc_half):
        if d_min is libtbx.Auto and (
            min_i_mean_over_sigma_mean is not None or min_cc_half is not None
        ):
            from dials.util import Resolutionizer

            rparams = Resolutionizer.phil_defaults.extract().resolutionizer
            rparams.nbins = 20
            rparams.plot = False
            resolutionizer = Resolutionizer.resolutionizer(self.intensities, rparams)
            d_min_isigi = 0
            d_min_cc_half = 0
            if min_i_mean_over_sigma_mean is not None:
                try:
                    d_min_isigi = resolutionizer.resolution_i_mean_over_sigma_mean(
                        min_i_mean_over_sigma_mean
                    )
                except RuntimeError as e:
                    logger.info(
                        "I/sigI resolution filter failed with the following error:"
                    )
                    logger.error(e)
                else:
                    logger.info(
                        "Resolution estimate from <I>/<sigI> > %.1f : %.2f"
                        % (min_i_mean_over_sigma_mean, d_min_isigi)
                    )
            if min_cc_half is not None:
                try:
                    d_min_cc_half = resolutionizer.resolution_cc_half(min_cc_half)
                except RuntimeError as e:
                    logger.info(
                        "CChalf resolution filter failed with the following error:"
                    )
                    logger.error(e)
                else:
                    logger.info(
                        "Resolution estimate from CC1/2 > %.2f: %.2f"
                        % (min_cc_half, d_min_cc_half)
                    )
            valid_d_mins = list({d_min_cc_half, d_min_isigi}.difference({0}))
            if valid_d_mins:
                d_min = min(valid_d_mins)
                logger.info("High resolution limit set to: %.2f" % d_min)
        if d_min is not None:
            sel = self.intensities.resolution_filter_selection(d_min=d_min)
            self.intensities = self.intensities.select(sel).set_info(
                self.intensities.info()
            )
            self.dataset_ids = self.dataset_ids.select(sel)
            logger.info(
                "Selecting %i reflections with d > %.2f"
                % (self.intensities.size(), d_min)
            )
