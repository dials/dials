"""Methods for symmetry determination.

This module provides a base class for symmetry determination algorithms.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

from io import StringIO

import libtbx
from cctbx import adptbx, sgtbx, uctbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling, matthews
from scitbx.array_family import flex

from dials.util import resolution_analysis
from dials.util.normalisation import quasi_normalisation


class symmetry_base:
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
        best_monoclinic_beta=True,
    ):
        """Initialise a symmetry_base object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            symmetry analysis.
          normalisation (str): The normalisation method to use. Possible choices are
            'kernel', 'quasi', 'ml_iso' and 'ml_aniso'. Set to None to switch off
            normalisation altogether.
          lattice_symmetry_max_delta (float): The maximum value of delta for
            determining the lattice symmetry using the algorithm of Le Page (1982).
          d_min (float): Optional resolution cutoff to be applied to the input
            intensities. If set to :data:`libtbx.Auto` then d_min will be
            automatically determined according to the parameters
            ``min_i_mean_over_sigma_mean`` and ``min_cc_half``.
          min_i_mean_over_sigma_mean (float): minimum value of :math:`|I|/|sigma(I)|` for
            automatic determination of resolution cutoff.
          min_cc_half (float): minimum value of CC½ for automatic determination of
            resolution cutoff.
          relative_length_tolerance (float): Relative length tolerance in checking
            consistency of input unit cells against the median unit cell.
          absolute_angle_tolerance (float): Absolute angle tolerance in checking
            consistency of input unit cells against the median unit cell.
          best_monoclinic_beta (bool): If True, then for monoclinic centered cells, I2
            will be preferred over C2 if it gives a less oblique cell (i.e. smaller
            beta angle).
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
        self.dataset_ids = flex.int(self.intensities.size(), 0)
        for i, d in enumerate(self.input_intensities[1:]):
            self.intensities = self.intensities.concatenate(
                d, assert_is_similar_symmetry=False
            )
            self.dataset_ids.extend(flex.int(d.size(), i + 1))
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
            best_monoclinic_beta=best_monoclinic_beta,
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
        logger.info("Patterson group: %s", self.patterson_group.info())

        sel = self.patterson_group.epsilon(self.intensities.indices()) == 1
        self.intensities = self.intensities.select(sel)
        self.dataset_ids = self.dataset_ids.select(sel)

        # Correct SDs by "typical" SD factors
        self._correct_sigmas(sd_fac=2.0, sd_b=0.0, sd_add=0.03)
        self._normalise(normalisation)
        self._resolution_filter(d_min, min_i_mean_over_sigma_mean, min_cc_half)

    def _check_unit_cell_consistency(
        self, relative_length_tolerance, absolute_angle_tolerance
    ):
        for d in self.input_intensities:
            if (
                relative_length_tolerance is not None
                and absolute_angle_tolerance is not None
            ):
                if not d.unit_cell().is_similar_to(
                    self.median_unit_cell,
                    relative_length_tolerance,
                    absolute_angle_tolerance,
                ):
                    raise ValueError(
                        f"Incompatible unit cell: {d.unit_cell()}\n"
                        + f"      median unit cell: {self.median_unit_cell}"
                    )

    def _normalise(self, method):
        if method is None:
            return
        elif method == "kernel":
            normalise = self.kernel_normalisation
        elif method == "quasi":
            normalise = quasi_normalisation
        elif method == "ml_iso":
            normalise = self.ml_iso_normalisation
        elif method == "ml_aniso":
            normalise = self.ml_aniso_normalisation

        for i in range(int(flex.max(self.dataset_ids) + 1)):
            logger.info("\n" + "-" * 80 + "\n")
            logger.info("Normalising intensities for dataset %i\n", i + 1)
            intensities = self.intensities.select(self.dataset_ids == i)
            try:
                intensities = normalise(intensities)
            # Catch any of the several errors that can occur when there are too few
            # reflections for the selected normalisation routine.
            except (AttributeError, IndexError, RuntimeError):
                logger.warning(
                    "A problem occurred when trying to normalise the intensities by "
                    "the %s method.\n"
                    "There may be too few unique reflections. Resorting to using "
                    "un-normalised intensities instead.",
                    method,
                    exc_info=True,
                )
                return
            if i == 0:
                normalised_intensities = intensities
            else:
                normalised_intensities = normalised_intensities.concatenate(intensities)

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
  %19.2f""",
                b_cart[0],
                b_cart[3],
                b_cart[4],
                b_cart[1],
                b_cart[5],
                b_cart[2],
            )
        else:
            logger.info("ML estimate of overall B value:")
            logger.info("   %5.2f A**2", normalisation.b_wilson)
        logger.info("ML estimate of  -log of scale factor:")
        logger.info("  %5.2f", normalisation.p_scale)

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
        logger.info("\n" + "-" * 80 + "\n")
        logger.info("Estimation of resolution for Laue group analysis\n")
        if d_min is libtbx.Auto and (
            min_i_mean_over_sigma_mean is not None or min_cc_half is not None
        ):
            d_min = resolution_filter_from_array(
                self.intensities, min_i_mean_over_sigma_mean, min_cc_half
            )
            if d_min is not None:
                logger.info("High resolution limit set to: %.2f", d_min)
            else:
                logger.info("High resolution limit set to: None")
        if d_min is not None:
            sel = self.intensities.resolution_filter_selection(d_min=d_min)
            self.intensities = self.intensities.select(sel).set_info(
                self.intensities.info()
            )
            self.dataset_ids = self.dataset_ids.select(sel)
            logger.info(
                "Selecting %i reflections with d > %.2f", self.intensities.size(), d_min
            )


def resolution_filter_from_array(intensities, min_i_mean_over_sigma_mean, min_cc_half):
    """Run the resolution filter using miller array data format."""
    rparams = resolution_analysis.phil_defaults.extract().resolution
    resolutionizer = resolution_analysis.Resolutionizer(intensities, rparams)
    return _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half)


def resolution_filter_from_reflections_experiments(
    reflections, experiments, min_i_mean_over_sigma_mean, min_cc_half
):
    """Run the resolution filter using native dials data formats."""
    rparams = resolution_analysis.phil_defaults.extract().resolution
    resolutionizer = (
        resolution_analysis.Resolutionizer.from_reflections_and_experiments(
            reflections, experiments, rparams
        )
    )
    return _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half)


def _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half):
    """Use a resolutionizer to perform filtering for symmetry analysis."""
    d_min_isigi = 0
    d_min_cc_half = 0
    if min_i_mean_over_sigma_mean is not None:
        try:
            d_min_isigi = resolutionizer.resolution(
                resolution_analysis.metrics.I_MEAN_OVER_SIGMA_MEAN,
                limit=min_i_mean_over_sigma_mean,
            ).d_min
        except RuntimeError as e:
            logger.info("I/σ(I) resolution filter failed with the following error:")
            logger.error(e)
        else:
            if d_min_isigi:
                logger.info(
                    "Resolution estimate from <I>/<σ(I)> > %.1f : %.2f",
                    min_i_mean_over_sigma_mean,
                    d_min_isigi,
                )
    if min_cc_half is not None:
        try:
            d_min_cc_half = resolutionizer.resolution(
                resolution_analysis.metrics.CC_HALF, limit=min_cc_half
            ).d_min
        except RuntimeError as e:
            logger.info("CC½ resolution filter failed with the following error:")
            logger.error(e)
        else:
            if d_min_cc_half:
                logger.info(
                    "Resolution estimate from CC½ > %.2f: %.2f",
                    min_cc_half,
                    d_min_cc_half,
                )
    valid = [d for d in (d_min_cc_half, d_min_isigi) if d]
    if valid:
        return min(valid)
