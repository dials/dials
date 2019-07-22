from __future__ import absolute_import, division, print_function

import math
from collections import namedtuple

from scitbx import matrix
from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    XYPhiPredictionParameterisation,
    SparseGradientVectorMixin,
)


class StateDerivativeCache(object):
    """Keep derivatives of the model states in a memory-efficient format
    by storing each derivative once alongside the indices of reflections affected
    by that derivative"""

    def __init__(self, parameterisations=None):

        if parameterisations is None:
            parameterisations = []
        self._cache = dict.fromkeys(parameterisations)

        self._Pair = namedtuple("Pair", ["derivative", "iselection"])

        # set up lists with the right number of elements
        self.clear()

        self._nref = 0

    def build_gradients(self, parameterisation, isel=None, imatch=None):
        """Return an object mimicking a list of flex arrays containing state
        gradients wrt each parameter of the parameterisation. In fact this is a
        generator so that iterating over the elements of the list will return
        control here so that the gradient array for a single parameter can be
        reconstructed on the fly"""

        # Get the data from the cache
        entry = self._cache[parameterisation]

        # Figure out the right flex array type from entries in the cache
        shape = None
        for e in entry:
            if e:
                shape = e[0].derivative.n
                break
        if shape is None:
            raise TypeError("No model state derivatives found")
        if shape == (3, 1):
            arr_type = flex.vec3_double
            null = (0, 0, 0)
        elif shape == (3, 3):
            arr_type = flex.mat3_double
            null = (0, 0, 0, 0, 0, 0, 0, 0, 0)
        else:
            raise TypeError("Unrecognised model state derivative type")

        # Loop over the data for each parameter
        for p_data in entry:

            # Build an empty array of the same length as the original reflection
            # list used when the cache was filled
            ds_dp = arr_type(self._nref, null)

            # Reconstitute full array from the cache
            for pair in p_data:
                ds_dp.set_selected(pair.iselection, pair.derivative)

            # First select only elements relevant to the current gradient calculation
            # block (i.e. if nproc > 1 or gradient_calculation_blocksize was set)
            if imatch is not None:
                ds_dp = ds_dp.select(imatch)

            # Now select only those reflections from the full list that are affected
            # by this parameterisation
            if isel is not None:
                ds_dp = ds_dp.select(isel)

            yield ds_dp

    def clear(self):
        """Clear all cached values"""

        for p in self._cache:
            self._cache[p] = [[] for i in range(p.num_free())]

    def append(self, parameterisation, iparam, derivative, iselection):
        """For a particular parameterisation and parameter number of the free
        parameters of that parameterisation, append a state derivative and the
        iselection of reflections it affects to the cache"""

        l1 = self._cache[parameterisation]
        l2 = l1[iparam]
        l2.append(self._Pair(derivative, iselection))

    @property
    def nref(self):
        """Get the length of the reflection list to which indices in the iselections
        refer"""

        return self._nref

    @nref.setter
    def nref(self, value):
        """Set the length of the reflection list to which indices in the iselections
        refer"""

        self._nref = value


class ScanVaryingPredictionParameterisation(XYPhiPredictionParameterisation):
    """An extension of the rotation scans version of the
    PredictionParameterisation class that supports model parameterisations that
    vary smoothly with the observed image number"""

    def __init__(
        self,
        experiments,
        detector_parameterisations=None,
        beam_parameterisations=None,
        xl_orientation_parameterisations=None,
        xl_unit_cell_parameterisations=None,
        goniometer_parameterisations=None,
    ):

        if detector_parameterisations is None:
            detector_parameterisations = []
        if beam_parameterisations is None:
            beam_parameterisations = []
        if xl_orientation_parameterisations is None:
            xl_orientation_parameterisations = []
        if xl_unit_cell_parameterisations is None:
            xl_unit_cell_parameterisations = []
        if goniometer_parameterisations is None:
            goniometer_parameterisations = []

        # determine once here which types of parameterisations are scan-varying
        self._varying_detectors = any(
            hasattr(p, "num_sets") for p in detector_parameterisations
        )
        self._varying_beams = any(
            hasattr(p, "num_sets") for p in beam_parameterisations
        )
        self._varying_xl_orientations = any(
            hasattr(p, "num_sets") for p in xl_orientation_parameterisations
        )
        self._varying_xl_unit_cells = any(
            hasattr(p, "num_sets") for p in xl_unit_cell_parameterisations
        )
        self._varying_goniometers = any(
            hasattr(p, "num_sets") for p in goniometer_parameterisations
        )

        to_cache = []
        if self._varying_detectors:
            to_cache.extend(detector_parameterisations)
        if self._varying_beams:
            to_cache.extend(beam_parameterisations)
        if self._varying_xl_orientations:
            to_cache.extend(xl_orientation_parameterisations)
        if self._varying_xl_unit_cells:
            to_cache.extend(xl_unit_cell_parameterisations)
        if self._varying_goniometers:
            to_cache.extend(goniometer_parameterisations)
        self._derivative_cache = StateDerivativeCache(to_cache)

        # set up base class
        super(ScanVaryingPredictionParameterisation, self).__init__(
            experiments,
            detector_parameterisations=detector_parameterisations,
            beam_parameterisations=beam_parameterisations,
            xl_orientation_parameterisations=xl_orientation_parameterisations,
            xl_unit_cell_parameterisations=xl_unit_cell_parameterisations,
            goniometer_parameterisations=goniometer_parameterisations,
        )

        # Avoid calculation in calculate_model_state_uncertainties unless this
        # is set to True
        self.set_scan_varying_errors = False

    def _get_xl_orientation_parameterisation(self, experiment_id):
        """Return the crystal orientation parameterisation for the requested
        experiment number (or None if the crystal orientation in that experiment
        is not parameterised)"""

        param_set = self._exp_to_param[experiment_id]
        xl_op = None
        if param_set.xl_ori_param is not None:
            xl_op = self._xl_orientation_parameterisations[param_set.xl_ori_param]

        return xl_op

    def _get_xl_unit_cell_parameterisation(self, experiment_id):
        """Return the crystal unit cell parameterisation for the requested
        experiment number (or None if the crystal unit cell in that experiment
        is not parameterised)"""

        param_set = self._exp_to_param[experiment_id]
        xl_ucp = None
        if param_set.xl_uc_param is not None:
            xl_ucp = self._xl_unit_cell_parameterisations[param_set.xl_uc_param]

        return xl_ucp

    def _get_beam_parameterisation(self, experiment_id):
        """Return the beam parameterisation for the requested experiment number
        (or None if the beam in that experiment is not parameterised)"""

        param_set = self._exp_to_param[experiment_id]
        bp = None
        if param_set.beam_param is not None:
            bp = self._beam_parameterisations[param_set.beam_param]

        return bp

    def _get_detector_parameterisation(self, experiment_id):
        """Return the detector parameterisation for the requested experiment number
        (or None if the detector in that experiment is not parameterised)"""

        param_set = self._exp_to_param[experiment_id]
        dp = None
        if param_set.det_param is not None:
            dp = self._detector_parameterisations[param_set.det_param]

        return dp

    def _get_goniometer_parameterisation(self, experiment_id):
        """Return the goniometer parameterisation for the requested experiment number
        (or None if the goniometer in that experiment is not parameterised)"""

        param_set = self._exp_to_param[experiment_id]
        gp = None
        if param_set.gonio_param is not None:
            gp = self._goniometer_parameterisations[param_set.gonio_param]

        return gp

    def _get_state_from_parameterisation(
        self, parameterisation, frame, multi_state_elt=None
    ):
        """Get the model state from the parameterisation at the specified frame,
        taking care of whether it is a scan-varying parameterisation or not"""

        if parameterisation is None:
            return None
        if (
            hasattr(parameterisation, "num_sets")
            and not self._current_frame.get(parameterisation) == frame
        ):
            parameterisation.compose(frame)
            self._current_frame[parameterisation] = frame
        if multi_state_elt is None:
            state = parameterisation.get_state()
        else:
            state = parameterisation.get_state(multi_state_elt=multi_state_elt)
        return state

    def _prepare_for_compose(self, reflections, skip_derivatives=False):
        """Add columns to the reflection table to hold the varying state matrices
        or vectors for the experimental models, if required. Also prepare the cache
        for the derivatives of states that are scan-varying"""

        nref = len(reflections)
        # set columns if needed
        if "u_matrix" not in reflections:
            reflections["u_matrix"] = flex.mat3_double(nref)
        if "b_matrix" not in reflections:
            reflections["b_matrix"] = flex.mat3_double(nref)
        if "s0_vector" not in reflections:
            reflections["s0_vector"] = flex.vec3_double(nref)
        if "d_matrix" not in reflections:
            reflections["d_matrix"] = flex.mat3_double(nref)
        if "D_matrix" not in reflections:
            reflections["D_matrix"] = flex.mat3_double(nref)
        if "S_matrix" not in reflections:
            reflections["S_matrix"] = flex.mat3_double(nref)

        # Clear the state derivative cache and set the number of reflections needed
        # to reconstruct the derivative arrays later
        self._derivative_cache.clear()
        self._derivative_cache.nref = nref

    def compose(self, reflections, skip_derivatives=False):
        """Compose scan-varying crystal parameterisations at the specified image
        number, for the specified experiment, for each image. Put the varying
        matrices in the reflection table, and cache the derivatives."""

        self._prepare_for_compose(reflections, skip_derivatives)

        for iexp, exp in enumerate(self._experiments):

            # select the reflections of interest
            sel = reflections["id"] == iexp
            isel = sel.iselection()

            blocks = reflections["block"].select(isel)

            # identify which parameterisations to use for this experiment
            xl_op = self._get_xl_orientation_parameterisation(iexp)
            xl_ucp = self._get_xl_unit_cell_parameterisation(iexp)
            bp = self._get_beam_parameterisation(iexp)
            dp = self._get_detector_parameterisation(iexp)
            gp = self._get_goniometer_parameterisation(iexp)

            # reset current frame cache for scan-varying parameterisations
            self._current_frame = {}

            # get state and derivatives for each block
            for block in range(flex.min(blocks), flex.max(blocks) + 1):

                # determine the subset of reflections this affects
                subsel = isel.select(blocks == block)
                if len(subsel) == 0:
                    continue

                # get the panels hit by these reflections
                panels = reflections["panel"].select(subsel)

                # get the integer frame number nearest the centre of that block
                frames = reflections["block_centre"].select(subsel)

                # can only be false if original block assignment has gone wrong
                assert frames.all_eq(
                    frames[0]
                ), "Failing: a block contains reflections that shouldn't be there"
                frame = int(math.floor(frames[0]))

                # model states at current frame
                U = self._get_state_from_parameterisation(xl_op, frame)
                if U is None:
                    U = matrix.sqr(exp.crystal.get_U())

                B = self._get_state_from_parameterisation(xl_ucp, frame)
                if B is None:
                    B = matrix.sqr(exp.crystal.get_B())

                s0 = self._get_state_from_parameterisation(bp, frame)
                if s0 is None:
                    s0 = matrix.col(exp.beam.get_s0())

                S = self._get_state_from_parameterisation(gp, frame)
                if S is None:
                    S = matrix.sqr(exp.goniometer.get_setting_rotation())

                # set states for crystal, beam and goniometer
                reflections["u_matrix"].set_selected(subsel, U.elems)
                reflections["b_matrix"].set_selected(subsel, B.elems)
                reflections["s0_vector"].set_selected(subsel, s0.elems)
                reflections["S_matrix"].set_selected(subsel, S.elems)

                # set states and derivatives for this detector
                if dp is not None:  # detector is parameterised
                    if dp.is_multi_state():  # parameterised detector is multi panel

                        # loop through the panels in this detector
                        for panel_id, _ in enumerate(exp.detector):

                            # get the right subset of array indices to set for this panel
                            subsel2 = subsel.select(panels == panel_id)
                            if len(subsel2) == 0:
                                # if no reflections intersect this panel, skip calculation
                                continue

                            dmat = self._get_state_from_parameterisation(
                                dp, frame, multi_state_elt=panel_id
                            )
                            if dmat is None:
                                dmat = exp.detector[panel_id].get_d_matrix()
                            Dmat = exp.detector[panel_id].get_D_matrix()
                            reflections["d_matrix"].set_selected(subsel2, dmat)
                            reflections["D_matrix"].set_selected(subsel2, Dmat)

                            if self._varying_detectors and not skip_derivatives:
                                for j, dd in enumerate(
                                    dp.get_ds_dp(
                                        multi_state_elt=panel_id, use_none_as_null=True
                                    )
                                ):
                                    if dd is None:
                                        continue
                                    self._derivative_cache.append(dp, j, dd, subsel)

                    else:  # parameterised detector is single panel
                        dmat = self._get_state_from_parameterisation(dp, frame)
                        if dmat is None:
                            dmat = exp.detector[0].get_d_matrix()
                        Dmat = exp.detector[0].get_D_matrix()
                        reflections["d_matrix"].set_selected(subsel, dmat)
                        reflections["D_matrix"].set_selected(subsel, Dmat)

                        if self._varying_detectors and not skip_derivatives:
                            for j, dd in enumerate(dp.get_ds_dp(use_none_as_null=True)):
                                if dd is None:
                                    continue
                                self._derivative_cache.append(dp, j, dd, subsel)

                else:  # set states for unparameterised detector (dp is None)
                    # loop through the panels in this detector
                    for panel_id, _ in enumerate(exp.detector):

                        # get the right subset of array indices to set for this panel
                        subsel2 = subsel.select(panels == panel_id)
                        if len(subsel2) == 0:
                            # if no reflections intersect this panel, skip to the next
                            continue

                        dmat = exp.detector[panel_id].get_d_matrix()
                        Dmat = exp.detector[panel_id].get_D_matrix()
                        reflections["d_matrix"].set_selected(subsel2, dmat)
                        reflections["D_matrix"].set_selected(subsel2, Dmat)

                # set derivatives of the states for crystal, beam and goniometer
                if not skip_derivatives:
                    if xl_op is not None and self._varying_xl_orientations:
                        for j, dU in enumerate(xl_op.get_ds_dp(use_none_as_null=True)):
                            if dU is None:
                                continue
                            self._derivative_cache.append(xl_op, j, dU, subsel)
                    if xl_ucp is not None and self._varying_xl_unit_cells:
                        for j, dB in enumerate(xl_ucp.get_ds_dp(use_none_as_null=True)):
                            if dB is None:
                                continue
                            self._derivative_cache.append(xl_ucp, j, dB, subsel)
                    if bp is not None and self._varying_beams:
                        for j, ds0 in enumerate(bp.get_ds_dp(use_none_as_null=True)):
                            if ds0 is None:
                                continue
                            self._derivative_cache.append(bp, j, ds0, subsel)
                    if gp is not None and self._varying_goniometers:
                        for j, dS in enumerate(gp.get_ds_dp(use_none_as_null=True)):
                            if dS is None:
                                continue
                            self._derivative_cache.append(gp, j, dS, subsel)

        # set the UB matrices for prediction
        reflections["ub_matrix"] = reflections["u_matrix"] * reflections["b_matrix"]

    # called by refiner.run for setting the crystal scan points
    def get_varying_UB(self, obs_image_numbers, experiment_id):
        """Extract the setting matrix from the contained scan-dependent crystal
        parameterisations at specified image number."""

        if not (self._varying_xl_unit_cells or self._varying_xl_orientations):
            return None

        # identify which crystal parameterisations to use for this experiment
        xl_op = self._get_xl_orientation_parameterisation(experiment_id)
        xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

        UB_list = []
        for i in obs_image_numbers:
            U = self._get_state_from_parameterisation(xl_op, i)
            B = self._get_state_from_parameterisation(xl_ucp, i)
            UB_list.append(U * B)

        return UB_list

    # called by refiner.run for setting the beam scan points
    def get_varying_s0(self, obs_image_numbers, experiment_id):
        """Extract the s0 vector from the contained scan-dependent beam
        parameterisation at specified image number."""

        if not self._varying_beams:
            return None

        # identify which beam parameterisation to use for this experiment
        bp = self._get_beam_parameterisation(experiment_id)

        s0_list = []
        for i in obs_image_numbers:
            s0 = self._get_state_from_parameterisation(bp, i)
            s0_list.append(s0)

        return s0_list

    # called by refiner.run for setting the goniometer scan points
    def get_varying_setting_rotation(self, obs_image_numbers, experiment_id):
        """Extract the S matrix from the contained scan-dependent goniometer
        parameterisation at specified image number."""

        if not self._varying_goniometers:
            return None

        # identify which goniometer parameterisation to use for this experiment
        gp = self._get_goniometer_parameterisation(experiment_id)

        S_list = []
        for i in obs_image_numbers:
            S = self._get_state_from_parameterisation(gp, i)
            S_list.append(S)

        return S_list

    # overloaded for the scan-varying case
    def _get_model_data_for_experiment(self, experiment, reflections):
        """helper function to return model data s0, U, B, D and S for a particular
        experiment. In this scan-varying overload this is trivial because these
        values are already set as arrays in the reflection table"""

        return {
            "s0": reflections["s0_vector"],
            "U": reflections["u_matrix"],
            "B": reflections["b_matrix"],
            "D": reflections["D_matrix"],
            "S": reflections["S_matrix"],
        }

    def _beam_derivatives(self, isel, parameterisation, reflections):
        """Determine whether ds0_dp was precalculated then call the base class
        method"""

        if self._varying_beams:
            if "imatch" in reflections:
                imatch = reflections["imatch"]
            else:
                imatch = None
            ds0_dxluc_p = self._derivative_cache.build_gradients(
                parameterisation=parameterisation, isel=isel, imatch=imatch
            )
        else:
            ds0_dxluc_p = None

        return super(ScanVaryingPredictionParameterisation, self)._beam_derivatives(
            isel, parameterisation, ds0_dxluc_p
        )

    def _xl_orientation_derivatives(self, isel, parameterisation, reflections):
        """Determine whether dU_dp was precalculated then call the base class
        method"""

        if self._varying_xl_orientations:
            if "imatch" in reflections:
                imatch = reflections["imatch"]
            else:
                imatch = None
            dU_dxlo_p = self._derivative_cache.build_gradients(
                parameterisation=parameterisation, isel=isel, imatch=imatch
            )
        else:
            dU_dxlo_p = None

        return super(
            ScanVaryingPredictionParameterisation, self
        )._xl_orientation_derivatives(isel, parameterisation, dU_dxlo_p)

    def _xl_unit_cell_derivatives(self, isel, parameterisation, reflections):
        """Determine whether dB_dp was precalculated then call the base class
        method"""

        if self._varying_xl_unit_cells:
            if "imatch" in reflections:
                imatch = reflections["imatch"]
            else:
                imatch = None
            dB_dxluc_p = self._derivative_cache.build_gradients(
                parameterisation=parameterisation, isel=isel, imatch=imatch
            )
        else:
            dB_dxluc_p = None

        return super(
            ScanVaryingPredictionParameterisation, self
        )._xl_unit_cell_derivatives(isel, parameterisation, dB_dxluc_p)

    def _detector_derivatives(self, isel, panel_id, parameterisation, reflections):
        """Determine whether dd_dp was precalculated then call the base class
        method"""

        if self._varying_detectors:
            if "imatch" in reflections:
                imatch = reflections["imatch"]
            else:
                imatch = None
            dd_ddet_p = self._derivative_cache.build_gradients(
                parameterisation=parameterisation, isel=isel, imatch=imatch
            )
        else:
            dd_ddet_p = None

        return super(ScanVaryingPredictionParameterisation, self)._detector_derivatives(
            isel, panel_id, parameterisation, dd_ddet_p
        )

    def _goniometer_derivatives(self, isel, parameterisation, reflections):
        """Determine whether dS_dp was precalculated then call the base class
        method"""

        if self._varying_goniometers:
            if "imatch" in reflections:
                imatch = reflections["imatch"]
            else:
                imatch = None
            dS_dgon_p = self._derivative_cache.build_gradients(
                parameterisation=parameterisation, isel=isel, imatch=imatch
            )
        else:
            dS_dgon_p = None

        return super(
            ScanVaryingPredictionParameterisation, self
        )._goniometer_derivatives(isel, parameterisation, dS_dgon_p)

    def calculate_model_state_uncertainties(
        self, var_cov=None, obs_image_number=None, experiment_id=None
    ):
        """Take a variance-covariance matrix of all free parameters (probably
        calculated by a minimisation engine). For each parameterisation in the
        global model, extract the subset of this matrix for the associated block
        of parameters. Pass this on to the relevant model parameterisation to
        calculate its own uncertainty of state.

        This scan-varying version should first be called with var_cov set but
        obs_image_number=None and experiment_id=None. This calls the scan-static
        version to do the calculation for the scan-static parameterisations and
        also caches the subsets of var_cov relevant for the scan-varying
        parameterisations. Subsequent calls should provide obs_image_number and
        experiment_id to calculate for a particular crystal at a particular
        scan-point"""

        # First call, only a variance-covariance matrix is supplied
        if var_cov is not None:
            assert [obs_image_number, experiment_id].count(None) == 2
            super(
                ScanVaryingPredictionParameterisation, self
            ).calculate_model_state_uncertainties(var_cov)
            return

        # Later calls, only an experiment and image number are supplied for
        # identify the crystal parameterisations for this experiment
        xl_op = self._get_xl_orientation_parameterisation(experiment_id)
        xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

        result = {}

        # compose at the requested image number and calculate using the cached
        # varcov matrices. Take the first elt of the list because the crystal
        # parameterisations are not multi-state
        if xl_op is not None:
            try:
                xl_op.compose(obs_image_number)
                result["U_cov"] = xl_op.calculate_state_uncertainties(var_cov=None)[0]
            except TypeError:
                pass

        if xl_ucp is not None:
            try:
                xl_ucp.compose(obs_image_number)
                result["B_cov"] = xl_ucp.calculate_state_uncertainties(var_cov=None)[0]
            except TypeError:
                pass

        return result

    def set_model_state_uncertainties(self, u_cov_list, b_cov_list, experiment_id=None):
        """Identify the crystal parameterisations and set the list of covariance
        matrices, if available. They will only be available if the parameterisation
        is a scan-varying type, otherwise they are None"""

        xl_op = self._get_xl_orientation_parameterisation(experiment_id)
        xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

        if u_cov_list:
            try:
                xl_op.set_state_uncertainties(u_cov_list)
            except AttributeError:
                pass

        if b_cov_list:
            try:
                xl_ucp.set_state_uncertainties(b_cov_list)
            except AttributeError:
                pass


class ScanVaryingPredictionParameterisationSparse(
    SparseGradientVectorMixin, ScanVaryingPredictionParameterisation
):
    pass
