from __future__ import division

import logging
from math import sqrt

from scitbx import matrix

from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.array_family import flex

logger = logging.getLogger("dials." + __name__)


class CrystalRefiner(object):
    def __init__(self, experiment, reflections, model):
        from dials.algorithms.refinement.parameterisation.crystal_parameters import (
            CrystalOrientationParameterisation,
            CrystalUnitCellParameterisation,
        )

        # Store the input
        self.experiment = experiment
        self.reflections = reflections
        self.model = model

        # Print RMSD
        Xobs, Yobs, _ = self.reflections["xyzobs.mm.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.mm"].parts()
        rmsd_x = sqrt(flex.sum((Xcal - Xobs) ** 2) / len(Xcal))
        rmsd_y = sqrt(flex.sum((Ycal - Yobs) ** 2) / len(Ycal))
        logger.info("Initial RMSD X, Y (mm): %f, %f" % (rmsd_x, rmsd_y))

        # Print RMSD
        Xobs, Yobs, _ = self.reflections["xyzobs.px.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.px"].parts()
        rmsd_x = sqrt(flex.sum((Xcal - Xobs) ** 2) / len(Xcal))
        rmsd_y = sqrt(flex.sum((Ycal - Yobs) ** 2) / len(Ycal))
        logger.info("Initial RMSD X, Y (px): %f, %f" % (rmsd_x, rmsd_y))

        # Get the crystal model and the parameterisation
        self.crystal = self.experiment.crystal
        self.cucp = CrystalUnitCellParameterisation(self.crystal)
        self.cop = CrystalOrientationParameterisation(self.crystal)

        # Get the current values and generate some offsets
        # values = flex.double(self.cucp.get_param_vals() + self.cop.get_param_vals())
        # offset = flex.double(
        #    [0.01 * v for v in self.cucp.get_param_vals()] + [0.1, 0.1, 0.1]
        # )

        # The optimization history
        self.history = []

        # Get the initial cell and initial score
        initial_cell = self.crystal.get_unit_cell()
        initial_orientation = self.crystal.get_U()
        # initial_score = self.target(values)

        # Perform the optimization
        # optimizer = SimpleSimplex(values, offset, self, 2000, tolerance=1e-3)
        # result = optimizer.get_solution()

        # Print some information
        fmt = "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f)"
        logger.info("Initial cell: %s" % initial_cell)
        logger.info("Final cell: %s" % self.crystal.get_unit_cell())
        logger.info("Initial orientation: %s" % (fmt % tuple(initial_orientation)))
        logger.info("Final orientation: %s" % (fmt % tuple(self.crystal.get_U())))

        # Print RMSD
        Xobs, Yobs, _ = self.reflections["xyzobs.mm.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.mm"].parts()
        rmsd_x = sqrt(flex.sum((Xcal - Xobs) ** 2) / len(Xcal))
        rmsd_y = sqrt(flex.sum((Ycal - Yobs) ** 2) / len(Ycal))
        logger.info("RMSD X, Y (mm): %f, %f" % (rmsd_x, rmsd_y))

        # Print RMSD
        Xobs, Yobs, _ = self.reflections["xyzobs.px.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.px"].parts()
        rmsd_x = sqrt(flex.sum((Xcal - Xobs) ** 2) / len(Xcal))
        rmsd_y = sqrt(flex.sum((Ycal - Yobs) ** 2) / len(Ycal))
        logger.info("RMSD X, Y (px): %f, %f" % (rmsd_x, rmsd_y))

    def target(self, vector):
        """
        The target function

        """
        from math import sqrt

        from dials.array_family import flex

        # Get the cell and orientation parameters
        cell_parms = self.cucp.get_param_vals()
        orientation_parms = self.cop.get_param_vals()
        assert len(vector) == len(cell_parms) + len(orientation_parms)

        # Update the cell and orientation parameters
        tst_cell = vector[: len(cell_parms)]
        tst_orientation = vector[
            len(cell_parms) : len(cell_parms) + len(orientation_parms)
        ]
        self.cucp.set_param_vals(tst_cell)
        self.cop.set_param_vals(tst_orientation)

        # Generate predicted positions
        s1_cal, s2_cal = self.generate_predictions(
            self.experiment, self.reflections, self.model
        )

        # Do the ray intersection
        self.reflections["s1"] = s1_cal
        self.reflections["s2"] = s2_cal
        self.reflections["xyzcal.mm"] = flex.vec3_double(
            [
                self.experiment.detector[0].get_ray_intersection(s1) + (0,)
                for s1 in s1_cal
            ]
        )
        self.reflections["xyzcal.px"] = flex.vec3_double(
            [
                self.experiment.detector[0].millimeter_to_pixel((mm[0], mm[1])) + (0,)
                for mm in self.reflections["xyzcal.mm"]
            ]
        )

        # Get predictions and observations
        Xobs, Yobs, _ = self.reflections["xyzobs.mm.value"].parts()
        Xcal, Ycal, _ = self.reflections["xyzcal.mm"].parts()

        # Compute the rmsd between observed and calculated
        score = flex.sum((Xobs - Xcal) ** 2 + (Yobs - Ycal) ** 2)

        # Append to the history
        self.history.append((tst_cell, tst_orientation, score))

        # Print some info
        logger.info(
            "Cell: %.3f %.3f %.3f %.3f %.3f %.3f; Phi: %.3f %.3f %.3f; RMSD: %.3f"
            % (
                tuple(self.crystal.get_unit_cell().parameters())
                + tuple(tst_orientation)
                + tuple((sqrt(score / len(Xobs)),))
            )
        )
        return score

    def generate_predictions(self, experiment, reflections, model):

        # The crystal A and beam s0
        A = matrix.sqr(experiment.crystal.get_A())
        s0 = matrix.col(experiment.beam.get_s0())
        s0_length = s0.length()

        # Compute all the vectors
        s1_cal = flex.vec3_double()
        s2_cal = flex.vec3_double()
        for i in range(len(reflections)):

            # Compute the reciprocal lattice vector
            h = matrix.col(reflections[i]["miller_index"])
            r = A * h
            s2 = s0 + r

            # Rotate the covariance matrix
            R = compute_change_of_basis_operation(s0, s2)
            S = R * model.sigma() * R.transpose()
            mu = R * s2
            assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7

            # Partition the mean vector
            mu1 = matrix.col((mu[0], mu[1]))
            mu2 = mu[2]

            # Partition the covariance matrix
            # S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
            S12 = matrix.col((S[2], S[5]))
            # S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            # Compute the conditional mean
            mubar = mu1 + S12 * (1 / S22) * (s0_length - mu2)

            # Compute the vector and rotate
            v = matrix.col((mubar[0], mubar[1], s0_length)).normalize() * s0_length
            s1 = R.transpose() * v

            # Append the 2 vectors
            s1_cal.append(s1)
            s2_cal.append(s2)

        # Return the predicted vectors
        return s1_cal, s2_cal
