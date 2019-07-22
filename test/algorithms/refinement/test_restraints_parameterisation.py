from __future__ import absolute_import, division, print_function

# Tests for RestraintsParameterisation and associated classes used in refinement

import os
import random

from dials.algorithms.refinement import RefinerFactory
from dials.array_family import flex
from dials.algorithms.refinement.restraints import RestraintsParameterisation
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.phil import parse
import pytest


def test_single_crystal_restraints_gradients():
    """Simple test with a single triclinic crystal restrained to a target unit cell"""

    from dials.test.algorithms.refinement.setup_geometry import Extract
    from dxtbx.model.experiment_list import ExperimentList, Experiment

    #### Import model parameterisations

    from dials.algorithms.refinement.parameterisation.prediction_parameters import (
        XYPhiPredictionParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.detector_parameters import (
        DetectorParameterisationSinglePanel,
    )
    from dials.algorithms.refinement.parameterisation.beam_parameters import (
        BeamParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalOrientationParameterisation,
        CrystalUnitCellParameterisation,
    )

    overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""

    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      """,
        process_includes=True,
    )

    models = Extract(master_phil, overrides)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 72 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 720),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(720)),
        deg=True,
    )

    # Create parameterisations of these models
    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    # Create an ExperimentList
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )

    # Build a prediction parameterisation
    pred_param = XYPhiPredictionParameterisation(
        experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
    )

    # Build a restraints parameterisation
    rp = RestraintsParameterisation(
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
    )

    # make a unit cell target
    sigma = 1.0
    uc = mycrystal.get_unit_cell().parameters()
    target_uc = [random.gauss(e, sigma) for e in uc]

    rp.add_restraints_to_target_xl_unit_cell(
        experiment_id=0, values=target_uc, sigma=[sigma] * 6
    )

    # get analytical values and gradients
    vals, grads, weights = rp.get_residuals_gradients_and_weights()
    assert len(vals) == rp.num_residuals

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grad = []

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        rev_state, foo, bar = rp.get_residuals_gradients_and_weights()
        rev_state = flex.double(rev_state)

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        fwd_state, foo, bar = rp.get_residuals_gradients_and_weights()
        fwd_state = flex.double(fwd_state)

        p_vals[i] = val

        fd = (fwd_state - rev_state) / delta
        fd_grad.append(fd)

    # for comparison, fd_grad is a list of flex.doubles, each of which corresponds
    # to a column of the sparse matrix grads.
    for i, fd in enumerate(fd_grad):
        # extract dense column from the sparse matrix
        an = grads.col(i).as_dense_vector()

        assert an == pytest.approx(fd, abs=1e-5)


def test_two_triclinic_crystals():
    """Simple test with two triclinic crystals restrained to a target unit cell"""

    from dials.test.algorithms.refinement.setup_geometry import Extract
    from dxtbx.model.experiment_list import ExperimentList, Experiment

    #### Import model parameterisations

    from dials.algorithms.refinement.parameterisation.prediction_parameters import (
        XYPhiPredictionParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.detector_parameters import (
        DetectorParameterisationSinglePanel,
    )
    from dials.algorithms.refinement.parameterisation.beam_parameters import (
        BeamParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalOrientationParameterisation,
        CrystalUnitCellParameterisation,
    )

    overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""

    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      """,
        process_includes=True,
    )

    models = Extract(master_phil, overrides)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    # duplicate the crystal
    from copy import deepcopy

    mycrystal2 = deepcopy(mycrystal)
    mybeam = models.beam

    # Build a mock scan for a 72 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 720),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(720)),
        deg=True,
    )

    # Create parameterisations of these models
    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)
    xluc_param2 = CrystalUnitCellParameterisation(mycrystal2, experiment_ids=[1])

    # Create an ExperimentList with the crystal duplicated
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal2,
            imageset=None,
        )
    )

    # Build a prediction parameterisation
    pred_param = XYPhiPredictionParameterisation(
        experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param, xluc_param2],
    )

    # Build a restraints parameterisation
    rp = RestraintsParameterisation(
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param, xluc_param2],
    )

    # make a unit cell target
    sigma = 1.0
    uc = mycrystal.get_unit_cell().parameters()
    target_uc = [random.gauss(e, sigma) for e in uc]

    rp.add_restraints_to_target_xl_unit_cell(
        experiment_id=0, values=target_uc, sigma=[sigma] * 6
    )
    rp.add_restraints_to_target_xl_unit_cell(
        experiment_id=1, values=target_uc, sigma=[sigma] * 6
    )

    # get analytical values and gradients
    vals, grads, weights = rp.get_residuals_gradients_and_weights()
    assert len(vals) == rp.num_residuals

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grad = []

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        rev_state, foo, bar = rp.get_residuals_gradients_and_weights()
        rev_state = flex.double(rev_state)

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        fwd_state, foo, bar = rp.get_residuals_gradients_and_weights()
        fwd_state = flex.double(fwd_state)

        p_vals[i] = val

        fd = (fwd_state - rev_state) / delta
        fd_grad.append(fd)

    # for comparison, fd_grad is a list of flex.doubles, each of which corresponds
    # to a column of the sparse matrix grads.
    for i, fd in enumerate(fd_grad):
        # extract dense column from the sparse matrix
        an = grads.col(i).as_dense_vector()
        assert an == pytest.approx(fd, abs=1e-5)


def test_10_crystals_with_stills_parameterisation(dials_regression):
    """Test with multiple crystals, and a stills refiner"""

    # The phil scope
    from dials.algorithms.refinement.refiner import phil_scope

    user_phil = parse(
        """
  refinement
  {
    parameterisation
    {
      crystal
      {
        unit_cell
        {
          restraints
          {
            tie_to_target
            {
              values=95,95,132,90,90,120
              sigmas=1,1,1,0,0,0
              id=0,1,2,3,4,5,6,7,8,9
            }
          }
        }
      }
    }
  }
  """
    )

    working_phil = phil_scope.fetch(source=user_phil)
    working_params = working_phil.extract()

    # use the multi stills test data
    data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
    experiments_path = os.path.join(data_dir, "combined_experiments.json")
    pickle_path = os.path.join(data_dir, "combined_reflections.pickle")

    experiments = ExperimentListFactory.from_json_file(
        experiments_path, check_format=False
    )
    reflections = flex.reflection_table.from_pickle(pickle_path)

    refiner = RefinerFactory.from_parameters_data_experiments(
        working_params, reflections, experiments
    )

    # hack to extract the objects needed from the Refiner
    rp = refiner._target._restraints_parameterisation
    pred_param = refiner._pred_param

    # get analytical values and gradients
    vals, grads, weights = rp.get_residuals_gradients_and_weights()
    assert len(vals) == rp.num_residuals

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grad = []

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        rev_state, foo, bar = rp.get_residuals_gradients_and_weights()
        rev_state = flex.double(rev_state)

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        fwd_state, foo, bar = rp.get_residuals_gradients_and_weights()
        fwd_state = flex.double(fwd_state)

        p_vals[i] = val

        fd = (fwd_state - rev_state) / delta
        fd_grad.append(fd)

    # for comparison, fd_grad is a list of flex.doubles, each of which corresponds
    # to a column of the sparse matrix grads.
    pnames = pred_param.get_param_names()
    for i, (pname, fd) in enumerate(zip(pnames, fd_grad)):
        # extract dense column from the sparse matrix
        an = grads.col(i).as_dense_vector()

        # print pname
        # print list(an.round(6))
        # print list(fd.round(6))
        # print
        assert an == pytest.approx(fd, abs=1e-5)


def test_group_restraint_with_multiple_crystals_and_a_stills_refiner(dials_regression):

    # The phil scope
    from dials.algorithms.refinement.refiner import phil_scope

    user_phil = parse(
        """
  refinement
  {
    parameterisation
    {
      crystal
      {
        unit_cell
        {
          restraints
          {
            tie_to_group
            {
              sigmas=1,0,2,0,0,0
            }
          }
        }
      }
    }
  }
  """
    )

    working_phil = phil_scope.fetch(source=user_phil)
    working_params = working_phil.extract()

    # use the multi stills test data
    data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
    experiments_path = os.path.join(data_dir, "combined_experiments.json")
    pickle_path = os.path.join(data_dir, "combined_reflections.pickle")

    experiments = ExperimentListFactory.from_json_file(
        experiments_path, check_format=False
    )
    reflections = flex.reflection_table.from_pickle(pickle_path)

    refiner = RefinerFactory.from_parameters_data_experiments(
        working_params, reflections, experiments
    )

    # hack to extract the objects needed from the Refiner
    rp = refiner._target._restraints_parameterisation
    pred_param = refiner._pred_param

    # get analytical values and gradients
    vals, grads, weights = rp.get_residuals_gradients_and_weights()
    assert len(vals) == rp.num_residuals

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grad = []

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        rev_state, foo, bar = rp.get_residuals_gradients_and_weights()
        rev_state = flex.double(rev_state)

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        fwd_state, foo, bar = rp.get_residuals_gradients_and_weights()
        fwd_state = flex.double(fwd_state)

        p_vals[i] = val

        fd = (fwd_state - rev_state) / delta
        fd_grad.append(fd)

    # for comparison, fd_grad is a list of flex.doubles, each of which corresponds
    # to the gradients of the residuals wrt to a single parameter.
    pnames = pred_param.get_param_names()
    for i, (pname, fd) in enumerate(zip(pnames, fd_grad)):
        # extract dense column from the sparse matrix
        an = grads.col(i).as_dense_vector()

        # print pname
        # print list(an.round(6))
        # print list(fd.round(6))
        # print
        assert an == pytest.approx(fd, abs=1e-5)
