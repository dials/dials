from __future__ import annotations

import bz2
import collections
import math
import pickle
import shutil

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

from dials.array_family import flex


@pytest.fixture(scope="module")
def data(dials_data, tmp_dir):  # read experiments and reflections
    directory = dials_data("integration_test_data", pathlib=True)
    experiments_filename = str(directory / "thaumatin_i04-integrated.expt")
    reflections_filename = str(directory / "thaumatin_i04-shoeboxes_0_0.refl.bz2")
    reference_filename = str(directory / "thaumatin_i04-reference_profiles.pickle.bz2")

    for f in [reflections_filename, reference_filename]:
        with bz2.BZ2File(f) as compr:
            with open(tmp_dir / f.name[:-4], "wb") as decompr:
                shutil.copyfileobj(compr, decompr)

    experiments = ExperimentListFactory.from_json_file(
        experiments_filename, check_format=False
    )
    reflections = flex.reflection_table.from_file(
        tmp_dir / reflections_filename.name[:-4]
    )
    with open(tmp_dir / reference_filename.name[:-4], "rb") as fh:
        reference = pickle.load(fh, encoding="bytes")

    Data = collections.namedtuple("Data", ["experiments", "reflections", "reference"])
    return Data(experiments=experiments, reflections=reflections, reference=reference)


def test_gaussianrs_mask_calculator(data):
    from dials.algorithms.integration.parallel_integrator import MaskCalculatorFactory

    algorithm = MaskCalculatorFactory.create(data.experiments)
    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    for r in reflections:
        algorithm(r, False)


def test_simple_background_calculator(data):
    from dials.algorithms.background.simple.algorithm import (
        SimpleBackgroundCalculatorFactory,
    )

    algorithm = SimpleBackgroundCalculatorFactory.create(data.experiments)
    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r)
        except Exception:
            errors += 1
    assert len(reflections) == 15193
    assert errors == 333


def test_glm_background_calculator(data):
    from dials.algorithms.background.glm.algorithm import GLMBackgroundCalculatorFactory

    algorithm = GLMBackgroundCalculatorFactory.create(data.experiments)
    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r)
        except Exception:
            errors += 1
    assert len(reflections) == 15193
    assert errors == 333


class IntensityCalculatorFactory:
    @staticmethod
    def create(data, detector_space=False, deconvolution=False):
        from dials.algorithms.integration.parallel_integrator import (
            GaussianRSMultiCrystalReferenceProfileData,
            GaussianRSReferenceProfileData,
            ReferenceProfileData,
        )
        from dials.algorithms.profile_model.gaussian_rs.algorithm import (
            GaussianRSIntensityCalculatorFactory,
        )
        from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
        from dials.algorithms.profile_model.modeller import CircleSampler

        reference = data.reference[0]
        experiments = data.experiments

        assert len(reference) % 9 == 0
        num_scan_points = len(reference) // 9

        data_spec = GaussianRSMultiCrystalReferenceProfileData()
        for e in experiments:
            sampler = CircleSampler(
                e.detector[0].get_image_size(),
                e.scan.get_array_range(),
                num_scan_points,
            )

            spec = TransformSpec(
                e.beam,
                e.detector,
                e.goniometer,
                e.scan,
                e.profile.sigma_b(deg=False),
                e.profile.sigma_m(deg=False),
                e.profile.n_sigma() * 1.5,
                5,
            )

            temp = reference

            reference = ReferenceProfileData()
            for d, m in temp:
                reference.append(d, m)

            spec = GaussianRSReferenceProfileData(reference, sampler, spec)

            data_spec.append(spec)

        return GaussianRSIntensityCalculatorFactory.create(
            data_spec, detector_space, deconvolution
        )


def test_gaussianrs_reciprocal_space_intensity_calculator(data):
    algorithm = IntensityCalculatorFactory.create(
        data, detector_space=False, deconvolution=False
    )

    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r, [])
        except Exception:
            errors += 1

    assert len(reflections) == 15193
    assert errors == 5296


def test_gaussianrs_detector_space_intensity_calculator(data):
    algorithm = IntensityCalculatorFactory.create(
        data, detector_space=True, deconvolution=False
    )

    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r, [])
            partiality_old = r.get("partiality_old")
            partiality_new = r.get("partiality")
        except Exception as e:
            if ("Partiality too small to fit" in str(e)) or (
                "DIALS_ASSERT(fit.niter() < 100) failure." in str(e)
            ):
                errors += 1
                continue
            raise

        assert partiality_old < 1.0 and partiality_old >= 0
        assert partiality_new < 1.0 and partiality_new >= 0

    assert len(reflections) == 15193
    assert errors == 4802


def test_gaussianrs_detector_space_with_deconvolution_intensity_calculator(data):
    algorithm = IntensityCalculatorFactory.create(
        data, detector_space=True, deconvolution=True
    )

    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r, [])
            partiality_old = r.get("partiality_old")
            partiality_new = r.get("partiality")
        except Exception as e:
            if ("Partiality too small to fit" in str(e)) or (
                "DIALS_ASSERT(fit.niter() < 100) failure." in str(e)
            ):
                errors += 1
                continue
            raise

        assert partiality_old < 1.0 and partiality_old >= 0
        assert partiality_new < 1.0 and partiality_new >= 0

    assert len(reflections) == 15193
    assert errors == 4802


def test_gaussianrs_detector_space_with_deconvolution_intensity_calculator2(data):
    from scitbx import matrix

    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    R = None
    for r in reflections:
        if r.get("partiality") > 0.9 and r.get("intensity.sum.value") > 100:
            R = r
            break
    assert R is not None

    px = R.get("xyzcal.px")
    mm = R.get("xyzcal.mm")

    px1 = (px[0] - 3, px[1] - 3, px[2])
    px2 = (px[0] + 3, px[1] + 3, px[2])
    mm1 = (mm[0] - 3 * 0.172, mm[1] - 3 * 0.172, mm[2])
    mm2 = (mm[0] + 3 * 0.172, mm[1] + 3 * 0.172, mm[2])

    s11 = matrix.col(
        data.experiments[0].detector[0].get_pixel_lab_coord(px1[0:2])
    ).normalize()
    s12 = matrix.col(
        data.experiments[0].detector[0].get_pixel_lab_coord(px2[0:2])
    ).normalize()

    R1 = R.copy()
    R2 = R.copy()

    R1.set_vec3_double("xyzcal.px", px1)
    R1.set_vec3_double("xyzcal.mm", mm1)
    R1.set_vec3_double("s1", s11)

    R2.set_vec3_double("xyzcal.px", px2)
    R2.set_vec3_double("xyzcal.mm", mm2)
    R2.set_vec3_double("s1", s12)

    compute_intensity = IntensityCalculatorFactory.create(
        data, detector_space=True, deconvolution=False
    )

    compute_intensity(R, [])
    R.set_double("intensity.prf.value.old", R.get("intensity.prf.value"))
    R.set_double("intensity.prf.variance.old", R.get("intensity.prf.variance"))

    R_no_deconvolution = R.copy()
    # print "Partiality", R.get("partiality")
    # print "Partiality Old", R.get("partiality_old")
    # print "Intensity", R.get("intensity.prf.value")
    # print "Intensity Old", R.get("intensity.prf.value.old")

    from dials.algorithms.integration.parallel_integrator import MaskCalculatorFactory

    mask_calculator = MaskCalculatorFactory.create(data.experiments)

    mask_calculator(R1, True)
    mask_calculator(R2, True)

    # print R.get("shoebox").mask.as_numpy_array()[0,:,:]

    compute_intensity = IntensityCalculatorFactory.create(
        data, detector_space=True, deconvolution=True
    )

    compute_intensity(R, [R1, R2])

    R_deconvolution = R.copy()

    partiality1 = R_no_deconvolution.get("partiality")
    partiality2 = R_deconvolution.get("partiality")
    intensity = R_deconvolution.get("intensity.prf.value")
    variance = R_deconvolution.get("intensity.prf.variance")

    assert partiality1 <= partiality2

    assert abs(intensity - 179.05997142249996) < 1e-7
    assert abs(variance - 203.56992949599677) < 1e-7

    # print "Partiality", R.get("partiality")
    # print "Partiality Old", R.get("partiality_old")
    # print "Intensity", R.get("intensity.prf.value")
    # print "Intensity Old", R.get("intensity.prf.value.old")
    # print "Variance", R.get("intensity.prf.variance")
    # print "Variance Old", R.get("intensity.prf.variance.old")

    # print R.get("shoebox").mask.all()
    # sbox = R.get("shoebox")
    # print sbox.mask.count(5) + sbox.mask.count(37) + sbox.mask.count(51)


def test_gaussianrs_profile_data_pickling(data):
    from dials.algorithms.integration.parallel_integrator import (
        GaussianRSMultiCrystalReferenceProfileData,
        GaussianRSReferenceProfileData,
        ReferenceProfileData,
    )
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
    from dials.algorithms.profile_model.modeller import CircleSampler

    reference = data.reference[0]
    experiments = data.experiments

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9

    data_spec = GaussianRSMultiCrystalReferenceProfileData()
    for e in experiments:
        sampler = CircleSampler(
            e.detector[0].get_image_size(), e.scan.get_array_range(), num_scan_points
        )

        spec = TransformSpec(
            e.beam,
            e.detector,
            e.goniometer,
            e.scan,
            e.profile.sigma_b(deg=False),
            e.profile.sigma_m(deg=False),
            e.profile.n_sigma() * 1.5,
            5,
        )

        temp = reference

        reference = ReferenceProfileData()
        for d, m in temp:
            reference.append(d, m)

        spec = GaussianRSReferenceProfileData(reference, sampler, spec)

        data_spec.append(spec)

    s = pickle.dumps(data_spec)

    pickle.loads(s)


def test_gaussianrs_reference_profile_calculator(data):
    from dials.algorithms.profile_model.gaussian_rs.algorithm import (
        GaussianRSReferenceCalculatorFactory,
    )

    algorithm = GaussianRSReferenceCalculatorFactory.create(data.experiments)
    reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

    errors = 0
    for r in reflections:
        try:
            algorithm(r)
        except Exception:
            errors += 1
    assert len(reflections) == 15193
    assert errors == 0

    profiles = algorithm.reference_profiles()

    count = 0
    for i in range(len(profiles)):
        p = profiles[i].reference()
        for j in range(len(p)):
            d = p.data(j)
            p.mask(j)
            if len(d) != 0:
                count += 1

    assert count == 9


def test_job_list():
    from dials.algorithms.integration.parallel_integrator import SimpleBlockList

    jobs = SimpleBlockList((0, 60), 20)

    assert jobs[0] == (0, 20)
    assert jobs[1] == (10, 30)
    assert jobs[2] == (20, 40)
    assert jobs[3] == (30, 50)
    assert jobs[4] == (40, 60)

    for frame in range(0, 15):
        assert jobs.block_index(frame) == 0
    for frame in range(16, 25):
        assert jobs.block_index(frame) == 1
    for frame in range(26, 35):
        assert jobs.block_index(frame) == 2
    for frame in range(36, 45):
        assert jobs.block_index(frame) == 3
    for frame in range(46, 60):
        assert jobs.block_index(frame) == 4


def test_reflection_manager(data):
    from dials.algorithms.integration.parallel_integrator import (
        SimpleBlockList,
        SimpleReflectionManager,
    )

    jobs = SimpleBlockList((0, 60), 20)
    reflections = data.reflections

    manager = SimpleReflectionManager(jobs, reflections, 5)

    def check_job(index):
        r = manager.split(index)
        selection = r.get_flags(r.flags.dont_integrate)
        r1 = r.select(~selection)

        j0, j1 = jobs[index]

        for rr in r1.rows():
            z0, z1 = rr["bbox"][4:6]
            zc = int(math.floor((z0 + z1) / 2.0))
            j = jobs.block_index(zc)
            assert j == index
            assert z0 >= j0 and z1 <= j1

    check_job(0)
    check_job(1)
    check_job(2)
    check_job(3)
    check_job(4)
