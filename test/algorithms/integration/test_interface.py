from __future__ import absolute_import, division, print_function

import random

import pytest


def test_split_blocks_1_frame():
    from dials.array_family import flex
    from dials.algorithms.integration.integrator import JobList

    r = flex.reflection_table()
    r["value1"] = flex.double()
    r["value2"] = flex.int()
    r["value3"] = flex.double()
    r["bbox"] = flex.int6()
    r["id"] = flex.int()
    expected = []
    for i in range(100):
        x0 = random.randint(0, 100)
        x1 = x0 + random.randint(1, 10)
        y0 = random.randint(0, 100)
        y1 = y0 + random.randint(1, 10)
        z0 = random.randint(0, 100)
        z1 = z0 + random.randint(1, 10)
        v1 = random.uniform(0, 100)
        v2 = random.randint(0, 100)
        v3 = random.uniform(0, 100)
        r.append(
            {
                "id": 0,
                "value1": v1,
                "value2": v2,
                "value3": v3,
                "bbox": (x0, x1, y0, y1, z0, z1),
            }
        )
        for z in range(z0, z1):
            expected.append(
                {
                    "id": 0,
                    "value1": v1,
                    "value2": v2,
                    "value3": v3,
                    "bbox": (x0, x1, y0, y1, z, z + 1),
                    "partial_id": i,
                }
            )

    jobs = JobList()
    jobs.add((0, 1), (0, 111), 1)

    jobs.split(r)
    assert len(r) == len(expected)
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
        assert r1["bbox"] == r2["bbox"]
        assert r1["partial_id"] == r2["partial_id"]
        assert abs(r1["value1"] - r2["value1"]) < EPS
        assert r1["value2"] == r2["value2"]
        assert abs(r1["value3"] - r2["value3"]) < EPS


def test_split_blocks_non_overlapping():
    from dials.array_family import flex
    from dials.algorithms.integration.integrator import JobList
    from scitbx.array_family import shared

    blocks = shared.tiny_int_2(
        [
            (0, 10),
            (10, 20),
            (20, 30),
            (30, 35),
            (35, 40),
            (40, 50),
            (50, 60),
            (60, 70),
            (70, 80),
            (80, 90),
            (90, 100),
            (100, 110),
        ]
    )

    jobs = JobList((0, 1), blocks)

    r = flex.reflection_table()
    r["value1"] = flex.double()
    r["value2"] = flex.int()
    r["value3"] = flex.double()
    r["bbox"] = flex.int6()
    r["id"] = flex.int()
    expected = []
    for i in range(100):
        x0 = random.randint(0, 100)
        x1 = x0 + random.randint(1, 10)
        y0 = random.randint(0, 100)
        y1 = y0 + random.randint(1, 10)
        z0 = random.randint(0, 100)
        z1 = z0 + random.randint(1, 10)
        v1 = random.uniform(0, 100)
        v2 = random.randint(0, 100)
        v3 = random.uniform(0, 100)
        r.append(
            {
                "id": 0,
                "value1": v1,
                "value2": v2,
                "value3": v3,
                "bbox": (x0, x1, y0, y1, z0, z1),
            }
        )

        for j in range(len(blocks)):
            b0 = blocks[j][0]
            b1 = blocks[j][1]
            if (
                (z0 >= b0 and z1 <= b1)
                or (z0 < b1 and z1 >= b1)
                or (z0 < b0 and z1 > b0)
            ):
                z00 = max(b0, z0)
                z11 = min(b1, z1)
                expected.append(
                    {
                        "id": 0,
                        "value1": v1,
                        "value2": v2,
                        "value3": v3,
                        "bbox": (x0, x1, y0, y1, z00, z11),
                        "partial_id": i,
                    }
                )

    jobs.split(r)
    assert len(r) == len(expected)
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
        assert r1["bbox"] == r2["bbox"]
        assert r1["partial_id"] == r2["partial_id"]
        assert abs(r1["value1"] - r2["value1"]) < EPS
        assert r1["value2"] == r2["value2"]
        assert abs(r1["value3"] - r2["value3"]) < EPS


def test_split_blocks_overlapping():
    from dials.array_family import flex
    from dials.algorithms.integration.integrator import JobList
    from scitbx.array_family import shared

    blocks = shared.tiny_int_2(
        [
            (0, 10),
            (5, 15),
            (10, 20),
            (15, 25),
            (20, 30),
            (25, 35),
            (30, 40),
            (35, 45),
            (40, 50),
            (45, 55),
            (50, 60),
            (55, 65),
            (60, 70),
            (65, 75),
            (70, 80),
            (75, 85),
            (80, 90),
            (85, 95),
            (90, 100),
            (95, 105),
            (100, 110),
        ]
    )

    jobs = JobList((0, 1), blocks)

    r = flex.reflection_table()
    r["value1"] = flex.double()
    r["value2"] = flex.int()
    r["value3"] = flex.double()
    r["bbox"] = flex.int6()
    r["id"] = flex.int()
    expected = []
    for i in range(100):
        x0 = random.randint(0, 100)
        x1 = x0 + random.randint(1, 10)
        y0 = random.randint(0, 100)
        y1 = y0 + random.randint(1, 10)
        z0 = random.randint(0, 90)
        z1 = z0 + random.randint(1, 20)
        v1 = random.uniform(0, 100)
        v2 = random.randint(0, 100)
        v3 = random.uniform(0, 100)
        r.append(
            {
                "id": 0,
                "value1": v1,
                "value2": v2,
                "value3": v3,
                "bbox": (x0, x1, y0, y1, z0, z1),
            }
        )
        expected.append(
            {
                "id": 0,
                "value1": v1,
                "value2": v2,
                "value3": v3,
                "bbox": (x0, x1, y0, y1, z0, z1),
            }
        )

    jobs.split(r)
    assert len(r) > 100
    for r1 in r:
        v1 = r1["value1"]
        v2 = r1["value2"]
        v3 = r1["value3"]
        bbox = r1["bbox"]
        pid = r1["partial_id"]

        z0 = bbox[4]
        z1 = bbox[5]
        success = False
        for i in range(len(blocks)):
            b0 = blocks[i][0]
            b1 = blocks[i][1]
            if z0 >= b0 and z1 <= b1:
                success = True
                break
        assert success

        v11 = expected[pid]["value1"]
        v22 = expected[pid]["value2"]
        v33 = expected[pid]["value3"]
        bb = expected[pid]["bbox"]
        assert v11 == v1
        assert v22 == v2
        assert v33 == v3
        assert bb[0] == bbox[0]
        assert bb[1] == bbox[1]
        assert bb[2] == bbox[2]
        assert bb[3] == bbox[3]


def test_reflection_manager():
    from dials.array_family import flex

    reflections = flex.reflection_table()
    reflections["panel"] = flex.size_t()
    reflections["bbox"] = flex.int6()
    reflections["miller_index"] = flex.miller_index()
    reflections["s1"] = flex.vec3_double()
    reflections["xyzcal.px"] = flex.vec3_double()
    reflections["xyzcal.mm"] = flex.vec3_double()
    reflections["entering"] = flex.bool()
    reflections["id"] = flex.int()
    reflections["flags"] = flex.size_t()

    width = 1000
    height = 1000
    nrefl = 10000
    array_range = (0, 130)
    block_size = 20

    random.seed(0)
    processed = [[] for i in range(12)]
    for i in range(nrefl):
        x0 = random.randint(0, width - 10)
        y0 = random.randint(0, height - 10)
        zs = random.randint(2, 9)
        x1 = x0 + random.randint(2, 10)
        y1 = y0 + random.randint(2, 10)
        for k, j in enumerate([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]):
            m = k + i * 12
            pos = random.choice(["left", "right", "centre"])
            if pos == "left":
                z0 = j - zs
                z1 = j
            elif pos == "right":
                z0 = j
                z1 = j + zs
            else:
                z0 = j - zs // 2
                z1 = j + zs // 2
            bbox = (x0, x1, y0, y1, z0, z1)
            reflections.append(
                {
                    "panel": random.randint(0, 1),
                    "bbox": bbox,
                    "flags": flex.reflection_table.flags.reference_spot,
                }
            )
            processed[k].append(m)

        # Add reflection to ignore
        # zc = random.choice([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120])
        # z0 = zc - 11
        # z1 = zc + 11
        # bbox = (x0, x1, y0, y1, z0, z1)
        # reflections.append({
        #   "panel" : randint(0,1),
        #   "bbox" : bbox,
        #   "flags" : flex.reflection_table.flags.reference_spot
        # })

    from dials.algorithms.integration.integrator import ReflectionManager
    from dials.algorithms.integration.integrator import JobList

    jobs = JobList()
    jobs.add((0, 1), array_range, block_size)

    # Create the executor
    executor = ReflectionManager(jobs, reflections)

    # Ensure the tasks make sense
    jobs = [executor.job(i) for i in range(len(executor))]
    assert len(executor) == 12
    assert not executor.finished()
    assert len(jobs) == 12
    assert jobs[0].frames() == (0, 20)
    assert jobs[1].frames() == (10, 30)
    assert jobs[2].frames() == (20, 40)
    assert jobs[3].frames() == (30, 50)
    assert jobs[4].frames() == (40, 60)
    assert jobs[5].frames() == (50, 70)
    assert jobs[6].frames() == (60, 80)
    assert jobs[7].frames() == (70, 90)
    assert jobs[8].frames() == (80, 100)
    assert jobs[9].frames() == (90, 110)
    assert jobs[10].frames() == (100, 120)
    assert jobs[11].frames() == (110, 130)

    # Get the task specs
    data0 = executor.split(0)
    data1 = executor.split(1)
    data2 = executor.split(2)
    data3 = executor.split(3)
    data4 = executor.split(4)
    data5 = executor.split(5)
    data6 = executor.split(6)
    data7 = executor.split(7)
    data8 = executor.split(8)
    data9 = executor.split(9)
    data10 = executor.split(10)
    data11 = executor.split(11)
    assert len(data0) == len(processed[0])
    assert len(data1) == len(processed[1])
    assert len(data2) == len(processed[2])
    assert len(data3) == len(processed[3])
    assert len(data4) == len(processed[4])
    assert len(data5) == len(processed[5])
    assert len(data6) == len(processed[6])
    assert len(data7) == len(processed[7])
    assert len(data8) == len(processed[8])
    assert len(data9) == len(processed[9])
    assert len(data10) == len(processed[10])
    assert len(data11) == len(processed[11])

    # Add some results
    data0["data"] = flex.double(len(data0), 1)
    data1["data"] = flex.double(len(data1), 2)
    data2["data"] = flex.double(len(data2), 3)
    data3["data"] = flex.double(len(data3), 4)
    data4["data"] = flex.double(len(data4), 5)
    data5["data"] = flex.double(len(data5), 6)
    data6["data"] = flex.double(len(data6), 7)
    data7["data"] = flex.double(len(data7), 8)
    data8["data"] = flex.double(len(data8), 9)
    data9["data"] = flex.double(len(data9), 10)
    data10["data"] = flex.double(len(data10), 11)
    data11["data"] = flex.double(len(data11), 12)

    # Accumulate the data again
    assert not executor.finished()
    executor.accumulate(0, data0)
    executor.accumulate(1, data1)
    executor.accumulate(2, data2)
    executor.accumulate(3, data3)
    executor.accumulate(4, data4)
    executor.accumulate(5, data5)
    executor.accumulate(6, data6)
    executor.accumulate(7, data7)
    executor.accumulate(8, data8)
    executor.accumulate(9, data9)
    executor.accumulate(10, data10)
    executor.accumulate(11, data11)
    assert executor.finished()

    # Get results and check they're as expected
    data = executor.data()
    result = data["data"]
    for i in range(len(processed)):
        for j in range(len(processed[i])):
            assert result[processed[i][j]] == i + 1

    # Test passed


@pytest.mark.parametrize("nproc", [1, 2])
def test_integrator_3d(dials_data, nproc):
    from dxtbx.model.experiment_list import ExperimentListFactory
    from dials.algorithms.profile_model.gaussian_rs import Model
    from dials.array_family import flex
    from math import pi

    path = dials_data("centroid_test_data").join("experiments.json").strpath

    exlist = ExperimentListFactory.from_json_file(path)
    exlist[0].profile = Model(
        None, n_sigma=3, sigma_b=0.024 * pi / 180.0, sigma_m=0.044 * pi / 180.0
    )

    rlist = flex.reflection_table.from_predictions(exlist[0])
    rlist["id"] = flex.int(len(rlist), 0)
    rlist.compute_bbox(exlist)
    rlist.compute_zeta_multi(exlist)
    rlist.compute_d(exlist)

    from dials.algorithms.integration.integrator import Integrator3D
    from dials.algorithms.integration.integrator import phil_scope
    from libtbx.phil import parse

    params = phil_scope.fetch(
        parse(
            """
    integration.block.size=%d
    integration.mp.nproc=%d
    integration.profile_fitting=False
  """
            % (5, nproc)
        )
    ).extract()

    integrator = Integrator3D(exlist, rlist, params)
    integrator.integrate()


def test_summation(dials_data):
    from dxtbx.model.experiment_list import ExperimentListFactory
    from dials.algorithms.profile_model.gaussian_rs import Model
    from dials.array_family import flex
    from math import pi

    path = dials_data("centroid_test_data").join("experiments.json").strpath

    exlist = ExperimentListFactory.from_json_file(path)
    exlist[0].profile = Model(
        None, n_sigma=3, sigma_b=0.024 * pi / 180.0, sigma_m=0.044 * pi / 180.0
    )

    rlist = flex.reflection_table.from_predictions(exlist[0])
    rlist["id"] = flex.int(len(rlist), 0)

    def integrate(integrator_type, rlist):
        from dials.algorithms.integration.integrator import IntegratorFactory
        from dials.algorithms.integration.integrator import (
            phil_scope as master_phil_scope,
        )
        from libtbx.phil import parse

        rlist = rlist.copy()

        phil_scope = parse(
            """
      integration.background.algorithm=null
      integration.intensity.algorithm=sum
      integration.intensity.sum.integrator=%s
      integration.block.size=0.5
      integration.profile_fitting=False
    """
            % integrator_type
        )

        params = master_phil_scope.fetch(source=phil_scope).extract()

        integrator = IntegratorFactory.create(params, exlist, rlist)

        result = integrator.integrate()
        return result

    from libtbx.test_utils import approx_equal

    def approx_equal_dict(a, b, k):
        return approx_equal(a[k], b[k])

    # Do summation by all different methods
    result1 = integrate("3d", rlist)
    result2 = integrate("flat3d", rlist)
    result3 = integrate("2d", rlist)
    result4 = integrate("single2d", rlist)
    assert len(result1) >= len(rlist)
    assert len(result2) >= len(rlist)
    assert len(result3) >= len(rlist)
    assert len(result4) >= len(rlist)

    # result1 and result2 should be the same
    assert len(result1) == len(result2)
    for r1, r2 in zip(result1, result2):
        assert r1["partial_id"] == r2["partial_id"]
        assert r1["bbox"] == r2["bbox"]
        assert r1["entering"] == r2["entering"]
        assert r1["flags"] == r2["flags"]
        assert r1["id"] == r2["id"]
        assert r1["miller_index"] == r2["miller_index"]
        assert r1["panel"] == r2["panel"]
        assert approx_equal_dict(r1, r2, "d")
        assert approx_equal_dict(r1, r2, "intensity.sum.value")
        assert approx_equal_dict(r1, r2, "intensity.sum.variance")
        assert approx_equal_dict(r1, r2, "lp")
        assert approx_equal_dict(r1, r2, "partiality")
        assert approx_equal_dict(r1, r2, "s1")
        assert approx_equal_dict(r1, r2, "xyzcal.mm")
        assert approx_equal_dict(r1, r2, "xyzcal.px")
        assert approx_equal_dict(r1, r2, "zeta")

    # result3 and result4 should be the same
    assert len(result3) == len(result4)
    for r3, r4 in zip(result3, result4):
        assert r3["partial_id"] == r4["partial_id"]
        assert r3["bbox"] == r4["bbox"]
        assert r3["entering"] == r4["entering"]
        assert r3["flags"] == r4["flags"]
        assert r3["id"] == r4["id"]
        assert r3["miller_index"] == r4["miller_index"]
        assert r3["panel"] == r4["panel"]
        assert approx_equal_dict(r3, r4, "d")
        assert approx_equal_dict(r3, r4, "intensity.sum.value")
        assert approx_equal_dict(r3, r4, "intensity.sum.variance")
        assert approx_equal_dict(r3, r4, "lp")
        assert approx_equal_dict(r3, r4, "partiality")
        assert approx_equal_dict(r3, r4, "s1")
        assert approx_equal_dict(r3, r4, "xyzcal.mm")
        assert approx_equal_dict(r3, r4, "xyzcal.px")
        assert approx_equal_dict(r3, r4, "xyzobs.px.value")
        assert approx_equal_dict(r3, r4, "xyzobs.px.variance")
        assert approx_equal_dict(r3, r4, "zeta")

    # result3 should add up to result1
    assert len(result3) >= len(result1)
    expected1 = rlist.copy()
    expected1["intensity.sum.value"] = flex.double(len(rlist), 0)
    expected1["intensity.sum.variance"] = flex.double(len(rlist), 0)
    for r1 in result1:
        pid = r1["partial_id"]
        r2 = expected1[pid]
        assert r1["entering"] == r2["entering"]
        assert r1["id"] == r2["id"]
        assert r1["miller_index"] == r2["miller_index"]
        assert r1["panel"] == r2["panel"]
        assert approx_equal_dict(r1, r2, "s1")
        assert approx_equal_dict(r1, r2, "xyzcal.mm")
        assert approx_equal_dict(r1, r2, "xyzcal.px")
        expected1["intensity.sum.value"][pid] += r1["intensity.sum.value"]
        expected1["intensity.sum.variance"][pid] += r1["intensity.sum.variance"]
    expected3 = rlist.copy()
    expected3["intensity.sum.value"] = flex.double(len(rlist), 0)
    expected3["intensity.sum.variance"] = flex.double(len(rlist), 0)
    for r1 in result3:
        pid = r1["partial_id"]
        r2 = expected3[pid]
        assert r1["entering"] == r2["entering"]
        assert r1["id"] == r2["id"]
        assert r1["miller_index"] == r2["miller_index"]
        assert r1["panel"] == r2["panel"]
        assert approx_equal_dict(r1, r2, "s1")
        assert approx_equal_dict(r1, r2, "xyzcal.mm")
        assert approx_equal_dict(r1, r2, "xyzcal.px")
        expected3["intensity.sum.value"][pid] += r1["intensity.sum.value"]
        expected3["intensity.sum.variance"][pid] += r1["intensity.sum.variance"]
    for r1, r3 in zip(expected1, expected3):
        assert approx_equal_dict(r1, r3, "intensity.sum.value")
        assert approx_equal_dict(r1, r3, "intensity.sum.variance")
