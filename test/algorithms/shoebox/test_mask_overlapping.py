from __future__ import absolute_import, division, print_function


def predict_reflections(sweep, crystal):
    from dials.algorithms import shoebox
    from dials.array_family import flex
    from dxtbx.model.experiment_list import ExperimentList
    from dxtbx.model.experiment_list import Experiment
    from dials.algorithms.profile_model.gaussian_rs import Model

    # Get models from the sweep
    beam = sweep.get_beam()
    detector = sweep.get_detector()
    gonio = sweep.get_goniometer()
    scan = sweep.get_scan()

    sigma_b = beam.get_sigma_divergence(deg=True)
    sigma_m = crystal.get_mosaicity(deg=True)

    exlist = ExperimentList()
    exlist.append(
        Experiment(
            imageset=sweep,
            beam=beam,
            detector=detector,
            goniometer=gonio,
            scan=scan,
            crystal=crystal,
            profile=Model(None, 3, sigma_b, sigma_m, deg=True),
        )
    )

    predicted = flex.reflection_table.from_predictions(exlist[0])
    predicted["id"] = flex.int(len(predicted), 0)
    predicted.compute_bbox(exlist)

    # Find overlapping reflections
    overlaps = shoebox.find_overlapping(predicted["bbox"])

    # Return the reflections and overlaps
    return predicted, overlaps


def test(dials_data):
    from dxtbx.serialize import load
    from dials.algorithms import shoebox
    from dials.array_family import flex

    # Load the sweep and crystal
    sweep = load.imageset(dials_data("centroid_test_data").join("sweep.json").strpath)
    crystal = load.crystal(
        dials_data("centroid_test_data").join("crystal.json").strpath
    )

    # Get models from the sweep
    beam = sweep.get_beam()
    detector = sweep.get_detector()
    gonio = sweep.get_goniometer()
    scan = sweep.get_scan()

    # Get the reflections and overlaps
    reflections, adjacency_list = predict_reflections(sweep, crystal)
    reflections["shoebox"] = flex.shoebox(reflections["panel"], reflections["bbox"])
    reflections["shoebox"].allocate_with_value(shoebox.MaskCode.Valid)

    # If the adjacency list is given, then create the reflection mask
    assert len(detector) == 1
    image_size = detector[0].get_image_size()
    shoeboxes = reflections["shoebox"]
    coords = reflections["xyzcal.px"]
    shoebox_masker = shoebox.MaskOverlapping()
    shoebox_masker(shoeboxes, coords, adjacency_list)

    # Loop through all edges
    overlapping = []
    for e in adjacency_list.edges():
        v1, v2 = adjacency_list.source(e), adjacency_list.target(e)
        overlapping.append(v1)
        overlapping.append(v2)

    # Ensure elements are unique
    overlapping = set(overlapping)

    # Ensure we have some overlaps
    assert len(overlapping) > 0

    # Get all non-overlapping reflections
    all_r = set(range(len(reflections)))
    non_overlapping = all_r.difference(overlapping)

    # Run the tests
    tst_non_overlapping(reflections, non_overlapping, detector[0].get_image_size())
    tst_overlapping(reflections, overlapping, adjacency_list, image_size)


def tst_non_overlapping(reflections, non_overlapping, image_size):
    """Ensure non-overlapping reflections have all their values 1."""
    from dials.algorithms import shoebox

    # Check that all elements in non_overlapping masks are 1
    shoeboxes = reflections["shoebox"]
    for i in non_overlapping:
        mask = shoeboxes[i].mask
        assert mask.all_eq(shoebox.MaskCode.Valid)


def tst_overlapping(reflections, overlapping, adjacency_list, image_size):
    """Ensure masks for overlapping reflections are set properly."""
    import numpy
    from scitbx import matrix
    from dials.algorithms import shoebox

    # Loop through all overlaps
    shoeboxes = reflections["shoebox"]
    coord = reflections["xyzcal.px"]
    for i in overlapping:
        r1 = shoeboxes[i]
        bbox_1 = r1.bbox
        r1_coord = matrix.col(coord[i])

        # Create a mask that we expect
        r1_size = (bbox_1[5] - bbox_1[4], bbox_1[3] - bbox_1[2], bbox_1[1] - bbox_1[0])
        expected_mask = numpy.zeros(shape=r1_size, dtype=numpy.int32)
        expected_mask[:, :, :] = shoebox.MaskCode.Valid

        # Loop through all reflections which this reflection overlaps
        for j in adjacency_list.adjacent_vertices(i):
            r2 = shoeboxes[j]
            bbox_2 = r2.bbox
            r2_coord = matrix.col(coord[j])

            # Get bounding box of intersection
            bbox_3 = (
                max(bbox_1[0], bbox_2[0]),
                min(bbox_1[1], bbox_2[1]),
                max(bbox_1[2], bbox_2[2]),
                min(bbox_1[3], bbox_2[3]),
                max(bbox_1[4], bbox_2[4]),
                min(bbox_1[5], bbox_2[5]),
            )

            # Check intersection is valid
            assert bbox_3[0] < bbox_3[1]
            assert bbox_3[2] < bbox_3[3]
            assert bbox_3[4] < bbox_3[5]

            # Get the coordinates are all mask values
            mask_coord = []
            for k in range(bbox_3[4], bbox_3[5]):
                for j in range(bbox_3[2], bbox_3[3]):
                    for i in range(bbox_3[0], bbox_3[1]):
                        mask_coord.append(matrix.col((i + 0.5, j + 0.5, k + 0.5)))

            def dist(a, m):
                return numpy.array([(a - b).length() for b in m])

            # Find the indices in the intersection area where r2 is closer to
            # the point than r1
            ind = numpy.where(dist(r1_coord, mask_coord) > dist(r2_coord, mask_coord))[
                0
            ]

            # Set the mask values for r1 where r2 is closer to 0
            k0, k1 = bbox_3[4] - bbox_1[4], bbox_3[5] - bbox_1[4]
            j0, j1 = bbox_3[2] - bbox_1[2], bbox_3[3] - bbox_1[2]
            i0, i1 = bbox_3[0] - bbox_1[0], bbox_3[1] - bbox_1[0]
            intersect_mask = expected_mask[k0:k1, j0:j1, i0:i1]
            intersect_mask_1d = intersect_mask.reshape(-1)
            intersect_mask_1d[ind] = 0
            intersect_mask[:, :] = intersect_mask_1d.reshape(intersect_mask.shape)
            expected_mask[k0:k1, j0:j1, i0:i1] = intersect_mask

        # Check the masks are the same
        calculated_mask = r1.mask.as_numpy_array()
        assert numpy.all(calculated_mask == expected_mask)
