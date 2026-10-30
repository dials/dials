from __future__ import annotations

import pickle
import random


def test_pickle():
    from scitbx.array_family import flex

    from dials.model.data import PixelList

    size = (100, 100)
    sf = 10
    image = flex.double(flex.grid(size))
    mask = flex.bool(flex.grid(size))
    for i in range(len(image)):
        image[i] = random.randint(0, 100)
        mask[i] = bool(random.randint(0, 1))
    pl = PixelList(sf, image, mask)
    assert pl.size() == size
    assert pl.frame() == sf

    obj = pickle.dumps(pl)
    pl2 = pickle.loads(obj)

    assert pl2.size() == size
    assert pl2.frame() == sf
    assert len(pl2) == len(pl)
    assert pl2.index().all_eq(pl.index())
    assert pl2.value().all_eq(pl.value())


def test_add_image():
    from scitbx.array_family import flex

    from dials.model.data import PixelList, PixelListLabeller

    size = (2000, 2000)
    sf = 10
    labeller = PixelListLabeller()

    count = 0
    for i in range(3):
        image = flex.random_int_gaussian_distribution(size[0] * size[1], 100, 5)
        mask = flex.random_bool(size[0] * size[1], 0.5)
        image.reshape(flex.grid(size))
        mask.reshape(flex.grid(size))
        pl = PixelList(sf + i, image, mask)
        count += len(mask.as_1d().select(mask.as_1d()))
        labeller.add(pl)
    assert len(labeller.values()) == count


def test_labels_3d():
    from scitbx.array_family import flex

    from dials.model.data import PixelList, PixelListLabeller

    size = (500, 500)
    sf = 0
    labeller = PixelListLabeller()

    count = 0
    mask_list = []
    for i in range(3):
        image = flex.random_int_gaussian_distribution(size[0] * size[1], 100, 5)
        mask = flex.random_bool(size[0] * size[1], 0.5)
        image.reshape(flex.grid(size))
        mask.reshape(flex.grid(size))
        pl = PixelList(sf + i, image, mask)
        count += len(mask.as_1d().select(mask.as_1d()))
        labeller.add(pl)
        mask_list.append(mask)

    coords = labeller.coords()
    labels = labeller.labels_3d()

    # Create a map of labels
    label_map = flex.int(flex.grid(3, size[0], size[1]))
    for c, l in zip(coords, labels):
        label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(3):
        for j in range(size[0]):
            for i in range(size[1]):
                if mask_list[k][j, i]:
                    l1 = labels[vi]
                    if k > 0 and mask_list[k - 1][j, i]:
                        l2 = label_map[k - 1, j, i]
                        assert l2 == l1
                    if j > 0 and mask_list[k][j - 1, i]:
                        l2 = label_map[k, j - 1, i]
                        assert l2 == l1
                    if i > 0 and mask_list[k][j, i - 1]:
                        l2 = label_map[k, j, i - 1]
                        assert l2 == l1
                    vi += 1


def test_labels_2d():
    from scitbx.array_family import flex

    from dials.model.data import PixelList, PixelListLabeller

    size = (500, 500)
    sf = 0
    labeller = PixelListLabeller()

    count = 0
    mask_list = []
    for i in range(3):
        image = flex.random_int_gaussian_distribution(size[0] * size[1], 100, 5)
        mask = flex.random_bool(size[0] * size[1], 0.5)
        image.reshape(flex.grid(size))
        mask.reshape(flex.grid(size))
        pl = PixelList(sf + i, image, mask)
        count += len(mask.as_1d().select(mask.as_1d()))
        labeller.add(pl)
        mask_list.append(mask)

    coords = labeller.coords()
    labels = labeller.labels_2d()

    # Create a map of labels
    label_map = flex.int(flex.grid(3, size[0], size[1]))
    for c, l in zip(coords, labels):
        label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(3):
        for j in range(size[0]):
            for i in range(size[1]):
                if mask_list[k][j, i]:
                    l1 = labels[vi]
                    if k > 0 and mask_list[k - 1][j, i]:
                        l2 = label_map[k - 1, j, i]
                        assert l2 != l1
                    if j > 0 and mask_list[k][j - 1, i]:
                        l2 = label_map[k, j - 1, i]
                        assert l2 == l1
                    if i > 0 and mask_list[k][j, i - 1]:
                        l2 = label_map[k, j, i - 1]
                        assert l2 == l1
                    vi += 1


def test_with_no_points():
    from scitbx.array_family import flex

    from dials.model.data import PixelList, PixelListLabeller

    size = (500, 500)
    sf = 0
    labeller = PixelListLabeller()

    count = 0
    mask_list = []
    for i in range(3):
        image = flex.random_int_gaussian_distribution(size[0] * size[1], 100, 5)
        mask = flex.bool(size[0] * size[0], False)
        image.reshape(flex.grid(size))
        mask.reshape(flex.grid(size))
        pl = PixelList(sf + i, image, mask)
        count += len(mask.as_1d().select(mask.as_1d()))
        labeller.add(pl)
        mask_list.append(mask)

    coords = labeller.coords()
    labels1 = labeller.labels_2d()
    labels2 = labeller.labels_2d()

    assert len(coords) == 0
    assert len(labels1) == 0
    assert len(labels2) == 0


def two_spots_on_two_frames():
    from dials.array_family import flex
    from dials.model.data import PixelList, PixelListLabeller

    size = (20, 20)
    labeller = PixelListLabeller()
    for frame in range(2):
        image = flex.double(flex.grid(size), 0)
        mask = flex.bool(flex.grid(size), False)
        for j, i in ((2, 2), (2, 3), (3, 2), (3, 3)):
            image[j, i] = 10 + frame
            mask[j, i] = True
            image[j + 8, i + 10] = 20 + frame
            mask[j + 8, i + 10] = True
        labeller.add(PixelList(frame, image, mask))
    return labeller


def check_spot_shoeboxes(shoeboxes):
    assert len(shoeboxes) == 2
    intensities = sorted(s.summed_intensity().observed.value for s in shoeboxes)
    assert intensities == [4 * 10 + 4 * 11, 4 * 20 + 4 * 21]
    for shoebox in shoeboxes:
        assert shoebox.zsize() == 2
        assert shoebox.is_data_allocated()
        assert not shoebox.is_background_allocated()
        assert shoebox.is_consistent()


def test_pixel_list_shoeboxes_have_no_background():
    from dials.array_family import flex

    creator = flex.PixelListShoeboxCreator(
        two_spots_on_two_frames(), min_pixels=1, max_pixels=100
    )
    check_spot_shoeboxes(creator.result())


def test_combined_strong_spots_have_no_background():
    from dials.algorithms.spot_finding import StrongSpotCombiner
    from dials.array_family import flex

    creator = flex.PixelListShoeboxCreator(
        two_spots_on_two_frames(), twod=True, min_pixels=1, max_pixels=100
    )
    shoeboxes = creator.result()
    first_frame = flex.bool(s.bbox[4] == 0 for s in shoeboxes)
    combiner = StrongSpotCombiner()
    combiner.add(shoeboxes.select(first_frame))
    combiner.add(shoeboxes.select(~first_frame))
    check_spot_shoeboxes(combiner.shoeboxes())
