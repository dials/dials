from __future__ import annotations

import pickle
import random
import subprocess
import sys

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.model.data import ImageVolume, MultiPanelImageVolume


def make_image_volume(frame0, frame1, height, width):
    volume = ImageVolume(frame0, frame1, height, width)
    for frame in range(frame0, frame1):
        # Integer-valued data so that comparisons after float storage are exact
        data = flex.double([random.randint(0, 100) for _ in range(height * width)])
        data.reshape(flex.grid(height, width))
        mask = flex.bool([random.random() < 0.8 for _ in range(height * width)])
        mask.reshape(flex.grid(height, width))
        volume.set_image(frame, data, mask)
    return volume


def assert_image_volumes_equal(a, b):
    assert b.frame0() == a.frame0()
    assert b.frame1() == a.frame1()
    for array in ("data", "background", "mask"):
        x, y = getattr(a, array)(), getattr(b, array)()
        assert y.all() == x.all()
        assert list(y) == list(x)
    assert b.is_consistent()


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_image_volume_pickle(protocol):
    volume = make_image_volume(5, 9, 6, 7)
    volume2 = pickle.loads(pickle.dumps(volume, protocol=protocol))
    assert_image_volumes_equal(volume, volume2)


def test_image_volume_pickle_is_independent_copy():
    volume = make_image_volume(0, 2, 3, 4)
    volume2 = pickle.loads(pickle.dumps(volume))

    # Modifying the original must not affect the unpickled copy
    data = flex.double(flex.grid(3, 4), 999)
    volume.set_image(0, data, flex.bool(flex.grid(3, 4), True))
    assert volume.data()[0] == 999
    assert volume2.data()[0] != 999


def test_unpickled_image_volume_is_usable():
    volume = make_image_volume(0, 3, 4, 5)
    volume2 = pickle.loads(pickle.dumps(volume))

    bbox = (1, 4, 0, 3, 1, 3)
    assert list(volume2.extract_data(bbox)) == list(volume.extract_data(bbox))

    data = flex.double(flex.grid(4, 5), 42)
    volume2.set_image(2, data, flex.bool(flex.grid(4, 5), True))
    assert list(volume2.data()[40:]) == [42] * 20
    assert volume2.is_consistent()


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_multi_panel_image_volume_pickle(protocol):
    multi = MultiPanelImageVolume()
    multi.add(make_image_volume(2, 5, 4, 5))
    multi.add(make_image_volume(2, 5, 6, 3))

    multi2 = pickle.loads(pickle.dumps(multi, protocol=protocol))
    assert len(multi2) == len(multi)
    assert multi2.frame0() == multi.frame0()
    assert multi2.frame1() == multi.frame1()
    for i in range(len(multi)):
        assert_image_volumes_equal(multi.get(i), multi2.get(i))


def test_empty_multi_panel_image_volume_pickle():
    multi2 = pickle.loads(pickle.dumps(MultiPanelImageVolume()))
    assert len(multi2) == 0


def test_image_volume_pickle_preserves_labels(dials_data):
    # Labels record which reflections own each foreground pixel, so are only
    # populated once a mask has been computed into the volume
    data_dir = dials_data("l_cysteine_dials_output")
    experiments = ExperimentListFactory.from_json_file(
        data_dir / "11_integrated.expt", check_format=False
    )
    reflections = flex.reflection_table.from_file(data_dir / "11_integrated.refl")
    reflections = reflections.select(reflections["id"] >= 0)

    frame0, frame1 = 0, 5
    nx, ny = experiments[0].detector[0].get_image_size()
    x0, x1, y0, y1, z0, z1 = reflections["bbox"].parts()
    reflections = reflections.select(
        (x0 >= 0) & (y0 >= 0) & (x1 <= nx) & (y1 <= ny) & (z0 >= 0) & (z1 <= frame1)
    )
    assert len(reflections) > 0

    volume = ImageVolume(frame0, frame1, ny, nx)
    for frame in range(frame0, frame1):
        volume.set_image(
            frame, flex.double(flex.grid(ny, nx), 1), flex.bool(flex.grid(ny, nx), True)
        )
    multi = MultiPanelImageVolume()
    multi.add(volume)
    reflections.compute_mask(experiments=experiments, image_volume=multi)

    volume2 = pickle.loads(pickle.dumps(volume))
    num_foreground = 0
    for index, bbox in enumerate(reflections["bbox"]):
        mask = volume.extract_mask(bbox, index)
        assert list(volume2.extract_mask(bbox, index)) == list(mask)
        num_foreground += sum(1 for m in mask if m & MaskCode.Foreground)
    assert num_foreground > 0


def test_image_volume_setstate_rejects_wrong_shapes():
    volume = make_image_volume(0, 2, 3, 4)
    state = make_image_volume(0, 2, 4, 3).__getstate__()
    with pytest.raises(RuntimeError):
        volume.__setstate__(state)


def test_unpickle_without_importing_dials(tmp_path):
    multi = MultiPanelImageVolume()
    multi.add(make_image_volume(0, 2, 3, 4))
    (tmp_path / "volume.pickle").write_bytes(pickle.dumps(multi))

    # Unpickling only imports the extension module, not dials.array_family
    script = (
        "import pickle, sys; "
        "v = pickle.load(open(sys.argv[1], 'rb')); "
        "assert 'dials.array_family' not in sys.modules; "
        "assert len(v) == 1 and v.get(0).is_consistent(); "
        "print(sum(v.get(0).data()))"
    )
    result = subprocess.run(
        [sys.executable, "-c", script, tmp_path / "volume.pickle"],
        capture_output=True,
        check=True,
    )
    assert float(result.stdout) == sum(multi.get(0).data())
