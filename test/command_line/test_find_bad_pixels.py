from __future__ import absolute_import, division, print_function

import os
import pickle

from dials.command_line import find_bad_pixels
from dials.command_line import dials_import
from dxtbx.model.experiment_list import ExperimentList


def test_x4wide(dials_data, run_in_tmpdir, capsys):
    filepaths = [f.strpath for f in dials_data("x4wide").listdir(sort=True)]
    find_bad_pixels.run(filepaths + ["nproc=2", "print_values=True"])
    assert os.path.exists("pixels.mask")
    with open("pixels.mask", "rb") as fh:
        mask = pickle.load(fh)
    out, _ = capsys.readouterr()
    assert "mask[1482, 1484] = 8" in out
    assert mask.all() == (2527, 2463)
    assert mask.count(False) == 27
    # Test that the output pixels.mask is compatible with dials.import
    dials_import.Script().run(filepaths + ["lookup.mask=pixels.mask"])
    expts = ExperimentList.from_file("imported.expt")
    assert expts[0].imageset.external_lookup.mask.filename == run_in_tmpdir.join(
        "pixels.mask"
    )
    assert len(expts[0].imageset.external_lookup.mask.data) == 1
    assert expts[0].imageset.external_lookup.mask.data[0].data() == mask


def test_x4wide_images(dials_data, run_in_tmpdir, capsys):
    filepaths = dials_data("x4wide").listdir(sort=True)
    find_bad_pixels.run(
        [f.strpath for f in filepaths]
        + ["images=1,10,20,30,40,50,60,70,80,90", "mask=mypixels.mask"]
    )
    assert os.path.exists("mypixels.mask")
    with open("mypixels.mask", "rb") as fh:
        mask = pickle.load(fh)
    out, _ = capsys.readouterr()
    assert "mask[" not in out
    assert "= 8" not in out
    assert mask.all() == (2527, 2463)
    assert mask.count(False) == 50


def test_x4wide_image_range(dials_data, run_in_tmpdir):
    filepaths = dials_data("x4wide").listdir(sort=True)
    find_bad_pixels.run([f.strpath for f in filepaths] + ["image_range=1,30"])
    assert os.path.exists("pixels.mask")
    with open("pixels.mask", "rb") as fh:
        mask = pickle.load(fh)
    assert mask.all() == (2527, 2463)
    assert mask.count(False) == 82
