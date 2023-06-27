from __future__ import annotations

from dxtbx.model.experiment_list import ExperimentListFactory
from iotbx.detectors import FlexImage

from dials.util import export_bitmaps


def test_imageset_as_flex_image(dials_data):
    expts = ExperimentListFactory.from_filenames(
        dials_data("centroid_test_data", pathlib=True).glob("centroid_000*.cbf")
    )
    imageset = expts[0].imageset
    flex_images = list(
        export_bitmaps.imageset_as_flex_image(imageset, images=(1, 3, 5))
    )
    assert len(flex_images) == 3
    for img in flex_images:
        assert isinstance(img, FlexImage)
