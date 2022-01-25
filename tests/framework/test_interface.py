from __future__ import annotations


def test_after_import_extensions():
    import dials.extensions
    from dials.extensions.dispersion_spotfinder_threshold_ext import (
        DispersionSpotFinderThresholdExt,
    )

    # Check the expected extensions are present for each interface
    from dials.extensions.gaussian_rs_profile_model_ext import GaussianRSProfileModelExt
    from dials.extensions.null_background_ext import NullBackgroundExt
    from dials.extensions.simple_background_ext import SimpleBackgroundExt
    from dials.extensions.simple_centroid_ext import SimpleCentroidExt

    extensions = dials.extensions.ProfileModel.extensions()
    assert GaussianRSProfileModelExt in extensions
    extensions = dials.extensions.SpotFinderThreshold.extensions()
    assert DispersionSpotFinderThresholdExt in extensions
    extensions = dials.extensions.Centroid.extensions()
    assert SimpleCentroidExt in extensions
    extensions = dials.extensions.Background.extensions()
    assert NullBackgroundExt in extensions
    assert SimpleBackgroundExt in extensions

    # Check phil scope
    phil_scope = dials.extensions.ProfileModel.phil_scope()
    assert phil_scope
    phil_scope = dials.extensions.SpotFinderThreshold.phil_scope()
    assert phil_scope
    phil_scope = dials.extensions.Centroid.phil_scope()
    assert phil_scope
    phil_scope = dials.extensions.Background.phil_scope()
    assert phil_scope
