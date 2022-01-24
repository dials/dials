from __future__ import annotations

import dials_algorithms_centroid_simple_ext


def centroid(experiments, reflections, image_volume=None):
    """
    Do the centroiding

    :param experiments: The experiment list
    :param reflections: The reflection list
    """

    # Create a centroider instance
    centroider = dials_algorithms_centroid_simple_ext.Centroider()

    # Add all the experiments
    for exp in experiments:
        if exp.scan is not None:
            centroider.add(exp.detector, exp.scan)
        else:
            centroider.add(exp.detector)

    if image_volume is None:
        return centroider(reflections)
    return centroider(reflections, image_volume)
