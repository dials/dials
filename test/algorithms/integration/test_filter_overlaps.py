from __future__ import absolute_import, division, print_function

import os


def test_for_overlaps(dials_regression):
    from cctbx.array_family import flex
    from dials.algorithms.shoebox import MaskCode

    code_fgd = MaskCode.Foreground | MaskCode.Valid

    def is_fgd(code):
        return (code & code_fgd) == code_fgd

    code_bgd = MaskCode.Background | MaskCode.Valid

    code_overlap = code_fgd | code_bgd

    def is_overlap(code):
        return (code & code_overlap) == code_overlap

    from dials.util.options import Importer, flatten_experiments, flatten_reflections

    # test data
    refl_path = os.path.join(
        dials_regression,
        "integration_test_data",
        "stills_PSII",
        "idx-20161021225550223_integrated.pickle",
    )
    expt_path = os.path.join(
        dials_regression,
        "integration_test_data",
        "stills_PSII",
        "idx-20161021225550223_refined_experiments.json",
    )

    importer = Importer(
        [refl_path, refl_path, expt_path, expt_path],
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    reflections = flatten_reflections(importer.reflections)
    experiments = flatten_experiments(importer.experiments)

    from dials.algorithms.integration.overlaps_filter import OverlapsFilterMultiExpt

    overlaps_filter = OverlapsFilterMultiExpt(reflections[0], experiments)
    overlaps_filter.remove_foreground_foreground_overlaps()
    overlaps_filter.remove_foreground_background_overlaps()
    reflections = [overlaps_filter.refl]

    for expt, refl in zip(experiments, reflections):
        det = expt.detector
        size_fast, size_slow = det[0].get_image_size()
        mask_array = flex.size_t(size_fast * size_slow)
        for obs in refl:
            shoebox = obs["shoebox"]
            fast_coords = range(shoebox.xsize())
            slow_coords = range(shoebox.ysize())
            for f, s in zip(fast_coords, slow_coords):
                f_abs = f + shoebox.bbox[0]  # relative to detector
                s_abs = s + shoebox.bbox[2]  # relative to detector
                posn = f_abs + s_abs * size_fast  # position in mask_array
                posn_in_shoebox = f + shoebox.xsize() * s  # position in shoebox
                assert not (
                    is_fgd(shoebox.mask[posn_in_shoebox]) and is_fgd(mask_array[posn])
                ), "Overlapping foreground found at indexed position (%d, %d), " % (
                    f_abs,
                    s_abs,
                ) + "observed centroid (%d, %d)" % (
                    obs["xyzcal.px"][0],
                    obs["xyzcal.px"][1],
                )
                try:
                    mask_array[posn] |= shoebox.mask[posn_in_shoebox]
                except IndexError:
                    continue
        for i, this_code in enumerate(mask_array):
            assert not is_overlap(this_code), (
                "Overlapping foreground and background found at (%d, %d)"
                % (i % shoebox.xsize(), i // shoebox.xsize())
            )
