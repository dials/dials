from __future__ import annotations

import pytest


@pytest.mark.xfail(
    reason="Test does not appear correct. See https://github.com/dials/dials/issues/2909"
)
def test_for_overlaps(dials_data):
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
    refl_path = str(
        dials_data("trypsin_multi_lattice") / "refined.refl",
    )
    expt_path = str(
        dials_data("trypsin_multi_lattice") / "refined.expt",
    )

    importer = Importer(
        [refl_path, expt_path],
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    reflections = flatten_reflections(importer.reflections)
    experiments = flatten_experiments(importer.experiments)

    from dials.algorithms.integration.overlaps_filter import OverlapsFilterMultiExpt

    # Take just first two experiments
    refl = reflections[0]
    refl = refl.select((refl["id"] == 0) | (refl["id"] == 1))
    experiments = experiments[0:2]

    overlaps_filter = OverlapsFilterMultiExpt(refl, experiments)
    overlaps_filter.remove_foreground_foreground_overlaps()
    overlaps_filter.remove_foreground_background_overlaps()
    reflections = [overlaps_filter.refl]

    for expt, refl in zip(experiments, reflections):
        det = expt.detector
        size_fast, size_slow = det[0].get_image_size()
        mask_array = flex.size_t(size_fast * size_slow)
        for obs in refl.rows():
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
                ) + "observed centroid (%d, %d) on image %d" % (
                    obs["xyzcal.px"][0],
                    obs["xyzcal.px"][1],
                    obs["xyzcal.px"][2],
                )
                try:
                    mask_array[posn] |= shoebox.mask[posn_in_shoebox]
                except IndexError:
                    continue
        for i, this_code in enumerate(mask_array):
            assert not is_overlap(this_code), (
                "Overlapping foreground and background found at (%d, %d)"
                % (
                    i % shoebox.xsize(),
                    i // shoebox.xsize(),
                )
            )
