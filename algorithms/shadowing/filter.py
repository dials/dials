from __future__ import absolute_import, division, print_function


def filter_shadowed_reflections(experiments, reflections, experiment_goniometer=False):
    from dxtbx.masking import is_inside_polygon
    from scitbx.array_family import flex

    shadowed = flex.bool(reflections.size(), False)
    for expt_id in range(len(experiments)):
        expt = experiments[expt_id]
        imageset = expt.imageset
        masker = imageset.masker()
        detector = expt.detector
        sel = reflections["id"] == expt_id
        isel = sel.iselection()
        x, y, z = reflections["xyzcal.px"].select(isel).parts()
        start, end = expt.scan.get_array_range()
        for i in range(start, end):
            shadow = masker.project_extrema(
                detector, expt.scan.get_angle_from_array_index(i)
            )
            img_sel = (z >= i) & (z < (i + 1))
            img_isel = img_sel.iselection()
            for p_id in range(len(detector)):
                panel = reflections["panel"].select(img_isel)
                if shadow[p_id].size() < 4:
                    continue
                panel_isel = img_isel.select(panel == p_id)
                inside = is_inside_polygon(
                    shadow[p_id],
                    flex.vec2_double(
                        x.select(isel.select(panel_isel)),
                        y.select(isel.select(panel_isel)),
                    ),
                )
                shadowed.set_selected(panel_isel, inside)

    return shadowed
