from __future__ import absolute_import, division, print_function

from dials.array_family import flex
from libtbx.phil import parse

phil_scope = parse(
    """
overlaps_filter {
  foreground_foreground {
    enable = False
      .type = bool
      .help = "Remove all spots in which neighbors' foreground"
              "impinges on the spot's foreground"
  }
  foreground_background {
    enable = False
      .type = bool
      .help = "Remove all spots in which neighbors' foreground"
              "impinges on the spot's background"
  }
}
""",
    process_includes=True,
)


class OverlapsFilter(object):
    from dials.algorithms.shoebox import MaskCode

    code_fgd = MaskCode.Foreground | MaskCode.Valid
    code_bgd = MaskCode.Background | MaskCode.Valid

    def is_fgd(self, code):
        return (code & self.code_fgd) == self.code_fgd

    def is_bgd(self, code):
        return (code & self.code_bgd) == self.code_bgd

    def __init__(self, refl, expt):
        self.refl = refl
        self.expt = expt
        self.masks = {}
        det = self.expt.detector
        assert len(det) == 1  # for now
        self.size_fast, self.size_slow = det[0].get_image_size()
        self.array_size = self.size_fast * self.size_slow

    def create_simple_mask(self):
        self.masks["simple_mask"] = flex.size_t(self.array_size)
        for obs in self.refl:
            shoebox = obs["shoebox"]
            fast_coords = range(shoebox.xsize())
            slow_coords = range(shoebox.ysize())
            for f, s in zip(fast_coords, slow_coords):
                f_abs = f + shoebox.bbox[0]  # relative to detector
                s_abs = s + shoebox.bbox[2]  # relative to detector
                posn = f_abs + s_abs * self.size_fast  # position in mask array
                posn_in_shoebox = f + shoebox.xsize() * s  # position in shoebox
                try:
                    self.masks["simple_mask"][posn] |= shoebox.mask[posn_in_shoebox]
                except IndexError:  # bbox may extend past detector limits
                    continue

    def create_referenced_mask(self, test_code, mask_name):
        self.masks[mask_name] = [flex.size_t() for _ in range(self.array_size)]
        for idx in range(len(self.refl)):
            obs = self.refl[idx]
            shoebox = obs["shoebox"]
            fast_coords = range(shoebox.xsize())
            slow_coords = range(shoebox.ysize())
            for f, s in zip(fast_coords, slow_coords):
                f_abs = f + shoebox.bbox[0]  # relative to detector
                s_abs = s + shoebox.bbox[2]  # relative to detector
                posn = f_abs + s_abs * self.size_fast  # position in mask array
                posn_in_shoebox = f + shoebox.xsize() * s  # position in shoebox
                if (shoebox.mask[posn_in_shoebox] & test_code) == test_code:
                    try:
                        self.masks[mask_name][posn].append(idx)
                    except IndexError:  # bbox may extend past detector limits
                        continue

    def filter_using_simple_mask(self, mask_lambda, shoebox_lambda=lambda x: True):
        """At each pixel, examine the simple mask to determine if contributing
        observations should be excluded. When this condition mask_lambda is met,
        for each contributing observation, use optional condition shoebox_lambda
        to determine if the observation should be excluded (e.g. to exclude only
        those reflections contributing foreground). Return the mask reflecting
        this filter.
        """
        keep_refl_bool = flex.bool(len(self.refl), True)
        for idx in range(len(self.refl)):
            obs = self.refl[idx]
            shoebox = obs["shoebox"]
            fast_coords = range(shoebox.xsize())
            slow_coords = range(shoebox.ysize())
            for f, s in zip(fast_coords, slow_coords):
                f_abs = f + shoebox.bbox[0]  # relative to detector
                s_abs = s + shoebox.bbox[2]  # relative to detector
                posn = f_abs + s_abs * self.size_fast  # position in mask array
                posn_in_shoebox = f + shoebox.xsize() * s  # position in shoebox
                try:
                    if mask_lambda(
                        self.masks["simple_mask"][posn]
                    ):  # condition met in simple mask
                        if shoebox_lambda(
                            shoebox.mask[posn_in_shoebox]
                        ):  # condition met in shoebox
                            keep_refl_bool[idx] = False
                except IndexError:  # bbox may extend past detector limits
                    continue
        return keep_refl_bool

    def filter_all_using_referenced_mask(self, mask_name):
        """Return the mask reflecting the exclusion of any reflections for which the
        mask condition is true (e.g. untrusted pixels).
        """
        keep_refl_bool = flex.bool(len(self.refl), True)
        for i in self.masks[mask_name]:
            if len(i) > 0:
                for ref in i:
                    keep_refl_bool[ref] = False
        return keep_refl_bool

    def filter_overlaps_using_referenced_mask(self, mask_name):
        """At each pixel, define an overlap to be more than one reference (to an
        observation) indicated in the mask. Return the mask reflecting the exclusion
        of any overlaps (e.g. foreground with foreground).
        """
        keep_refl_bool = flex.bool(len(self.refl), True)
        for i in self.masks[mask_name]:
            if len(i) > 1:
                for ref in i:
                    keep_refl_bool[ref] = False
        return keep_refl_bool

    def remove_foreground_foreground_overlaps(self):
        self.create_referenced_mask(self.code_fgd, "foreground")
        self.refl = self.refl.select(
            self.filter_overlaps_using_referenced_mask("foreground")
        )

    def remove_foreground_background_overlaps(self):
        self.create_simple_mask()

        def is_overlap(code):
            return self.is_fgd(code) and self.is_bgd(code)

        self.refl = self.refl.select(
            self.filter_using_simple_mask(mask_lambda=is_overlap)
        )


class OverlapsFilterMultiExpt(object):
    def __init__(self, refl, expt):
        self.filters = [
            OverlapsFilter(r, e) for (r, e) in zip(refl.split_by_experiment_id(), expt)
        ]

    def remove_foreground_foreground_overlaps(self):
        for f in self.filters:
            f.remove_foreground_foreground_overlaps()

    def remove_foreground_background_overlaps(self):
        for f in self.filters:
            f.remove_foreground_background_overlaps()

    @property
    def refl(self):
        rlist = [f.refl for f in self.filters]
        r0 = flex.reflection_table()
        for r in rlist:
            r0.extend(r)
        return r0

    @property
    def expt(self):
        elist = [f.expt for f in self.filters]
        from dxtbx.model.experiment_list import ExperimentList

        e0 = ExperimentList()
        for e in elist:
            e0.append(e)
        return e0
