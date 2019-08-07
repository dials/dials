from __future__ import absolute_import, division, print_function

from dials_pychef_ext import *

from cctbx.array_family import flex
from iotbx.data_plots import table_data
from libtbx import phil

dose_phil_str = """\
dose {
  remove_gaps = True
    .type = bool

  batch
    .multiple = True
  {
    range = None
      .type = ints(value_min=0, size=2)
    dose_start = None
      .type = float(value_min=0)
    dose_step = None
      .type = float(value_min=0)
  }
}
"""

phil_scope = phil.parse(
    """\
d_min = None
  .type = float(value_min=0)
d_max = None
  .type = float(value_min=0)
min_completeness = None
  .type = float(value_min=0, value_max=1)
  .help = "Minimum value of completeness in outer resolution shell used to "
          "determine suitable resolution cutoff for analysis"
resolution_bins = 8
  .type = int
anomalous = False
  .type = bool
range {
  width = 1
    .type = float(value_min=0)
  min = None
    .type = float(value_min=0)
  max = None
    .type = float(value_min=0)
}
%s
"""
    % dose_phil_str
)


class Statistics(object):
    def __init__(
        self, intensities, dose, n_bins=8, range_min=None, range_max=None, range_width=1
    ):

        if isinstance(dose, flex.double):
            sorted_dose = flex.sorted(dose)
            dd = sorted_dose[1:] - sorted_dose[:-1]
            dd = dd.select(dd > 0)
            if len(dd):
                step_size = flex.min(dd)
                dose /= step_size
            dose = dose.iround()
            if flex.min(dose) == 0:
                dose += 1

        # fix for completeness > 1 if screw axes present
        intensities = intensities.customized_copy(
            space_group_info=intensities.space_group()
            .build_derived_reflection_intensity_group(
                anomalous_flag=intensities.anomalous_flag()
            )
            .info(),
            info=intensities.info(),
        )

        self.intensities = intensities
        self.dose = dose
        self.n_bins = n_bins
        self.range_min = range_min
        self.range_max = range_max
        self.range_width = range_width
        assert self.range_width > 0

        if self.range_min is None:
            self.range_min = flex.min(self.dose) - self.range_width
        if self.range_max is None:
            self.range_max = flex.max(self.dose)
        self.n_steps = 2 + int((self.range_max - self.range_min) - self.range_width)

        sel = (self.dose.as_double() <= self.range_max) & (
            self.dose.as_double() >= self.range_min
        )
        self.dose = self.dose.select(sel)

        self.intensities = self.intensities.select(sel)
        self.d_star_sq = self.intensities.d_star_sq().data()

        self.binner = self.intensities.setup_binner_d_star_sq_step(
            d_star_sq_step=(flex.max(self.d_star_sq) - flex.min(self.d_star_sq) + 1e-8)
            / self.n_bins
        )

        # self.dose /= range_width
        self.dose -= int(self.range_min)

        self.dose = flex.size_t(list(self.dose))

        binner_non_anom = intensities.as_non_anomalous_array().use_binning(self.binner)
        n_complete = flex.size_t(binner_non_anom.counts_complete()[1:-1])

        chef_stats = ChefStatistics(
            intensities.indices(),
            intensities.data(),
            intensities.sigmas(),
            intensities.d_star_sq().data(),
            self.dose,
            n_complete,
            self.binner,
            intensities.space_group(),
            intensities.anomalous_flag(),
            self.n_steps,
        )

        self.iplus_comp_bins = chef_stats.iplus_completeness_bins()
        self.iminus_comp_bins = chef_stats.iminus_completeness_bins()
        self.ieither_comp_bins = chef_stats.ieither_completeness_bins()
        self.iboth_comp_bins = chef_stats.iboth_completeness_bins()
        self.iplus_comp_overall = chef_stats.iplus_completeness()
        self.iminus_comp_overall = chef_stats.iminus_completeness()
        self.ieither_comp_overall = chef_stats.ieither_completeness()
        self.iboth_comp_overall = chef_stats.iboth_completeness()
        self.rcp_bins = chef_stats.rcp_bins()
        self.rcp = chef_stats.rcp()
        self.scp_bins = chef_stats.scp_bins()
        self.scp = chef_stats.scp()
        self.rd = chef_stats.rd()

    def completeness_vs_dose_str(self):

        anomalous = self.intensities.anomalous_flag()

        title = "Completeness vs. dose:"
        graph_names = ["Completeness", "Completeness in resolution shells"]

        if anomalous:
            column_labels = (
                ["Dose"]
                + [
                    "%.2f-%.2f(A)" % self.binner.bin_d_range(i + 1)
                    for i in range(self.n_bins)
                ]
                + ["I+", "I-", "I", "dI"]
            )
            column_formats = (
                ["%8.1f"]
                + ["%5.3f" for i in range(self.n_bins)]
                + ["%5.3f", "%5.3f", "%5.3f", "%5.3f"]
            )
            # graph_columns = [[0,1,2,3,4]]
            graph_columns = [
                [0] + list(range(self.n_bins + 1, self.n_bins + 5)),
                list(range(self.n_bins + 1)),
            ]
        else:
            column_labels = (
                ["Dose"]
                + [
                    "%.2f-%.2f(A)" % self.binner.bin_d_range(i + 1)
                    for i in range(self.n_bins)
                ]
                + ["I"]
            )
            column_formats = (
                ["%8.1f"] + ["%5.3f" for i in range(self.n_bins)] + ["%5.3f"]
            )
            graph_columns = [[0, self.n_bins + 1], list(range(self.n_bins + 1))]

        table_completeness = table_data(
            title=title,
            column_labels=column_labels,
            column_formats=column_formats,
            graph_names=graph_names,
            graph_columns=graph_columns,
        )
        for i in range(self.n_steps):
            if anomalous:
                row = (
                    [i * self.range_width + self.range_min]
                    + [self.ieither_comp_bins[i_bin, i] for i_bin in range(self.n_bins)]
                    + [
                        self.iplus_comp_overall[i],
                        self.iminus_comp_overall[i],
                        self.ieither_comp_overall[i],
                        self.iboth_comp_overall[i],
                    ]
                )
            else:
                row = (
                    [i * self.range_width + self.range_min]
                    + [self.ieither_comp_bins[i_bin, i] for i_bin in range(self.n_bins)]
                    + [self.ieither_comp_overall[i]]
                )
            table_completeness.add_row(row)

        return table_completeness.format_loggraph()

    def rcp_vs_dose_str(self):
        title = "Cumulative radiation damage analysis:"
        column_labels = (
            ["Dose"]
            + [
                "%.2f-%.2f(A)" % self.binner.bin_d_range(i + 1)
                for i in range(self.n_bins)
            ]
            + ["Rcp(d)"]
        )
        column_formats = ["%8.1f"] + ["%7.4f" for i in range(self.n_bins + 1)]
        graph_names = ["Rcp(d)", "Rcp(d), in resolution shells"]
        graph_columns = [[0, self.n_bins + 1], list(range(self.n_bins + 1))]

        table_rcp = table_data(
            title=title,
            column_labels=column_labels,
            column_formats=column_formats,
            graph_names=graph_names,
            graph_columns=graph_columns,
        )
        for i in range(self.n_steps):
            row = (
                [i * self.range_width + self.range_min]
                + [self.rcp_bins[j, i] for j in range(self.binner.n_bins_used())]
                + [self.rcp[i]]
            )
            table_rcp.add_row(row)

        return table_rcp.format_loggraph()

    def scp_vs_dose_str(self):
        title = "Normalised radiation damage analysis:"
        column_labels = (
            ["Dose"]
            + [
                "%.2f-%.2f(A)" % self.binner.bin_d_range(i + 1)
                for i in range(self.n_bins)
            ]
            + ["Rcp(d)"]
        )
        column_formats = ["%8.1f"] + ["%7.4f" for i in range(self.n_bins + 1)]
        graph_names = ["Scp(d)", "Scp(d), in resolution shells"]
        graph_columns = [[0, self.n_bins + 1], list(range(self.n_bins + 1))]

        table_scp = table_data(
            title=title,
            column_labels=column_labels,
            column_formats=column_formats,
            graph_names=graph_names,
            graph_columns=graph_columns,
        )
        for i in range(self.n_steps):
            row = (
                [i * self.range_width + self.range_min]
                + [self.scp_bins[j, i] for j in range(self.binner.n_bins_used())]
                + [self.scp[i]]
            )
            table_scp.add_row(row)

        return table_scp.format_loggraph()

    def rd_vs_dose_str(self):

        title = "R vs. BATCH difference:"
        column_labels = ["BATCH", "Rd"]
        column_formats = ["%8.1f", "%5.3f"]
        graph_names = ["Rd"]
        graph_columns = [[0, 1]]

        table_rd = table_data(
            title=title,
            column_labels=column_labels,
            column_formats=column_formats,
            graph_names=graph_names,
            graph_columns=graph_columns,
        )
        for i in range(self.n_steps):
            row = [i * self.range_width + self.range_min, self.rd[i]]
            table_rd.add_row(row)

        return table_rd.format_loggraph()

    def to_dict(self):
        x = [i * self.range_width + self.range_min for i in range(self.n_steps)]
        scp_data = []
        rcp_data = []
        rd_data = []
        completeness_data = []

        scp_data.append(
            {
                "x": x,
                "y": list(self.scp),
                "type": "scatter",
                "name": "Scp overall",
                "line": {"width": 3},
            }
        )
        rcp_data.append(
            {
                "x": x,
                "y": list(self.rcp),
                "type": "scatter",
                "name": "Rcp overall",
                "line": {"width": 3},
            }
        )
        rd_data.append({"x": x, "y": list(self.rd), "type": "scatter", "name": "Rd"})

        anomalous = self.intensities.anomalous_flag()
        completeness_data.append(
            {
                "x": x,
                "y": list(self.ieither_comp_overall),
                "type": "scatter",
                "name": "I",
                "line": {"width": 3},
            }
        )
        if anomalous:
            completeness_data.append(
                {
                    "x": x,
                    "y": list(self.iboth_comp_overall),
                    "type": "scatter",
                    "name": "dI",
                    "line": {"width": 3},
                }
            )
            completeness_data.append(
                {
                    "x": x,
                    "y": list(self.iplus_comp_overall),
                    "type": "scatter",
                    "name": "I+",
                    "line": {"width": 3},
                }
            )
            completeness_data.append(
                {
                    "x": x,
                    "y": list(self.iminus_comp_overall),
                    "type": "scatter",
                    "name": "I-",
                    "line": {"width": 3},
                }
            )

        if self.binner.n_bins_used() > 1:
            for j in range(self.binner.n_bins_used()):
                bin_range_suffix = " (%.2f-%.2f A)" % self.binner.bin_d_range(j + 1)
                scp_data.append(
                    {
                        "x": x,
                        "y": list(self.scp_bins[j : j + 1, :].as_1d()),
                        "type": "scatter",
                        "name": "Scp" + bin_range_suffix,
                        "line": {"width": 1, "dash": "dot"},
                    }
                )
                rcp_data.append(
                    {
                        "x": x,
                        "y": list(self.rcp_bins[j : j + 1, :].as_1d()),
                        "type": "scatter",
                        "name": "Rcp" + bin_range_suffix,
                        "line": {"width": 1, "dash": "dot"},
                    }
                )

                completeness_data.append(
                    {
                        "x": x,
                        "y": list(self.ieither_comp_bins[j : j + 1, :].as_1d()),
                        "type": "scatter",
                        "name": "I" + bin_range_suffix,
                        "line": {"width": 1, "dash": "dot"},
                    }
                )

        d = {
            "scp_vs_dose": {
                "data": scp_data,
                "layout": {
                    "title": "Scp vs dose",
                    "xaxis": {"title": "Dose"},
                    "yaxis": {"title": "Scp", "rangemode": "tozero"},
                },
            },
            "rcp_vs_dose": {
                "data": rcp_data,
                "layout": {
                    "title": "Rcp vs dose",
                    "xaxis": {"title": "Dose"},
                    "yaxis": {"title": "Rcp", "rangemode": "tozero"},
                },
            },
            "rd_vs_batch_difference": {
                "data": rd_data,
                "layout": {
                    "title": "Rd vs batch difference",
                    "xaxis": {"title": "Batch difference"},
                    "yaxis": {"title": "Rd", "rangemode": "tozero"},
                },
            },
            "completeness_vs_dose": {
                "data": completeness_data,
                "layout": {
                    "title": "Completeness vs dose",
                    "xaxis": {"title": "Dose"},
                    "yaxis": {"title": "Completeness", "rangemode": "tozero"},
                },
            },
        }
        return d


def run(args):
    from iotbx.reflection_file_reader import any_reflection_file

    interp = phil_scope.command_line_argument_interpreter()
    params, unhandled = interp.process_and_fetch(
        args, custom_processor="collect_remaining"
    )
    params = params.extract()

    args = unhandled

    intensities = None
    batches = None
    dose = None

    reader = any_reflection_file(args[0])
    assert reader.file_type() == "ccp4_mtz"
    arrays = reader.as_miller_arrays(merge_equivalents=False)
    for ma in arrays:
        if ma.info().labels == ["BATCH"]:
            batches = ma
        elif ma.info().labels == ["DOSE"]:
            dose = ma
        elif ma.info().labels == ["I", "SIGI"]:
            intensities = ma
        elif ma.info().labels == ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]:
            intensities = ma

    assert intensities is not None
    assert batches is not None
    mtz_object = reader.file_content()

    indices = mtz_object.extract_original_index_miller_indices()
    intensities = intensities.customized_copy(indices=indices)
    batches = batches.customized_copy(indices=indices)

    if params.anomalous:
        intensities = intensities.as_anomalous_array()
        batches = batches.as_anomalous_array()

    if dose is None:
        dose = batches_to_dose(batches.data(), params.dose)
    else:
        dose = dose.data()

    sel = dose > -1
    intensities = intensities.select(sel)
    dose = dose.select(sel)

    if not params.d_min and params.min_completeness:
        params.d_min = resolution_limit(
            mtz_file=args[0],
            min_completeness=params.min_completeness,
            n_bins=params.resolution_bins,
        )
        print("Estimated d_min: %.2f" % params.d_min)
    if params.d_min or params.d_max:
        sel = flex.bool(intensities.size(), True)
        d_spacings = intensities.d_spacings().data()
        if params.d_min:
            sel &= d_spacings >= params.d_min
        if params.d_max:
            sel &= d_spacings <= params.d_max
        intensities = intensities.select(sel)
        dose = dose.select(sel)

    stats = Statistics(
        intensities,
        dose,
        n_bins=params.resolution_bins,
        range_min=params.range.min,
        range_max=params.range.max,
        range_width=params.range.width,
    )

    print(stats.completeness_vs_dose_str())
    print(stats.rcp_vs_dose_str())
    print(stats.scp_vs_dose_str())
    print(stats.rd_vs_dose_str())


def batches_to_dose(batches, params):
    if len(params.batch):
        dose = flex.double(batches.size(), -1)
        for batch in params.batch:
            start = batch.dose_start
            step = batch.dose_step
            for i in range(batch.range[0], batch.range[1] + 1):
                # inclusive range
                dose.set_selected(batches == i, start + step * (i - batch.range[0]))
        dose = dose.iround()
    elif params.remove_gaps:
        dose = remove_batch_gaps(batches)
    else:
        dose = batches
    return dose


def remove_batch_gaps(batches):
    perm = flex.sort_permutation(batches)
    new_batches = flex.int(batches.size(), -1)
    sorted_batches = batches.select(perm)
    curr_batch = -1
    new_batch = -1
    for i, b in enumerate(sorted_batches):
        if b != curr_batch:
            curr_batch = b
            new_batch += 1
        new_batches[perm[i]] = new_batch
    return new_batches


def resolution_limit(mtz_file, min_completeness, n_bins):
    from xia2.Modules.Resolutionizer import resolutionizer, phil_defaults

    params = phil_defaults.extract().resolutionizer
    params.nbins = n_bins
    r = resolutionizer(mtz_file, params)
    return r.resolution_completeness(limit=min_completeness)
