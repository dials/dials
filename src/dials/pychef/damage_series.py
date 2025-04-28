from __future__ import annotations

import logging
from copy import copy
from math import ceil, floor

import iotbx
from cctbx import crystal, miller, uctbx
from dxtbx.model import ExperimentList
from scitbx.array_family import flex

from dials.command_line.symmetry import median_unit_cell
from dials.pychef import interpret_images_to_doses_options
from dials.report.plots import d_star_sq_to_d_ticks
from dials.util.export_mtz import MTZWriterBase
from dials.util.filter_reflections import filter_reflection_table

logger = logging.getLogger("dials.command_line.damage_analysis")


def _generate_blank_plots_template(label):
    d = {
        "i_over_sigma": {
            "data": [],
            "layout": {
                "title": "<I/σ(I)> vs resolution<br>" + f"{label}",
                "xaxis": {
                    "title": "Resolution (Å)",
                    "tickvals": None,
                    "ticktext": None,
                },
                "yaxis": {"title": "<I/σ(I)>", "rangemode": "tozero"},
            },
        },
        "r_pim": {
            "data": [],
            "layout": {
                "title": "R<sub>pim</sub> vs resolution<br>" + f"{label}",
                "xaxis": {
                    "title": "Resolution (Å)",
                    "tickvals": None,
                    "ticktext": None,
                },
                "yaxis": {"title": "R<sub>pim</sub>", "rangemode": "tozero"},
            },
        },
        "completeness": {
            "data": [],
            "layout": {
                "title": "Completeness vs resolution<br>" + f"{label}",
                "xaxis": {
                    "title": "Resolution (Å)",
                    "tickvals": None,
                    "ticktext": None,
                },
                "yaxis": {"title": "Completeness", "range": (0, 1)},
            },
        },
        "cc_half": {
            "data": [],
            "layout": {
                "title": "CC<sub>½</sub> vs resolution<br>" + f"{label}",
                "xaxis": {
                    "title": "Resolution (Å)",
                    "tickvals": None,
                    "ticktext": None,
                },
                "yaxis": {"title": "CC<sub>½</sub>", "range": (0, 1)},
            },
        },
    }
    return d


class DamageSeriesPlots:
    def __init__(self, d_max=None, d_min=None):
        self._d_star_sq_tickvals = None
        self._d_star_sq_ticktext = None
        self.d_star_sq_bins = None
        self.d_max = d_max
        self.d_min = d_min
        self._damage_series = _generate_blank_plots_template("damage group series")
        self._accumulation_series = _generate_blank_plots_template(
            "damage accumulation series"
        )

    @property
    def damage_series(self):
        return self._damage_series

    @property
    def accumulation_series(self):
        return self._accumulation_series

    def add_to_damage_series(self, data, low, upper):
        self._add_data_to_plots_dict(self._damage_series, data, low, upper)

    def add_to_accumulation_series(self, data, low, upper):
        self._add_data_to_plots_dict(self._accumulation_series, data, low, upper)

    def _add_data_to_plots_dict(self, plots_dict, data, low, upper):
        label = f"{low} <= dose < {upper}"
        try:
            result = iotbx.merging_statistics.dataset_statistics(
                i_obs=data,
                n_bins=20,
                anomalous=False,
                use_internal_variance=False,
                eliminate_sys_absent=False,
                d_min=self.d_min,
                d_max=self.d_max,
            )
        except Exception:
            pass
        else:
            if not self.d_star_sq_bins:
                # when the first data is added, use this to define the overall data range
                self.d_star_sq_bins = [
                    0.5
                    * (uctbx.d_as_d_star_sq(b.d_max) + uctbx.d_as_d_star_sq(b.d_min))
                    for b in result.bins
                ]
                self.d_star_sq_tickvals, self.d_star_sq_ticktext = d_star_sq_to_d_ticks(
                    self.d_star_sq_bins, nticks=5
                )
                for k in self._damage_series.values():
                    k["layout"]["xaxis"]["tickvals"] = self.d_star_sq_tickvals
                    k["layout"]["xaxis"]["ticktext"] = self.d_star_sq_ticktext
                for k in self._accumulation_series.values():
                    k["layout"]["xaxis"]["tickvals"] = self.d_star_sq_tickvals
                    k["layout"]["xaxis"]["ticktext"] = self.d_star_sq_ticktext
                if not self.d_max:
                    self.d_max = result.bins[0].d_max
                if not self.d_min:
                    self.d_min = result.bins[-1].d_min
            d_star_sq_bins_this = [
                0.5 * (uctbx.d_as_d_star_sq(b.d_max) + uctbx.d_as_d_star_sq(b.d_min))
                for b in result.bins
            ]
            plots_dict["r_pim"]["data"].append(
                {
                    "x": d_star_sq_bins_this,
                    "y": [b.r_pim for b in result.bins],
                    "type": "scatter",
                    "name": label,
                    "mode": "lines",
                }
            )
            plots_dict["cc_half"]["data"].append(
                {
                    "x": d_star_sq_bins_this,
                    "y": [b.cc_one_half for b in result.bins],
                    "type": "scatter",
                    "name": label,
                    "mode": "lines",
                }
            )
            plots_dict["i_over_sigma"]["data"].append(
                {
                    "x": d_star_sq_bins_this,
                    "y": [b.i_over_sigma_mean for b in result.bins],
                    "type": "scatter",
                    "name": label,
                    "mode": "lines",
                }
            )
            plots_dict["completeness"]["data"].append(
                {
                    "x": d_star_sq_bins_this,
                    "y": [b.completeness for b in result.bins],
                    "type": "scatter",
                    "name": label,
                    "mode": "lines",
                }
            )


def _refl_to_miller_array(reflection_table, new_experiments):
    logger2 = logging.getLogger("dials")
    logger2.disabled = True
    reflection_table = filter_reflection_table(
        reflection_table, intensity_choice=["scale"], partiality_threshold=0.4
    )
    logger2.disabled = False
    intensities = miller.array(
        miller.set(
            crystal.symmetry(
                unit_cell=median_unit_cell(new_experiments),
                space_group=new_experiments[0].crystal.get_space_group(),
            ),
            indices=reflection_table["miller_index"],
            anomalous_flag=False,
        ),
        data=reflection_table["intensity.scale.value"],
        sigmas=flex.sqrt(reflection_table["intensity.scale.variance"]),
    )
    intensities.set_observation_type_xray_intensity()
    return intensities


def _write_mtz(sel_intensities, sel_doses, fname):
    writer = MTZWriterBase(sel_intensities.space_group())
    writer.add_crystal(unit_cell=sel_intensities.unit_cell())
    writer.add_empty_dataset(wavelength=sel_intensities.info().wavelength)

    nref = sel_intensities.size()
    writer.mtz_file.adjust_column_array_sizes(nref)
    writer.mtz_file.set_n_reflections(nref)
    type_table = {"H": "H", "K": "H", "L": "H", "M_ISYM": "Y"}
    # assign H, K, L, M_ISYM space
    for column in "H", "K", "L", "M_ISYM":
        writer.current_dataset.add_column(column, type_table[column]).set_values(
            flex.double(nref, 0.0).as_float()
        )
    writer.mtz_file.replace_original_index_miller_indices(sel_intensities.indices())
    writer.current_dataset.add_column("I", "J").set_values(
        sel_intensities.data().as_float()
    )
    writer.current_dataset.add_column("SIGI", "Q").set_values(
        sel_intensities.sigmas().as_float()
    )
    writer.current_dataset.add_column("DOSE", "R").set_values(
        sel_doses.as_double().as_float()
    )
    writer.mtz_file.write(fname)
    logger.info(f"Saved {nref} reflections to {fname}")


def generate_damage_series_mtz(params, doses, intensities):
    plots = DamageSeriesPlots(d_max=params.d_max, d_min=params.d_min)
    group_size = params.damage_series.dose_group_size
    assert group_size > 0.0

    max_dose = flex.max(doses)
    min_dose = flex.min(doses)
    dose_range = max_dose - min_dose + 1
    n_output = ceil(dose_range / group_size)

    dose_boundaries = [ceil(n * group_size) for n in range(n_output + 1)]

    for n in range(n_output):
        lower_dose_boundary = dose_boundaries[n]
        upper_dose_boundary = dose_boundaries[n + 1]
        fname = f"damage_series_{lower_dose_boundary}_{upper_dose_boundary}.mtz"

        sel = (doses >= lower_dose_boundary) & (doses < upper_dose_boundary)
        sel_intensities = intensities.select(sel)
        sel_intensities.set_info(intensities.info())
        plots.add_to_damage_series(
            sel_intensities, lower_dose_boundary, upper_dose_boundary
        )
        if params.output.damage_series:
            _write_mtz(sel_intensities, doses.select(sel), fname)

    lower_dose_boundary = 0
    for n in range(n_output):
        upper_dose_boundary = dose_boundaries[n + 1]
        fname = f"damage_series_{lower_dose_boundary}_{upper_dose_boundary}.mtz"
        sel = (doses >= lower_dose_boundary) & (doses < upper_dose_boundary)
        sel_intensities = intensities.select(sel)
        sel_intensities.set_info(intensities.info())
        plots.add_to_accumulation_series(
            sel_intensities, lower_dose_boundary, upper_dose_boundary
        )
        if params.output.accumulation_series:
            if n == 0 and params.output.damage_series:
                continue
            _write_mtz(sel_intensities, doses.select(sel), fname)
    return plots


def generate_damage_series(params, experiments, reflection_table):
    # first set up plotting stuff
    plots = DamageSeriesPlots(d_max=params.d_max, d_min=params.d_min)

    group_size = params.damage_series.dose_group_size
    assert group_size > 0.0

    reflection_table = reflection_table.select(experiments)

    doses = flex.double()
    start_doses, doses_per_image = interpret_images_to_doses_options(
        experiments,
        params.dose.experiments.dose_per_image,
        params.dose.experiments.starting_doses,
        params.dose.experiments.shared_crystal,
    )

    for expt, starting_dose, dose_per_img in zip(
        experiments, start_doses, doses_per_image
    ):
        refls = reflection_table.select(expt)
        imgno = flex.floor(refls["xyzobs.px.value"].parts()[2])
        dose = (imgno * dose_per_img) + starting_dose
        doses.extend(dose)

    # e.g. start dose of 0, images have doses of 0 -> 499 for the purposes
    # of this calculation

    doses = doses.iround()

    max_dose = flex.max(doses)
    min_dose = flex.min(doses)
    dose_range = max_dose - min_dose + 1
    n_output = ceil(dose_range / group_size)

    def select_data_in_dose_range(lower_dose_boundary, upper_dose_boundary):
        sel = (doses >= lower_dose_boundary) & (doses < upper_dose_boundary)
        refl = reflection_table.select(sel)

        sel = (doses >= lower_dose_boundary) & (doses < upper_dose_boundary)
        refl = reflection_table.select(sel)

        # now need to get the subsets of the experiments corresponding to that
        # dose range

        new_expts = ExperimentList()
        for expt, start_dose, dose_per_image in zip(
            experiments, start_doses, doses_per_image
        ):
            imgrange = expt.scan.get_image_range()
            n_images = (
                imgrange[1] - imgrange[0] + 1
            )  # i.e. 1->500 image range has 500 images
            end_dose = start_dose + (dose_per_image * (n_images - 1))
            if end_dose < lower_dose_boundary:
                continue
            elif start_dose >= upper_dose_boundary:
                continue
            else:
                first_image_in_dose_range = max(
                    [
                        floor((lower_dose_boundary - start_dose) / dose_per_image),
                        imgrange[0] - 1,
                    ]
                )
                last_image_in_dose_range_plus_one = min(
                    [
                        ceil((upper_dose_boundary - start_dose) / dose_per_image),
                        imgrange[1],
                    ]
                )  # use this to slice as upper half-bound
                new_expt = copy(expt)
                new_expt.scan = new_expt.scan[
                    first_image_in_dose_range:last_image_in_dose_range_plus_one
                ]
                # if new_expt.crystal.scan_varying_model
                if new_expt.scaling_model:
                    new_expt.scaling_model.limit_image_range(
                        (first_image_in_dose_range, last_image_in_dose_range_plus_one)
                    )
                new_expts.append(new_expt)

        return new_expts, refl

    dose_boundaries = [ceil(n * group_size) for n in range(n_output + 1)]

    for n in range(n_output):
        lower_dose_boundary = dose_boundaries[n]
        upper_dose_boundary = dose_boundaries[n + 1]

        name = f"damage_series_{lower_dose_boundary}_{upper_dose_boundary}"
        reflections_filename = name + ".refl"
        experiments_filename = name + ".expt"

        new_expts, refl = select_data_in_dose_range(
            lower_dose_boundary,
            upper_dose_boundary,
        )

        assert len(new_expts) == len(refl.experiment_identifiers().keys()), (
            f"{len(new_expts)} != {list(refl.experiment_identifiers().keys())}"
        )
        if params.output.damage_series:
            logger.info(
                f"Saving experimental data for range {lower_dose_boundary} <= dose < {upper_dose_boundary}"
            )
            refl.as_file(reflections_filename)
            logger.info(
                f"Saved {refl.size()} reflections from {len(new_expts)} experiments to {reflections_filename}"
            )
            new_expts.as_json(experiments_filename)
            logger.info(f"Saved corresponding experiments to {experiments_filename}")

        intensities = _refl_to_miller_array(refl, new_expts)
        # now add to plots data
        plots.add_to_damage_series(
            intensities, lower_dose_boundary, upper_dose_boundary
        )

    lower_dose_boundary = 0
    for n in range(n_output):
        upper_dose_boundary = dose_boundaries[n + 1]
        name = f"damage_series_{lower_dose_boundary}_{upper_dose_boundary}"
        reflections_filename = name + ".refl"
        experiments_filename = name + ".expt"

        new_expts, refl = select_data_in_dose_range(
            lower_dose_boundary,
            upper_dose_boundary,
        )

        assert len(new_expts) == len(refl.experiment_identifiers().keys()), (
            f"{len(new_expts)} != {list(refl.experiment_identifiers().keys())}"
        )
        if params.output.accumulation_series:
            if n == 0 and params.output.damage_series:
                continue
            logger.info(
                f"Saving experimental data for range {lower_dose_boundary} <= dose < {upper_dose_boundary}"
            )
            refl.as_file(reflections_filename)
            logger.info(
                f"Saved {refl.size()} reflections from {len(new_expts)} experiments to {reflections_filename}"
            )
            new_expts.as_json(experiments_filename)
            logger.info(f"Saved corresponding experiments to {experiments_filename}")

        intensities = _refl_to_miller_array(refl, new_expts)
        plots.add_to_accumulation_series(
            intensities, lower_dose_boundary, upper_dose_boundary
        )

    return plots
