from __future__ import annotations

from os.path import join

from dxtbx.model.experiment_list import ExperimentListFactory

from dials.array_family import flex
from dials_tof_scaling_ext import (
    TOFCorrectionsData,
    tof_extract_shoeboxes_to_reflection_table,
)


def test_tof_extract_shoeboxes(dials_data):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_data("isis_sxd_nacl_processed", pathlib=True), "imported.expt")
    )
    reflections = flex.reflection_table.from_msgpack_file(
        join(dials_data("isis_sxd_nacl_processed", pathlib=True), "strong.refl")
    )

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )

    expt_data = experiments[0].imageset

    ## Shoeboxes with no corrections

    tof_extract_shoeboxes_to_reflection_table(
        reflections, experiments[0], expt_data, False
    )

    ## Shoeboxes with Lorentz correction

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    tof_extract_shoeboxes_to_reflection_table(
        reflections, experiments[0], expt_data, True
    )

    ## Shoeboxes with incident/empty run normalisation

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    experiment_cls = experiments[0].imageset.get_format_class()
    incident_run_file = join(
        dials_data("isis_sxd_example_data", pathlib=True), "sxd_vanadium_run.nxs"
    )
    empty_run_file = join(
        dials_data("isis_sxd_example_data", pathlib=True), "sxd_empty_run.nxs"
    )
    incident_fmt_class = experiment_cls.get_instance(incident_run_file)
    empty_fmt_class = experiment_cls.get_instance(empty_run_file)

    incident_data = experiment_cls(incident_run_file).get_imageset(incident_run_file)
    empty_data = experiment_cls(empty_run_file).get_imageset(empty_run_file)
    incident_proton_charge = incident_fmt_class.get_proton_charge()
    empty_proton_charge = empty_fmt_class.get_proton_charge()
    expt_proton_charge = experiment_cls.get_instance(
        experiments[0].imageset.paths()[0],
        **experiments[0].imageset.data().get_params(),
    ).get_proton_charge()

    tof_extract_shoeboxes_to_reflection_table(
        reflections,
        experiments[0],
        expt_data,
        incident_data,
        empty_data,
        expt_proton_charge,
        incident_proton_charge,
        empty_proton_charge,
        False,
    )

    ## Shoeboxes with incident/empty run normalisation
    ## and Lorentz correction
    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    tof_extract_shoeboxes_to_reflection_table(
        reflections,
        experiments[0],
        expt_data,
        incident_data,
        empty_data,
        expt_proton_charge,
        incident_proton_charge,
        empty_proton_charge,
        True,
    )

    ## Shoeboxes with incident/empty run normalisation
    ## and spherical absorption correction
    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    target_spectrum_sample_number_density = 0.0223
    target_spectrum_sample_radius = 0.3
    target_spectrum_scattering_x_section = 10.040
    target_spectrum_absorption_x_section = 17.015
    incident_spectrum_sample_number_density = 0.0722
    incident_spectrum_sample_radius = 0.3
    incident_spectrum_scattering_x_section = 5.158
    incident_spectrum_absorption_x_section = 4.4883

    corrections_data = TOFCorrectionsData(
        expt_proton_charge,
        incident_proton_charge,
        empty_proton_charge,
        target_spectrum_sample_radius,
        target_spectrum_scattering_x_section,
        target_spectrum_absorption_x_section,
        target_spectrum_sample_number_density,
        incident_spectrum_sample_radius,
        incident_spectrum_scattering_x_section,
        incident_spectrum_absorption_x_section,
        incident_spectrum_sample_number_density,
    )

    tof_extract_shoeboxes_to_reflection_table(
        reflections,
        experiments[0],
        expt_data,
        incident_data,
        empty_data,
        corrections_data,
        False,
    )

    ## Shoeboxes with incident/empty run normalisation,
    ## spherical absorption correction and Lorentz correction
    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    tof_extract_shoeboxes_to_reflection_table(
        reflections,
        experiments[0],
        expt_data,
        incident_data,
        empty_data,
        corrections_data,
        True,
    )
