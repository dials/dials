#
# report.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division, print_function

import collections
from dials.array_family import flex
from dials.array_family.flex import Binner
from dials.util.report import Array, Report, Table

import six


def flex_ios(val, var):
    """
    Compute I/sigma or return zero for each element.

    """
    assert len(val) == len(var)
    result = flex.double(len(val), 0)
    indices = flex.size_t(range(len(val))).select(var > 0)
    val = val.select(indices)
    var = var.select(indices)
    assert var.all_gt(0)
    result.set_selected(indices, val / flex.sqrt(var))
    return result


def generate_integration_report(experiment, reflections, n_resolution_bins=20):
    """
    Generate the integration report

    """
    from dials.algorithms.statistics import pearson_correlation_coefficient
    from dials.algorithms.statistics import spearman_correlation_coefficient
    from cctbx import miller, crystal

    def overall_report(data):

        # Start by adding some overall numbers
        report = collections.OrderedDict()
        report["n"] = len(reflections)
        report["n_full"] = data["full"].count(True)
        report["n_partial"] = data["full"].count(False)
        report["n_overload"] = data["over"].count(True)
        report["n_ice"] = data["ice"].count(True)
        report["n_summed"] = data["sum"].count(True)
        report["n_fitted"] = data["prf"].count(True)
        report["n_integated"] = data["int"].count(True)
        report["n_invalid_bg"] = data["ninvbg"].count(True)
        report["n_invalid_fg"] = data["ninvfg"].count(True)
        report["n_failed_background"] = data["fbgd"].count(True)
        report["n_failed_summation"] = data["fsum"].count(True)
        report["n_failed_fitting"] = data["fprf"].count(True)

        # Compute mean background
        try:
            report["mean_background"] = flex.mean(
                data["background.mean"].select(data["int"])
            )
        except Exception:
            report["mean_background"] = 0.0

        # Compute mean I/Sigma summation
        try:
            report["ios_sum"] = flex.mean(data["intensity.sum.ios"].select(data["sum"]))
        except Exception:
            report["ios_sum"] = 0.0

        # Compute mean I/Sigma profile fitting
        try:
            report["ios_prf"] = flex.mean(data["intensity.prf.ios"].select(data["prf"]))
        except Exception:
            report["ios_prf"] = 0.0

        # Compute the mean profile correlation
        try:
            report["cc_prf"] = flex.mean(
                data["profile.correlation"].select(data["prf"])
            )
        except Exception:
            report["cc_prf"] = 0.0

        # Compute the correlations between summation and profile fitting
        try:
            mask = data["sum"] & data["prf"]
            Isum = data["intensity.sum.value"].select(mask)
            Iprf = data["intensity.prf.value"].select(mask)
            report["cc_pearson_sum_prf"] = pearson_correlation_coefficient(Isum, Iprf)
            report["cc_spearman_sum_prf"] = spearman_correlation_coefficient(Isum, Iprf)
        except Exception:
            report["cc_pearson_sum_prf"] = 0.0
            report["cc_spearman_sum_prf"] = 0.0

        # Return the overall report
        return report

    def binned_report(binner, index, data):

        # Create the indexers
        indexer_all = binner.indexer(index)
        indexer_sum = binner.indexer(index.select(data["sum"]))
        indexer_prf = binner.indexer(index.select(data["prf"]))
        indexer_int = binner.indexer(index.select(data["int"]))

        # Add some stats by resolution
        report = collections.OrderedDict()
        report["bins"] = list(binner.bins())
        report["n_full"] = list(indexer_all.sum(data["full"]))
        report["n_partial"] = list(indexer_all.sum(~data["full"]))
        report["n_overload"] = list(indexer_all.sum(data["over"]))
        report["n_ice"] = list(indexer_all.sum(data["ice"]))
        report["n_summed"] = list(indexer_all.sum(data["sum"]))
        report["n_fitted"] = list(indexer_all.sum(data["prf"]))
        report["n_integrated"] = list(indexer_all.sum(data["int"]))
        report["n_invalid_bg"] = list(indexer_all.sum(data["ninvbg"]))
        report["n_invalid_fg"] = list(indexer_all.sum(data["ninvfg"]))
        report["n_failed_background"] = list(indexer_all.sum(data["fbgd"]))
        report["n_failed_summation"] = list(indexer_all.sum(data["fsum"]))
        report["n_failed_fitting"] = list(indexer_all.sum(data["fprf"]))

        # Compute mean background
        try:
            report["mean_background"] = list(
                indexer_int.mean(data["background.mean"].select(data["int"]))
            )
        except Exception:
            report["mean_background"] = [0.0] * len(binner)

        # Compute mean I/Sigma summation
        try:
            report["ios_sum"] = list(
                indexer_sum.mean(data["intensity.sum.ios"].select(data["sum"]))
            )
        except Exception:
            report["ios_sum"] = [0.0] * len(binner)

        # Compute mean I/Sigma profile fitting
        try:
            report["ios_prf"] = list(
                indexer_prf.mean(data["intensity.prf.ios"].select(data["prf"]))
            )
        except Exception:
            report["ios_prf"] = [0.0] * len(binner)

        # Compute the mean profile correlation
        try:
            report["cc_prf"] = list(
                indexer_prf.mean(data["profile.correlation"].select(data["prf"]))
            )
        except Exception:
            report["cc_prf"] = [0.0] * len(binner)

        try:
            report["rmsd_xy"] = list(
                indexer_sum.mean(data["xyz.rmsd"].select(data["int"]))
            )
        except Exception:
            report["rmsd_xy"] = [0.0] * len(binner)

        # Return the binned report
        return report

    def resolution_bins(experiment, hkl, nbins):

        # Create the crystal symmetry object
        cs = crystal.symmetry(
            space_group=experiment.crystal.get_space_group(),
            unit_cell=experiment.crystal.get_unit_cell(),
        )

        # Create the resolution binner object
        ms = miller.set(cs, hkl)
        ms.setup_binner(n_bins=nbins)
        binner = ms.binner()
        brange = list(binner.range_used())
        bins = [binner.bin_d_range(brange[0])[0]]
        for i in brange:
            bins.append(binner.bin_d_range(i)[1])
        return flex.double(reversed(bins))

    def select(data, indices):
        # Select rows from columns
        result = {key: value.select(indices) for key, value in six.iteritems(data)}
        return result

    # Check the required columns are there
    assert "miller_index" in reflections
    assert "d" in reflections
    assert "flags" in reflections
    assert "bbox" in reflections
    assert "xyzcal.px" in reflections
    assert "partiality" in reflections
    assert "intensity.sum.value" in reflections
    assert "intensity.sum.variance" in reflections

    # Get the flag enumeration
    flags = flex.reflection_table.flags

    # Get some keys from the data
    data = {}
    for key in [
        "miller_index",
        "xyzcal.px",
        "xyzobs.px.value",
        "d",
        "bbox",
        "background.mean",
        "partiality",
        "intensity.sum.value",
        "intensity.sum.variance",
        "intensity.prf.value",
        "intensity.prf.variance",
        "profile.correlation",
    ]:
        if key in reflections:
            data[key] = reflections[key]

    # Compute some flag stuff
    data["full"] = data["partiality"] > 0.997300203937
    data["over"] = reflections.get_flags(flags.overloaded)
    data["ice"] = reflections.get_flags(flags.in_powder_ring)
    data["sum"] = reflections.get_flags(flags.integrated_sum)
    data["prf"] = reflections.get_flags(flags.integrated_prf)
    data["int"] = reflections.get_flags(flags.integrated, all=False)
    data["ninvbg"] = reflections.get_flags(flags.background_includes_bad_pixels)
    data["ninvfg"] = reflections.get_flags(flags.foreground_includes_bad_pixels)
    data["fbgd"] = reflections.get_flags(flags.failed_during_background_modelling)
    data["fsum"] = reflections.get_flags(flags.failed_during_summation)
    data["fprf"] = reflections.get_flags(flags.failed_during_profile_fitting)

    # Try to calculate the i over sigma for summation
    data["intensity.sum.ios"] = flex_ios(
        data["intensity.sum.value"], data["intensity.sum.variance"]
    )

    # Try to calculate the i over sigma for profile fitting
    try:
        data["intensity.prf.ios"] = flex_ios(
            data["intensity.prf.value"], data["intensity.prf.variance"]
        )
    except Exception:
        pass

    # Try to calculate the rmsd between observation and prediction
    try:
        xcal, ycal, zcal = data["xyzcal.px"].parts()
        xobs, yobs, zobs = data["xyzobs.px.value"].parts()
        data["xyz.rmsd"] = flex.sqrt(flex.pow2(xcal - xobs) + flex.pow2(ycal - yobs))
    except Exception:
        pass

    # Create the resolution binner
    resolution_binner = Binner(
        resolution_bins(experiment, data["miller_index"], n_resolution_bins)
    )

    # Create the frame binner object
    try:
        array_range = experiment.imageset.get_array_range()
    except Exception:
        array_range = (0, len(experiment.imageset))
    frame_binner = Binner(
        flex.int(range(array_range[0], array_range[1] + 1)).as_double()
    )

    # Create the overall report
    overall = overall_report(data)

    # Create high/low resolution reports
    hl_binner = resolution_binner.indexer(data["d"])
    high_summary = overall_report(select(data, hl_binner.indices(0)))
    low_summary = overall_report(select(data, hl_binner.indices(n_resolution_bins - 1)))
    high_summary["dmin"] = resolution_binner.bins()[0]
    high_summary["dmax"] = resolution_binner.bins()[1]
    low_summary["dmin"] = resolution_binner.bins()[n_resolution_bins - 1]
    low_summary["dmax"] = resolution_binner.bins()[n_resolution_bins]
    overall["dmin"] = high_summary["dmin"]
    overall["dmax"] = low_summary["dmax"]

    # Create the overall report
    summary = collections.OrderedDict(
        [("overall", overall), ("low", low_summary), ("high", high_summary)]
    )

    # Create a report binned by resolution
    resolution = binned_report(resolution_binner, data["d"], data)

    # Create the report binned by image
    image = binned_report(frame_binner, data["xyzcal.px"].parts()[2], data)

    # Return the report
    return collections.OrderedDict(
        [("summary", summary), ("resolution", resolution), ("image", image)]
    )


class IntegrationReport(Report):
    """
    A class to store the integration report

    """

    def __init__(self, experiments, reflections):
        """
        Create the integration report

        :param experiments: The experiment list
        :param reflections: The reflection table

        """
        # Initialise the report class
        super(IntegrationReport, self).__init__()

        # Split the tables by experiment id
        tables = reflections.split_by_experiment_id()
        assert len(tables) == len(experiments)

        # Initialise the dictionary
        report_list = []

        # Generate an integration report for each experiment
        for i, (expr, data) in enumerate(zip(experiments, tables)):
            report_list.append(generate_integration_report(expr, data))

        # Construct the per image table
        table = Table()
        table.name = "integration.image.summary"
        table.title = "Summary vs image number"
        table.cols.append(("id", "ID"))
        table.cols.append(("image", "Image"))
        table.cols.append(("n_full", "# full"))
        table.cols.append(("n_part", "# part"))
        table.cols.append(("n_over", "# over"))
        table.cols.append(("n_ice", "# ice"))
        table.cols.append(("n_sum", "# sum"))
        table.cols.append(("n_prf", "# prf"))
        table.cols.append(("ibg", "Ibg"))
        table.cols.append(("ios_sum", "I/sigI\n (sum)"))
        table.cols.append(("ios_prf", "I/sigI\n (prf)"))
        table.cols.append(("cc_prf", "CC prf"))
        table.cols.append(("rmsd_xy", "RMSD XY"))
        for j, report in enumerate(report_list):
            report = report["image"]
            for i in range(len(report["bins"]) - 1):
                table.rows.append(
                    [
                        "%d" % j,
                        "%d" % report["bins"][i],
                        "%d" % report["n_full"][i],
                        "%d" % report["n_partial"][i],
                        "%d" % report["n_overload"][i],
                        "%d" % report["n_ice"][i],
                        "%d" % report["n_summed"][i],
                        "%d" % report["n_fitted"][i],
                        "%.2f" % report["mean_background"][i],
                        "%.2f" % report["ios_sum"][i],
                        "%.2f" % report["ios_prf"][i],
                        "%.2f" % report["cc_prf"][i],
                        "%.2f" % report["rmsd_xy"][i],
                    ]
                )
        self.add_table(table)

        # Construct the per resolution table
        table = Table()
        table.name = "integration.resolution.summary"
        table.title = "Summary vs resolution"
        table.cols.append(("id", "ID"))
        table.cols.append(("dmin", "d min"))
        table.cols.append(("n_full", "# full"))
        table.cols.append(("n_part", "# part"))
        table.cols.append(("n_over", "# over"))
        table.cols.append(("n_ice", "# ice"))
        table.cols.append(("n_sum", "# sum"))
        table.cols.append(("n_prf", "# prf"))
        table.cols.append(("ibg", "Ibg"))
        table.cols.append(("ios_sum", "I/sigI\n (sum)"))
        table.cols.append(("ios_prf", "I/sigI\n (prf)"))
        table.cols.append(("cc_prf", "CC prf"))
        table.cols.append(("rmsd_xy", "RMSD XY"))
        for j, report in enumerate(report_list):
            report = report["resolution"]
            for i in range(len(report["bins"]) - 1):
                table.rows.append(
                    [
                        "%d" % j,
                        "%.2f" % report["bins"][i],
                        "%d" % report["n_full"][i],
                        "%d" % report["n_partial"][i],
                        "%d" % report["n_overload"][i],
                        "%d" % report["n_ice"][i],
                        "%d" % report["n_summed"][i],
                        "%d" % report["n_fitted"][i],
                        "%.2f" % report["mean_background"][i],
                        "%.2f" % report["ios_sum"][i],
                        "%.2f" % report["ios_prf"][i],
                        "%.2f" % report["cc_prf"][i],
                        "%.2f" % report["rmsd_xy"][i],
                    ]
                )
        self.add_table(table)

        # Create the overall table
        for j, report in enumerate(report_list):
            report = report["summary"]
            summary = report["overall"]
            high = report["high"]
            low = report["low"]

            table = Table()
            table.name = "integration.overall.summary"
            table.title = "Summary for experiment %d" % j
            table.cols.append(("item", "Item"))
            table.cols.append(("overall", "Overall"))
            table.cols.append(("low", "Low"))
            table.cols.append(("high", "High"))
            desc_fmt_key = [
                ("dmin", "%.2f", "dmin"),
                ("dmax", "%.2f", "dmax"),
                ("number fully recorded", "%d", "n_full"),
                ("number partially recorded", "%d", "n_partial"),
                ("number with invalid background pixels", "%d", "n_invalid_bg"),
                ("number with invalid foreground pixels", "%d", "n_invalid_fg"),
                ("number with overloaded pixels", "%d", "n_overload"),
                ("number in powder rings", "%d", "n_ice"),
                ("number processed with summation", "%d", "n_summed"),
                ("number processed with profile fitting", "%d", "n_fitted"),
                ("number failed in background modelling", "%d", "n_failed_background"),
                ("number failed in summation", "%d", "n_failed_summation"),
                ("number failed in profile fitting", "%d", "n_failed_fitting"),
                ("ibg", "%.2f", "mean_background"),
                ("i/sigi (summation)", "%.2f", "ios_sum"),
                ("i/sigi (profile fitting)", "%.2f", "ios_prf"),
                ("cc prf", "%.2f", "cc_prf"),
                ("cc_pearson sum/prf", "%.2f", "cc_pearson_sum_prf"),
                ("cc_spearman sum/prf", "%.2f", "cc_spearman_sum_prf"),
            ]
            for desc, fmt, key in desc_fmt_key:
                table.rows.append(
                    [desc, fmt % summary[key], fmt % low[key], fmt % high[key]]
                )
            self.add_table(table)


class ProfileModelReport(Report):
    """
    A class to store the profile model report

    """

    def __init__(self, experiments, fitter, reflections):
        """
        Create the integration report

        :param experiments: The experiment list
        :param profile_model: The profile model
        :param reflections: The reflection table

        """
        # Initialise the report class
        super(ProfileModelReport, self).__init__()

        # Create the table
        table = Table()

        # Set the title
        table.name = "profile.summary"
        table.title = "Summary of profile model"

        # Add the columns
        table.cols.append(("id", "ID"))
        table.cols.append(("profile", "Profile"))
        table.cols.append(("created", "Created"))
        table.cols.append(("x", "X (px)"))
        table.cols.append(("y", "Y (px)"))
        table.cols.append(("z", "Z (im)"))
        table.cols.append(("n_reflections", "# reflections"))

        # Create the summary for each profile model
        for i in range(len(fitter)):
            model = fitter[i]
            for j in range(len(model)):
                table.rows.append(
                    [
                        "%d" % i,
                        "%d" % j,
                        "%s" % model.valid(j),
                        "%.2f" % model.coord(j)[0],
                        "%.2f" % model.coord(j)[1],
                        "%.2f" % model.coord(j)[2],
                        "%d" % model.n_reflections(j),
                    ]
                )

        # Add the table
        self.add_table(table)

        # Add the profiles
        for i in range(len(fitter)):
            model = fitter[i]
            for j in range(len(model)):
                if model.valid(j):
                    array = Array()
                    array.name = "profile.model.%d.%d" % (i, j)
                    array.title = "Profile model (id: %d, profile: %d)" % (i, j)
                    array.data = model.data(j)
                    self.add_array(array)


class ProfileModelReport2(Report):
    """
    A class to store the profile model report

    """

    def __init__(self, experiments, reference, reflections):
        """
        Create the integration report

        :param experiments: The experiment list
        :param reference: The reference profiles
        :param reflections: The reflection table

        """
        # Initialise the report class
        super(ProfileModelReport2, self).__init__()

        # Create the table
        table = Table()

        # Set the title
        table.name = "profile.summary"
        table.title = "Summary of profile model"

        # Add the columns
        table.cols.append(("id", "ID"))
        table.cols.append(("profile", "Profile"))
        table.cols.append(("created", "Created"))
        table.cols.append(("x", "X (px)"))
        table.cols.append(("y", "Y (px)"))
        table.cols.append(("z", "Z (im)"))
        table.cols.append(("n_reflections", "# reflections"))

        # Create the summary for each profile model
        for i in range(len(reference)):
            model = reference[i]
            for j in range(len(model)):
                table.rows.append(
                    [
                        "%d" % i,
                        "%d" % j,
                        "%s" % True,
                        "%.2f" % model.sampler().coord(j)[0],
                        "%.2f" % model.sampler().coord(j)[1],
                        "%.2f" % model.sampler().coord(j)[2],
                        "%d" % model.n_reflections(j),
                    ]
                )

        # Add the table
        self.add_table(table)

        # Add the profiles
        for i in range(len(fitter)):
            model = fitter[i]
            for j in range(len(model)):
                if model.valid(j):
                    array = Array()
                    array.name = "profile.model.%d.%d" % (i, j)
                    array.title = "Profile model (id: %d, profile: %d)" % (i, j)
                    array.data = model.data(j)
                    self.add_array(array)


class ProfileValidationReport(Report):
    """
    A class to store the profile validation report

    """

    def __init__(self, experiments, profile_fitter, reflections, num_folds):
        """
        Create the integration report

        :param experiments: The experiment list
        :param profile_model: The profile model
        :param reflections: The reflection table

        """
        # Initialise the report class
        super(ProfileValidationReport, self).__init__()

        # Create the table
        table = Table()

        # Set the title
        table.name = "validation.summary"
        table.title = "Summary of profile validation "

        # Add the columns
        table.cols.append(("id", "ID"))
        table.cols.append(("subsample", "Sub-sample"))
        table.cols.append(("n_valid", "# validated"))
        table.cols.append(("cc", "<CC>"))
        table.cols.append(("nrmsd", "<NRMSD>"))

        # Split the reflections
        reflection_tables = reflections.split_by_experiment_id()
        assert len(reflection_tables) == len(experiments)
        assert len(profile_fitter) == num_folds

        # Create the summary for each profile model
        for i in range(len(reflection_tables)):
            reflection_table = reflection_tables[i]
            reflection_table = reflection_table.select(
                reflection_table.get_flags(reflection_table.flags.integrated_prf)
            )
            index = reflection_table["profile.index"]
            cc = reflection_table["profile.correlation"]
            nrmsd = reflection_table["profile.rmsd"]
            for j in range(num_folds):
                mask = index == j
                num_validated = mask.count(True)
                if num_validated == 0:
                    mean_cc = 0
                    mean_nrmsd = 0
                else:
                    mean_cc = flex.mean(cc.select(mask))
                    mean_nrmsd = flex.mean(nrmsd.select(mask))
                table.rows.append(
                    [
                        "%d" % i,
                        "%d" % j,
                        "%d" % num_validated,
                        "%.2f" % mean_cc,
                        "%.2f" % mean_nrmsd,
                    ]
                )

        # Add the table
        self.add_table(table)
