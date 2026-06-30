"""Python interface to refstat, a tool to check space group symmetry based on
systematic absences in all space groups. This code is based on functionality
in Olex2 from OlexSys Ltd."""

from __future__ import annotations

import functools
import logging

import cctbx.miller
from cctbx import crystal

from dials.util import tabulate
from dials_algorithms_symmetry_refstat_ext import extinctions_registry, merge_test

logger = logging.getLogger(__name__)


class registry(extinctions_registry):
    """Useful for testing"""

    def __init__(self):
        extinctions_registry.__init__(self)
        self.elements = {}
        for e in self:
            self.elements[e.name] = e

    def describe(self):
        """Returns a string giving information about the extinction elements and their relationship"""
        lines = []
        for k, e in self.elements.items():
            lines.append(str(k))
            for rtm in e.rmx():
                lines.append(f"\t{rtm}")
            sl = []
            for i in range(e.shadowed_by_size()):
                sl.append(e.get_shadowed_by(i).name)
            if len(sl) > 0:
                lines.append("  Shadowed by: %s" % (" ".join(sl)))
        return "\n".join(lines)

    def show_extinctions_for(self, sg_name):
        """Returns a string giving the extinction elements for the given space group name"""
        sg = self.find_sg(sg_name)
        return "Extinction elements for %s: %s" % (
            sg_name,
            " ".join([self.__getitem__(i).name for i in self.get_extinctions(sg)]),
        )


class extinctions(extinctions_registry):
    def __init__(self, miller_array, sigma_level=5):
        """Sigma_level is used to identify systematic absences, some datasets may require
        lower values
        """
        extinctions_registry.__init__(self)
        self.all_elements = list(self)
        self.present, self.unique = [], []
        self.miller_array = miller_array
        self.sigma_level = sigma_level

    def analyse(self, scale_I_to=None):
        """Analyses the given miller array and collects statistics on present and
        unique extinction elements. Unique elements exclude 'shadowed' elements.
        Use scale_I_to to 'normalise' view between different datasets.
        """
        if scale_I_to is None:
            scale_I_to = 0.0
        if self.has_omp:
            self.process_omp(
                self.miller_array.indices(),
                self.miller_array.data(),
                self.miller_array.sigmas(),
                scale=scale_I_to,
            )
        else:
            self.process(
                self.miller_array.indices(),
                self.miller_array.data(),
                self.miller_array.sigmas(),
                scale=scale_I_to,
            )
        self.meanI = self.sumI / self.ref_count
        self.mean_sig = (self.sum_sig_sq / self.ref_count) ** 0.5

        present, unique = [], []
        for x in self.all_elements:
            if not x.count:
                x.meanI, x.sig = 0, 0
                present.append(x)
                continue
            x.meanI, x.sig = x.sumI / x.count, (x.sumS_sq**0.5) / x.count
            # check the overall intensity as well
            if (
                x.meanI < self.sigma_level * x.sig
                or x.meanI * self.sigma_level * 2 < self.meanI
            ):
                present.append(x)
        self.present = present
        unique = []
        for x in self.all_elements:
            if not x.count:
                continue
            if x in present and not x.is_shadowed_by(present):
                unique.append(x)
        self.unique = unique

    def get_symm_element_table(self):
        """Return the header and rows for the symmetry element table, must be called after 'analyse'"""
        headers = ["Element", "Count", "<I>", "<σ(I)>", "Present"]
        rows = []
        for x in self.all_elements:
            if x in self.present:
                flag = "+"
                if x not in self.unique:
                    flag += "-"
            else:
                flag = "-"

            rows.append([x.name, f"{x.count}", f"{x.meanI:.2f}", f"{x.sig:.2f}", flag])
        return rows, headers

    def show_stats(self):
        """Return a string giving extinction elements statistics, must be called after 'analyse'"""
        lines = []
        for x in self.all_elements:
            if x in self.present:
                flag = "+"
                if x not in self.unique:
                    flag += "-"
            else:
                flag = "-"

            lines.append(
                "%-4s (%5s): %16.2f(%6.2f) %s" % (x.name, x.count, x.meanI, x.sig, flag)
            )
        return "\n".join(lines)

    def get_all_matching_space_groups(self, centering="P"):
        """Returns a tuple(space_group, fraction_of_matching_elements). The list
        contains all of the space groups that match the given centring and unique
        extinctions elements, must be called after 'analyse'
        """
        sgs, all_sg = [], []
        u_s = {x.id for x in self.unique}
        p_s = {x.id for x in self.present}
        if not p_s:
            return sgs
        for i in range(self.sg_count()):
            sg = self.get_space_group(i)
            if not sg.name.startswith(centering):
                continue
            sge = set(self.get_extinctions(i))
            matches_n = len(sge & u_s)
            if matches_n == len(u_s):
                sgs.append((sg, len(sge & p_s) / len(p_s)))
            elif matches_n:
                all_sg.append((sg, len(sge & p_s) / len(p_s)))
        if not sgs and u_s:
            all_sg.sort(key=lambda x: x[1])
            sgs = all_sg
        return sgs

    def get_filtered_matching_space_groups(
        self, matches=None, cell_compatible_only=True
    ):
        """
        Filter the list of matching space groups, selecting only those that are
        compatible with the unit cell of the miller array (if cell_compatible_only
        is True) and have acceptable systematic absence statistics. The filtered
        list is then sorted based on Rint, where for elements with similar Rint
        the number of matching present elements is used to break ties. For
        elements with the same number of matching present elements, the number
        of weak reflections is used to break ties.

        Parameters:
        matches (list): The list of matching space groups to filter.
        cell_compatible_only (bool): If true, only consider space groups that are
            compatible with the unit cell of the given miller array.

        Returns:
        list: A list of matching space groups that pass the filter criteria.
        """
        if matches is None:
            matches = self.get_all_matching_space_groups()
        sgs = []
        if cell_compatible_only:
            uc = self.miller_array.crystal_symmetry().unit_cell()
        for sg, mp in matches:
            if cell_compatible_only and not sg.is_compatible_unit_cell(uc):
                continue
            mt = merge_test(
                self.miller_array.indices(),
                self.miller_array.data(),
                self.miller_array.sigmas(),
            )
            sysabs_stats = mt.sysabs_test(sg, self.scale)
            if sysabs_stats.weak_count:
                wI = sysabs_stats.weak_I_sum / sysabs_stats.weak_count
                wIs = (sysabs_stats.weak_sig_sq_sum / sysabs_stats.weak_count) ** 0.5
                if wI > self.sigma_level * wIs:
                    continue
                if (
                    sysabs_stats.strong_count
                    and wI > sysabs_stats.strong_I_sum / sysabs_stats.strong_count
                ):
                    continue
            merge_stats = mt.merge_test(sg)
            sgs.append(
                (
                    merge_stats.r_int,
                    sg,
                    mp,
                    sysabs_stats.weak_count,
                    merge_stats.inconsistent_count,
                )
            )

        if len(sgs) == 0:
            return sgs
        if len(sgs) > 1:

            def cmp(a, b):
                if abs(a[0] - b[0]) < 0.5e-2:  # check for Rint
                    if a[2] == b[2]:  # check for number of matches
                        return b[3] - a[3]  # check the number of weak refs
                    return b[2] - a[2]
                return a[0] - b[0]

            sgs.sort(key=functools.cmp_to_key(cmp))
            rv = [sgs[0][1]]
            r_int_th = sgs[0][0] * 1.5
            i_eq_th = sgs[0][4]
            for i in range(1, len(sgs)):
                if sgs[i][0] > r_int_th:
                    continue
                if sgs[i][4] > i_eq_th:
                    continue
                rv.append(sgs[i][1])
            return rv
        else:
            return [sgs[0][1]]


def check_reflections(experiments, reflections, params):
    cs = crystal.symmetry(experiments[0].crystal.get_unit_cell(), "P1")
    centering = (
        experiments[0]
        .crystal.get_space_group()
        .match_tabulated_settings()
        .hermann_mauguin()[0]
    )
    miller_set = cctbx.miller.set(
        crystal_symmetry=cs,
        indices=reflections["miller_index"],
        anomalous_flag=False,
    )
    i_obs = cctbx.miller.array(miller_set, data=reflections["intensity"])
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(cctbx.array_family.flex.sqrt(reflections["variance"]))
    i_obs.set_info(
        cctbx.miller.array_info(source="DIALS", source_type="reflection_tables")
    )

    xr = registry()

    miller_array = i_obs.merge_equivalents(algorithm="gaussian").array()
    data = miller_array.data()
    sigmas = miller_array.sigmas()
    logger.info("Uniq in P 1: %s" % (len(data)))
    xr.process(miller_array.indices(), data, sigmas)
    xr.reset()

    sa = extinctions(miller_array, sigma_level=params.systematic_absences.sigma_level)
    sa.analyse(scale_I_to=1)
    rows, headers = sa.get_symm_element_table()
    logger.info(tabulate(rows, headers))

    logger.info(
        "<I>: %.3f and <σ(I)>: %.2f for %s unique reflections"
        % (sa.meanI, sa.mean_sig, sa.ref_count)
    )
    matches = sa.get_all_matching_space_groups(centering=centering)
    filtered_matches = sa.get_filtered_matching_space_groups(matches=matches)

    # merge_test object
    t = merge_test(miller_array.indices(), data, sigmas)
    rows = []
    matches = {sg.name: mp for sg, mp in matches}
    # Loop through the filtered matches
    for sg in filtered_matches:
        mp = matches[sg.name]
        weak_stats = t.sysabs_test(sg, sa.scale)
        if weak_stats.weak_count:
            wI = weak_stats.weak_I_sum / weak_stats.weak_count
            wIs = (weak_stats.weak_sig_sq_sum / weak_stats.weak_count) ** 0.5
            if wI > 5 * wIs:
                continue
        else:
            wI, wIs = 0, 0
        if weak_stats.strong_count:
            sI = weak_stats.strong_I_sum / weak_stats.strong_count
            sIs = (weak_stats.strong_sig_sq_sum / weak_stats.strong_count) ** 0.5
        else:
            sI, sIs = 0, 0
        merge_stats = t.merge_test(sg)
        rows.append(
            [
                sg.name,
                "✔" if sg.is_centric() else "",
                f"{int(mp * 100)}",
                merge_stats.inconsistent_count,
                f"{merge_stats.r_int * 100:.3f}",
                weak_stats.weak_count,
                f"{wI / wIs if wIs else 0:.2f}",
                weak_stats.strong_count,
                f"{sI / sIs if sIs else 0:.2f}",
            ]
        )
    if rows:
        logger.info(
            tabulate(
                rows,
                headers=[
                    "Space group",
                    "Centric",
                    "Matches (%)",
                    "Incons.\nequiv.",
                    "Rint",
                    "#Weak",
                    "Weak I/σ(I)",
                    "#Strong",
                    "Strong I/σ(I)",
                ],
            )
        )

    # Further filtering to select only those SGs in the same Laue group
    space_group = experiments[0].crystal.get_space_group()
    laue_group = str(space_group.build_derived_patterson_group().info())
    filtered_matches = [
        sg
        for sg in filtered_matches
        if str(sg.build_derived_patterson_group().info()) == laue_group
    ]

    # Set the space group for all experiments to the first match
    if filtered_matches:
        logger.info(
            f"Selected results in Laue group {laue_group}: "
            + ", ".join([sg.name for sg in filtered_matches])
        )
        for experiment in experiments:
            experiment.crystal.set_space_group(filtered_matches[0])
    else:
        logger.info(f"No space groups found in Laue group {laue_group}")
