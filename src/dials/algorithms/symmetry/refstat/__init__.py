# Python interface to refstat, (c) OVD, OlexSys, 2025
from __future__ import annotations

from dials_algorithms_symmetry_refstat_ext import extinctions_registry, merge_test


class registry(extinctions_registry):
    """Useful for testing"""

    def __init__(self):
        extinctions_registry.__init__(self)
        self.elements = {}
        for e in self:
            self.elements[e.name] = e

    def describe(self):
        """Prints information about the extinction elements and their relationship"""
        for k, e in self.elements.items():
            print(k)
            for rtm in e.rmx():
                print(f"\t{rtm}")
            sl = []
            for i in range(e.shadowed_by_size()):
                sl.append(e.get_shadowed_by(i).name)
            if len(sl) > 0:
                print("  Shadowed by: %s" % (" ".join(sl)))

    def show_extinctions_for(self, sg_name):
        sg = self.find_sg(sg_name)
        print(
            "Extinction elements for %s: %s"
            % (
                sg_name,
                " ".join([self.__getitem__(i).name for i in self.get_extinctions(sg)]),
            )
        )


class extinctions(extinctions_registry):
    def __init__(self, miller_array, sigma_level=5):
        """Sigma_level is used to identify systematic absences, some datasets may require
        lower values
        """
        extinctions_registry.__init__(self)
        self.all_elements = list([e for e in self])
        self.present, self.unique = [], []
        self.miller_array = miller_array
        self.sigma_level = sigma_level

    def analyse(self, scale_I_to=None):
        """Analyses the given miller array and collects statistics on present and
        unique extinction elements. Unique elements exclude 'shadowed' elements.
        Use scake_I_to to 'normalise' view between different datasets.
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
                continue
            x.meanI, x.sig = x.sumI / x.count, (x.sumS_sq**0.5) / x.count
            if x.meanI < self.sigma_level * x.sig:
                present.append(x)
        self.present = present
        unique = []
        for x in self.all_elements:
            if not x.count:
                continue
            if x in present and not x.is_shadowed_by(present):
                unique.append(x)
        self.unique = unique

    def print_stats(self):
        """Prints extinction elements statistics, must be called after 'analyse'"""
        for x in self.all_elements:
            if not x.count:
                continue
            if x in self.present:
                flag = "+"
                if x not in self.unique:
                    flag += "-"
            else:
                flag = "-"

            print(
                "%-4s (%5s): %16.2f(%6.2f) %s" % (x.name, x.count, x.meanI, x.sig, flag)
            )

    def get_all_matching_space_groups(self, centering="P"):
        """Returns a tuple(space_group, fraction_of_matching_elements). The list
        contains all of the space groups that match the given centring and unique
        extinctions elements, must be called after 'analyse'
        """
        sgs = []
        if not self.unique:
            return sgs
        u_s = set(x.id for x in self.unique)
        p_s = set(x.id for x in self.present)
        for i in range(self.sg_count()):
            sg = self.get_space_group(i)
            if not sg.name.startswith(centering):
                continue
            sge = set(self.get_extinctions(i))
            if len(sge & u_s) == len(u_s):
                sgs.append((sg, len(sge & p_s) / len(p_s)))
        return sgs

    def get_filtered_matching_space_groups(
        self, matches=None, cell_compatible_only=True
    ):
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
            weak_stats = mt.sysabs_test(sg, self.scale)
            wI = weak_stats.weak_I_sum / weak_stats.weak_count
            wIs = (weak_stats.weak_sig_sq_sum / weak_stats.weak_count) ** 0.5
            if wI > self.sigma_level * wIs:
                continue
            merge_stats = mt.merge_test(sg)
            sgs.append((merge_stats.r_int + wI / wIs / 10, sg))

        if len(sgs) == 0:
            return sgs
        if len(sgs) > 1:
            sgs.sort(key=lambda x: x[0])
            rv = [sgs[0][1]]
            r_int_th = sgs[0][0] * 2
            for i in range(1, len(sgs)):
                if sgs[i][0] > r_int_th:
                    break
                rv.append(sgs[i][1])
            return rv
        else:
            return [sgs[0][1]]
