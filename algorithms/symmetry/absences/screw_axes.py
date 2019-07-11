"""Definitions of screw axes with methods for scoring against data."""
from __future__ import absolute_import, division, print_function
import math
import logging
from scitbx.array_family import flex
from dials.algorithms.symmetry.absences.plots import plot_screw_axes
from dials.util.observer import Observer, Subject, singleton
from jinja2 import Environment, ChoiceLoader, PackageLoader


logger = logging.getLogger("dials.space_group")


@singleton
class ScrewAxisObserver(Observer):

    """Observer to record data used in screw axis analysis."""

    def update(self, screw_axis):
        self.data[screw_axis.name] = {
            "miller_axis_vals": screw_axis.miller_axis_vals,
            "i_over_sigma": screw_axis.i_over_sigma,
            "intensities": screw_axis.intensities,
            "sigmas": screw_axis.sigmas,
        }

    def generate_html_report(self, filename):
        """Generate a html report using the data."""
        screw_axes_graphs = plot_screw_axes(self.data)
        self.data["screw_axes"] = screw_axes_graphs
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("systematic_absences_report.html")
        html = template.render(
            page_title="DIALS systematic absences report",
            screw_axes_graphs=self.data["screw_axes"],
        )
        with open(filename, "wb") as f:
            f.write(html.encode("ascii", "xmlcharrefreplace"))


class ScrewAxis(Subject):

    """Definition of a generic screw axis."""

    axis_idx = None  # x=0, y=1, z=2
    axis_repeat = None  # repeat of present reflections e.g =4 for 41, =2 for 42
    name = None

    def __init__(self):
        super(ScrewAxis, self).__init__(events=["selected data for scoring"])
        self.equivalent_axes = []
        self.n_refl_used = (0.0, 0.0)
        self.miller_axis_vals = []
        self.i_over_sigma = []
        self.intensities = []
        self.sigmas = []
        self.mean_I_sigma_abs = 0.0
        self.mean_I_sigma = 0.0
        self.mean_I_abs = 0.0
        self.mean_I = 0.0

    def add_equivalent_axis(self, equivalent):
        """Add a symmetry equivalent axis."""
        self.equivalent_axes.append(equivalent)

    def select_axial_reflections(self, miller_indices):
        """Select reflections along the screw axis."""
        h, k, l = miller_indices.as_vec3_double().parts()
        if self.axis_idx == 0:
            selection = (k == 0) & (l == 0)
        elif self.axis_idx == 1:
            selection = (h == 0) & (l == 0)
        else:
            selection = (h == 0) & (k == 0)
        return selection

    @Subject.notify_event(event="selected data for scoring")
    def get_all_suitable_reflections(self, reflection_table):
        """Select suitable reflections for testing the screw axis."""
        refl = reflection_table
        sel = self.select_axial_reflections(refl["miller_index"])
        miller_idx = refl["miller_index"].select(sel)
        self.miller_axis_vals = miller_idx.as_vec3_double().parts()[self.axis_idx]
        self.intensities = refl["intensity"].select(sel)
        self.sigmas = refl["variance"].select(sel) ** 0.5
        self.i_over_sigma = self.intensities / self.sigmas

        if self.equivalent_axes:
            for a in self.equivalent_axes:
                sel = a.select_axial_reflections(refl["miller_index"])
                miller_idx = refl["miller_index"].select(sel)
                intensities = refl["intensity"].select(sel)
                sigmas = refl["variance"].select(sel) ** 0.5
                self.i_over_sigma.extend(intensities / sigmas)
                self.miller_axis_vals.extend(
                    miller_idx.as_vec3_double().parts()[a.axis_idx]
                )
                self.intensities.extend(intensities)
                self.sigmas.extend(sigmas)

    def score_axis(self, reflection_table, significance_level=0.95):
        """Score the axis give a reflection table of data."""
        assert significance_level in [0.95, 0.975, 0.99]
        self.get_all_suitable_reflections(reflection_table)

        expected_sel = self.miller_axis_vals.iround() % self.axis_repeat == 0

        expected = self.i_over_sigma.select(expected_sel)
        expected_abs = self.i_over_sigma.select(~expected_sel)
        self.n_refl_used = (expected.size(), expected_abs.size())
        # Limit to best #n reflections to avoid weak at high res - use wilson B?
        if not expected or not expected_abs:
            return 0.0

        # z = (sample mean - population mean) / standard error
        S_E_abs = 1.0  # errors probably correlated so say standard error = 1
        S_E_pres = 1.0  # / expected.size() ** 0.5

        self.mean_I_sigma_abs = flex.mean(expected_abs)
        self.mean_I_sigma = flex.mean(expected)

        self.mean_I = flex.mean(self.intensities.select(expected_sel))
        self.mean_I_abs = flex.mean(self.intensities.select(~expected_sel))

        z_score_absent = self.mean_I_sigma_abs / S_E_abs
        z_score_present = self.mean_I_sigma / S_E_pres

        # get a p-value for z > z_score
        P_absent = 0.5 * (1.0 + math.erf(z_score_absent / (2 ** 0.5)))
        P_present = 0.5 * (1.0 + math.erf(z_score_present / (2 ** 0.5)))

        # sanity check - is most of intensity in 'expected' channel?
        intensity_test = self.mean_I_sigma > (20.0 * self.mean_I_sigma_abs)

        cutoffs = {0.95: 1.645, 0.975: 1.960, 0.99: 2.326}
        cutoff = cutoffs[significance_level]

        if z_score_absent > cutoff and not intensity_test:
            # z > 1.65 in only 5% of cases for normal dist
            # significant nonzero intensity where expected absent.
            return (1.0 - P_absent) * P_present
        elif z_score_absent > cutoff:
            # results appear inconsistent - significant i_over_sigma_abs, but this
            # is still low compared to i_over_sigma_expected
            # try removing the highest absent reflection in case its an outlier
            sel = flex.sort_permutation(expected_abs)
            sorted_exp_abs = expected_abs.select(sel)
            mean_i_sigma_abs = flex.mean(sorted_exp_abs[:-1])
            if (mean_i_sigma_abs / S_E_abs) > cutoff:
                # Still looks like reflections in expected absent
                logger.info(
                    """Test results for %s appear inconsistent (significant nonzero intensity for
'absent' reflections, but majority of intensity in reflection condition).
There may be a set of weak reflections due to pseudosymmetry.""",
                    self.name,
                )
                # Still high intensity of absent, so return as before
                return (1.0 - P_absent) * P_present
            # Looks like there was an outlier, now 'absent' reflections ~ 0.
            self.mean_I_sigma_abs = mean_i_sigma_abs
            self.mean_I_abs = flex.mean(
                self.intensities.select(~expected_sel).select(sel)[:-1]
            )
            if z_score_present > cutoff:
                logger.info(
                    """Screw axis %s could only be assigned after removing a suspected outlier
from the expected 'absent' reflections.""",
                    self.name,
                )
                return P_present
            # else in the uncertain case
        elif z_score_present > cutoff:  # evidence with confidence
            return P_present
        logger.info(
            """No evidence to suggest a screw axis for %s, but insufficient
evidence to rule out completely, possibly due to limited data.""",
            self.name,
        )
        return 0.0  # should this be zero or a small number?


class ScrewAxis21c(ScrewAxis):

    """Definition of a 21c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "21c"


class ScrewAxis21b(ScrewAxis):

    """Definition of a 21b screw axis"""

    axis_idx = 1
    axis_repeat = 2
    name = "21b"


class ScrewAxis21a(ScrewAxis):

    """Definition of a 21a screw axis"""

    axis_idx = 0
    axis_repeat = 2
    name = "21a"


class ScrewAxis41c(ScrewAxis):

    """Definition of a 41c screw axis"""

    axis_idx = 2
    axis_repeat = 4
    name = "41c"


class ScrewAxis42c(ScrewAxis):

    """Definition of a 42c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "42c"


class ScrewAxis41b(ScrewAxis):

    """Definition of a 41b screw axis"""

    axis_idx = 1
    axis_repeat = 4
    name = "41b"


class ScrewAxis41a(ScrewAxis):

    """Definition of a 41a screw axis"""

    axis_idx = 0
    axis_repeat = 4
    name = "41a"


class ScrewAxis31c(ScrewAxis):

    """Definition of a 31c screw axis"""

    axis_idx = 2
    axis_repeat = 3
    name = "31c"


class ScrewAxis61c(ScrewAxis):

    """Definition of a 61c screw axis"""

    axis_idx = 2
    axis_repeat = 6
    name = "61c"


class ScrewAxis62c(ScrewAxis):

    """Definition of a 62c screw axis"""

    axis_idx = 2
    axis_repeat = 3
    name = "62c"


class ScrewAxis63c(ScrewAxis):

    """Definition of a 63c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "63c"
