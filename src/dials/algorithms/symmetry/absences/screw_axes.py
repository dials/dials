"""Definitions of screw axes with methods for scoring against data."""

from __future__ import annotations

import logging
import math

import numpy as np
from jinja2 import ChoiceLoader, Environment, PackageLoader
from scipy.stats import norm

from scitbx.array_family import flex

from dials.algorithms.symmetry.absences.plots import plot_screw_axes
from dials.util.observer import Observer, Subject, singleton

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
            "fourier_space_data": screw_axis.fourier_space_data,
            "axis_repeat": screw_axis.axis_repeat,
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
        template = env.get_template("simple_report.html")
        html = template.render(
            page_title="DIALS systematic absences report",
            panel_title="Screw axes analysis",
            panel_id="screw_axes",
            graphs=self.data["screw_axes"],
        )
        with open(filename, "wb") as f:
            f.write(html.encode("utf-8", "xmlcharrefreplace"))


class ScrewAxis(Subject):
    """Definition of a generic screw axis."""

    axis_idx = None  # x=0, y=1, z=2
    axis_repeat = None  # repeat of present reflections e.g =4 for 41, =2 for 42
    name = None

    def __init__(self):
        super().__init__(events=["scored axis"])
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
        self.fourier_space_data = {}

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

    def get_all_suitable_reflections(self, reflection_table):
        """Select suitable reflections for testing the screw axis."""
        refl = reflection_table
        sel = self.select_axial_reflections(refl["miller_index"])
        miller_idx = refl["miller_index"].select(sel)
        self.miller_axis_vals = miller_idx.as_vec3_double().parts()[self.axis_idx]
        self.intensities = refl["intensity"].select(sel)
        self.sigmas = flex.sqrt(refl["variance"].select(sel))
        self.i_over_sigma = self.intensities / self.sigmas

        if self.equivalent_axes:
            for a in self.equivalent_axes:
                sel = a.select_axial_reflections(refl["miller_index"])
                miller_idx = refl["miller_index"].select(sel)
                intensities = refl["intensity"].select(sel)
                sigmas = flex.sqrt(refl["variance"].select(sel))
                self.i_over_sigma.extend(intensities / sigmas)
                self.miller_axis_vals.extend(
                    miller_idx.as_vec3_double().parts()[a.axis_idx]
                )
                self.intensities.extend(intensities)
                self.sigmas.extend(sigmas)

    @Subject.notify_event(event="scored axis")
    def score_axis(self, reflection_table, significance_level=0.95, method="direct"):
        """Score the axis given a reflection table of data."""
        if method == "direct":
            return self.score_axis_direct(reflection_table, significance_level)
        else:
            return self.score_axis_fourier(reflection_table, significance_level)

    @staticmethod
    def _score_axis_fourier(miller_index, i_over_sigma, axis_repeat):
        """
        Score screw axis based on selected miller indices and I over SigI values

        Parameters
        ----------
        miller_index : np.array
            Vector of miller indices along the chosen axis.
        i_over_sigma : np.array
            Vector of signal to noise ratios with the same shape as miller_indices.
        axis_repeat : int
            The expected periodicity of the screw axis. Must be one of {2, 3, 4, 6}.
        """
        if axis_repeat not in {2, 3, 4, 6}:
            raise ValueError(
                f"Received axis_repeat, {axis_repeat}, but expected an integer in {2, 3, 4, 6}"
            )

        # We must take care to make sure the length of the fourier transformed vector is divisible by 6
        miller_index = miller_index - miller_index.min()
        n = miller_index.max() + 1
        n = 6 * ((n // 6) + 1)
        direct_space = np.zeros(n)
        direct_space[miller_index] = i_over_sigma

        # Fourier transform i over sigma
        fourier_space = np.abs(np.fft.fft(direct_space))

        # Indices for Fourier frequencies which may correspond to screw periodicities
        # These correspond to absences every 0, 2, 3, 4, and 6 reflections.
        screw_idx = [0]
        if axis_repeat == 2:
            screw_idx += [n // 2]
        elif axis_repeat == 3:
            screw_idx += [n // 3, -n // 3]
        elif axis_repeat == 4:
            screw_idx += [n // 4, n // 2, -n // 4]
        else:
            screw_idx += [n // 6, n // 3, n // 2, -n // 3, -n // 6]

        null_idx = np.ones_like(fourier_space, dtype=bool)
        null_idx[screw_idx] = False

        # To determine the probability of a screw axis, use the frequencies which do not
        # correspond to any screw axis periodicity to form a null model. The ask what
        # the probability of the candidate frequency is under the null model.
        # In this case, we will parameterize the null frequencies by a normal distribution.

        mean = fourier_space[null_idx].mean()
        std = fourier_space[null_idx].std()
        p_screw = norm.cdf(fourier_space[n // axis_repeat], loc=mean, scale=std)
        fourier_space_data = {"fourier_space": fourier_space, "n": n}
        return p_screw, fourier_space_data

    def score_axis_fourier(self, reflection_table, significance_level=0.95):
        """Estimate the probability of a screw axis using Fourier analysis."""

        self.get_all_suitable_reflections(reflection_table)
        expected_sel = self.miller_axis_vals.iround() % self.axis_repeat == 0

        expected = self.i_over_sigma.select(expected_sel)
        expected_abs = self.i_over_sigma.select(~expected_sel)
        self.n_refl_used = (expected.size(), expected_abs.size())

        if not expected or not expected_abs:
            return 0.0

        # Log these quantities for reporting.
        self.mean_I_sigma_abs = flex.mean(expected_abs)
        self.mean_I_sigma = flex.mean(expected)
        self.mean_I = flex.mean(self.intensities.select(expected_sel))
        self.mean_I_abs = flex.mean(self.intensities.select(~expected_sel))

        i_over_sigma = np.array(self.i_over_sigma)
        miller_index = np.array(self.miller_axis_vals.iround())
        p_screw, fourier_space_data = self._score_axis_fourier(
            miller_index, i_over_sigma, self.axis_repeat
        )
        # Record the fourier data for reporting.
        self.fourier_space_data = fourier_space_data
        # Zero out the probability if it is less than the significance_level
        if p_screw < significance_level:
            return 0.0

        return p_screw

    def score_axis_direct(self, reflection_table, significance_level=0.95):
        """Score the axis given a reflection table of data."""
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
        P_absent = 0.5 * (1.0 + math.erf(z_score_absent / (2**0.5)))
        P_present = 0.5 * (1.0 + math.erf(z_score_present / (2**0.5)))

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
            outlier_msg = (
                """Screw axis %s could only be assigned after removing a suspected outlier
from the expected 'absent' reflections."""
                % self.name
            )
            if expected_abs.size() <= 1:
                logger.info(outlier_msg)
                return P_present
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
                logger.info(outlier_msg)
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
