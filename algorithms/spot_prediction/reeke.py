"""Provides a class for producing efficient looping limits for reflection
prediction based on the Reeke algorithm (see Mosflm)"""

from __future__ import absolute_import, division, print_function

import math

from scitbx import matrix


def solve_quad(a, b, c):
    """Robust solution, for real roots only, of a quadratic in the form
    (ax^2 + bx + c)."""

    discriminant = b ** 2 - 4 * a * c

    if discriminant > 0:
        sign = cmp(b, 0)
        if sign == 0:
            sign = 1.0
        q = -0.5 * (b + sign * math.sqrt(discriminant))
        x1 = q / a if a != 0 else None
        x2 = c / q if q != 0 else None
        return [x1, x2]

    elif discriminant == 0:
        return [(-b) / (2 * a)] * 2

    else:
        return [None]


class reeke_model:
    """Model and methods for the Reeke algorithm"""

    def __init__(self, ub_beg, ub_end, axis, s0, dmin, margin=3):

        # the source vector and wavelength
        self._source = -s0
        self._wavelength = 1 / math.sqrt(s0.dot(s0))
        self._wavelength_sq = self._wavelength ** 2

        # the rotation axis
        self._axis = axis

        # the resolution limit
        self._dstarmax = 1 / dmin
        self._dstarmax2 = self._dstarmax ** 2

        # Margin by which to expand limits. Mosflm uses 3.
        self._margin = int(margin)

        # Determine the permutation order of columns of the setting matrix. Use
        # the setting from the beginning for this.
        # As a side-effect set self._permutation.
        col1, col2, col3 = self._permute_axes(ub_beg)

        # Thus set the reciprocal lattice axis vectors, in permuted order
        # p, q and r for both orientations
        rl_vec = [
            ub_beg.extract_block(start=(0, 0), stop=(3, 1)),
            ub_beg.extract_block(start=(0, 1), stop=(3, 2)),
            ub_beg.extract_block(start=(0, 2), stop=(3, 3)),
        ]
        self._rlv_beg = [rl_vec[col1], rl_vec[col2], rl_vec[col3]]
        rl_vec = [
            ub_end.extract_block(start=(0, 0), stop=(3, 1)),
            ub_end.extract_block(start=(0, 1), stop=(3, 2)),
            ub_end.extract_block(start=(0, 2), stop=(3, 3)),
        ]
        self._rlv_end = [rl_vec[col1], rl_vec[col2], rl_vec[col3]]

        # Set permuted setting matrices
        self._p_beg = matrix.sqr(
            self._rlv_beg[0].elems + self._rlv_beg[1].elems + self._rlv_beg[2].elems
        ).transpose()
        self._p_end = matrix.sqr(
            self._rlv_end[0].elems + self._rlv_end[1].elems + self._rlv_end[2].elems
        ).transpose()

        ## Define a new coordinate system concentric with the Ewald sphere.
        ##
        ## X' = X - source_x
        ## Y' = Y - source_y
        ## Z' = Z - source_z
        ##
        ## X = P' h'
        ## -   =  -
        ##                                    / p11 p12 p13 -source_X \
        ## where h' = (p, q, r, 1)^T and P' = | p21 p22 p23 -source_y |
        ##       -                       =    \ p31 p32 p33 -source_z /
        ##

        # Calculate P' matrices for the beginning and end settings
        pp_beg = matrix.rec(
            self._p_beg.elems[0:3]
            + (-1.0 * self._source[0],)
            + self._p_beg.elems[3:6]
            + (-1.0 * self._source[1],)
            + self._p_beg.elems[6:9]
            + (-1.0 * self._source[2],),
            n=(3, 4),
        )
        pp_end = matrix.rec(
            self._p_end.elems[0:3]
            + (-1.0 * self._source[0],)
            + self._p_end.elems[3:6]
            + (-1.0 * self._source[1],)
            + self._p_end.elems[6:9]
            + (-1.0 * self._source[2],),
            n=(3, 4),
        )

        # Various quantities of interest are obtained from the reciprocal metric
        # tensor T of P'. These quantities are to be used (later) for solving
        # the intersection of a line of constant p, q index with the Ewald
        # sphere. It is efficient to calculate these before the outer loop. So,
        # calculate T for both beginning and end settings

        t_beg = (pp_beg.transpose() * pp_beg).as_list_of_lists()
        t_end = (pp_end.transpose() * pp_end).as_list_of_lists()

        # quantities that are constant with p, beginning orientation
        self._cp_beg = [
            t_beg[2][2],
            t_beg[2][3] ** 2,
            t_beg[0][2] * t_beg[2][3] - t_beg[0][3] * t_beg[2][2],
            t_beg[0][2] ** 2 - t_beg[0][0] * t_beg[2][2],
            t_beg[1][2] * t_beg[2][3] - t_beg[1][3] * t_beg[2][2],
            t_beg[0][2] * t_beg[1][2] - t_beg[0][1] * t_beg[2][2],
            t_beg[1][2] ** 2 - t_beg[1][1] * t_beg[2][2],
            2.0 * t_beg[0][2],
            2.0 * t_beg[1][2],
            t_beg[0][0],
            t_beg[1][1],
            2.0 * t_beg[0][1],
            2.0 * t_beg[2][3],
            2.0 * t_beg[1][3],
            2.0 * t_beg[0][3],
        ]

        # quantities that are constant with p, end orientation
        self._cp_end = [
            t_end[2][2],
            t_end[2][3] ** 2,
            t_end[0][2] * t_end[2][3] - t_end[0][3] * t_end[2][2],
            t_end[0][2] ** 2 - t_end[0][0] * t_end[2][2],
            t_end[1][2] * t_end[2][3] - t_end[1][3] * t_end[2][2],
            t_end[0][2] * t_end[1][2] - t_end[0][1] * t_end[2][2],
            t_end[1][2] ** 2 - t_end[1][1] * t_end[2][2],
            2.0 * t_end[0][2],
            2.0 * t_end[1][2],
            t_end[0][0],
            t_end[1][1],
            2.0 * t_end[0][1],
            2.0 * t_end[2][3],
            2.0 * t_end[1][3],
            2.0 * t_end[0][3],
        ]

        ## The following are set during the generation of indices

        # planes of constant p tangential to the Ewald sphere
        self._ewald_p_lim_beg = None
        self._ewald_p_lim_end = None

        # planes of constant p touching the circle of intersection between the
        # Ewald and resolution limiting spheres
        self._res_p_lim_beg = None
        self._res_p_lim_end = None

        # looping p limits
        self._p_lim = None

    def get_source(self):
        return self._source

    def get_ub(self):
        return self._ub

    def get_axis(self):
        return self._axis

    def get_all_p_limits(self):
        """Get both the Ewald and limiting sphere limits for planes of p.
        This is useful for plotting the planes, for example."""

        return (
            self._ewald_p_lim_beg,
            self._ewald_p_lim_end,
            self._res_p_lim_beg,
            self._res_p_lim_end,
        )

    def _permute_axes(self, ub):
        """Find permutation of the columns of an orientation matrix so that
        column p is closest to the source direction, column r is
        closest of q and r to the spindle axis and column q is the remaining
        direction."""

        # Extract the reciprocal lattice directions from the columns of UB
        rl_dirs = [matrix.col(v).normalize() for v in ub.transpose().as_list_of_lists()]

        # Find reciprocal lattice axis closest to source direction
        along_beam = [math.fabs(rl_dirs[j].dot(self._source)) for j in range(3)]
        index_of_p = along_beam.index(max(along_beam))

        # Swap order to put the 'p' axis first
        rl_dirs[0], rl_dirs[index_of_p] = rl_dirs[index_of_p], rl_dirs[0]
        indices = list(range(3))
        indices[0], indices[index_of_p] = indices[index_of_p], indices[0]

        # Now find which of the two remaining reciprocal lattice axes is
        # closest to the rotation axis.
        along_spindle = [math.fabs(rl_dirs[j].dot(self._axis)) for j in (1, 2)]

        index_of_r = along_spindle.index(max(along_spindle)) + 1
        index_of_r = indices[index_of_r]

        # Which is the remaining column index?
        index_of_q = [j for j in range(3) if j not in (index_of_p, index_of_r)][0]

        # permutation matrix such that h, k, l = M * (p, q, r)
        elems = [int(0)] * 9
        elems[3 * index_of_p] = int(1)
        elems[3 * index_of_q + 1] = int(1)
        elems[3 * index_of_r + 2] = int(1)
        self._permutation = matrix.sqr(elems)

        # Return the permuted order of the columns

        return index_of_p, index_of_q, index_of_r

    def _p_limits(self):
        """
        Calculate the values of p at which planes of constant p are tangential
        to the Ewald sphere, and values of p at which planes of constant p touch
        the circle of intersection between the Ewald and resolution limiting
        sphere.

        Note p is the reciprocal cell axis given by the first column of the
        permuted orientation matrix. Set the limits as attributes and return a
        single set of overall limits.
        """

        # Calculate unit vectors normal to planes of constant p, ensuring
        # they point in the direction of increasing p.
        v_beg = self._rlv_beg[1].cross(self._rlv_beg[2]).normalize()

        if self._rlv_beg[0].dot(v_beg) < 0:
            v_beg = -1 * v_beg

        v_end = self._rlv_end[1].cross(self._rlv_end[2]).normalize()

        if self._rlv_end[0].dot(v_end) < 0:
            v_end = -1 * v_end

        # Find distance between the planes of p
        p_dist = abs(self._rlv_beg[0].dot(v_beg))

        # Find distances between p = 0 and the plane passing through the
        # centre of the Ewald sphere
        dp_beg = abs(v_beg.dot(self._source))
        dp_end = abs(v_end.dot(self._source))

        # There are two planes of constant p that are tangential to the Ewald
        # sphere, on either side of the sphere. The smaller in magnitude of p
        # is the number of planes that fit in one radius of the Ewald sphere
        # minus the number of planes between the centre of the Ewald sphere
        # and the p=0 plane (a diagram helps!). The larger is the number of
        # planes in one radius of the Ewald sphere *plus* the the number of
        # planes between the centre of the Ewald sphere and p = 0.
        #
        # The correct sign is determined by whether the plane normal vector is
        # more closely parallel or antiparallel to the beam direction.

        sign = cmp(v_beg.dot(self._source), 0)

        limits = [
            (sign * s * (self._source.length() + s * dp_beg) / p_dist) for s in (-1, 1)
        ]

        self._ewald_p_lim_beg = tuple(sorted(limits))

        sign = cmp(v_end.dot(self._source), 0)

        limits = [
            (sign * s * (self._source.length() + s * dp_end) / p_dist) for s in (-1, 1)
        ]

        self._ewald_p_lim_end = tuple(sorted(limits))

        # Now determine limits for the planes of p that touch the circle of
        # intersection between the Ewald and resolution limiting spheres

        # FIXME is there a more efficient way to get sin_2theta?
        sin_theta = 0.5 * self._wavelength * self._dstarmax
        assert abs(sin_theta) <= 1.0  # sanity check
        sin_2theta = math.sin(2.0 * math.asin(sin_theta))

        e = 2.0 * sin_theta ** 2 * dp_beg
        f = sin_2theta * math.sqrt(max(1.0 / self._wavelength_sq - dp_beg ** 2, 0.0))
        limits = [(sign * e + s * f) / p_dist for s in (-1, 1)]

        self._res_p_lim_beg = tuple(sorted(limits))

        e = 2.0 * sin_theta ** 2 * dp_end
        f = sin_2theta * math.sqrt(max(1.0 / self._wavelength_sq - dp_end ** 2, 0))
        limits = [(sign * e + s * f) / p_dist for s in (-1, 1)]

        self._res_p_lim_end = tuple(sorted(limits))

        # select between Ewald and resolution limits on the basis of sign
        if sign < 0:  # p axis aligned with beam, against source

            p_min_beg = max(min(self._res_p_lim_beg), min(self._ewald_p_lim_beg))
            p_min_end = max(min(self._res_p_lim_end), min(self._ewald_p_lim_end))

            p_max_beg = max(max(self._res_p_lim_beg), max(self._ewald_p_lim_beg))
            p_max_end = max(max(self._res_p_lim_end), max(self._ewald_p_lim_end))

        else:  # p axis aligned with source, against beam

            p_min_beg = min(min(self._res_p_lim_beg), min(self._ewald_p_lim_beg))
            p_min_end = min(min(self._res_p_lim_end), min(self._ewald_p_lim_end))

            p_max_beg = min(max(self._res_p_lim_beg), max(self._ewald_p_lim_beg))
            p_max_end = min(max(self._res_p_lim_end), max(self._ewald_p_lim_end))

        p_lim_beg = (p_min_beg, p_max_beg)
        p_lim_end = (p_min_end, p_max_end)
        # p_lim_beg = sorted(self._ewald_p_lim_beg + self._res_p_lim_beg)[1:3]
        # p_lim_end = sorted(self._ewald_p_lim_end + self._res_p_lim_end)[1:3]

        # single set of limits covering overall range
        p_lim = sorted(p_lim_beg + p_lim_end)[0::3]
        p_lim[0] = int(p_lim[0]) - self._margin
        p_lim[1] = int(p_lim[1]) + self._margin

        return p_lim

    def _q_limits(self, p):
        """Calculate the values of q at which lines of constant p, q are
        tangential to the circle intersecting the Ewald sphere at plane p,
        and values of q at which lines of constant p, q are tangential to
        the circle intersecting the resolution limiting sphere at plane p.i
        Return the appropriate overall limits."""

        # First the resolution limits. Set up the quadratic to solve
        a = self._cp_beg[6]
        b = 2.0 * p * self._cp_beg[5]
        c = p ** 2 * self._cp_beg[3] + self._cp_beg[0] * self._dstarmax2

        res_q_lim = solve_quad(a, b, c)
        res_q_lim = sorted([item for item in res_q_lim if item is not None])
        if len(res_q_lim) == 0:
            return None

        # Extend limits by the margin, ensuring there is a range even for
        # a single quadratic root
        res_q_lim = [
            int(res_q_lim[0]) - max(self._margin, 1),
            int(res_q_lim[-1]) + max(self._margin, 1),
        ]

        # Ewald sphere limits for the beginning setting
        b = 2.0 * (self._cp_beg[4] + p * self._cp_beg[5])
        c = self._cp_beg[1] + p * (2 * self._cp_beg[2] + p * self._cp_beg[3])

        ewald_q_lim_beg = solve_quad(a, b, c)

        # Ewald sphere limits for the end setting
        a = self._cp_end[6]
        b = 2.0 * (self._cp_end[4] + p * self._cp_end[5])
        c = self._cp_end[1] + p * (2 * self._cp_end[2] + p * self._cp_end[3])

        ewald_q_lim_end = solve_quad(a, b, c)

        # Determine the overall Ewald limits
        ewald_q_lim = sorted(
            [item for item in ewald_q_lim_beg + ewald_q_lim_end if item is not None]
        )
        if len(ewald_q_lim) > 0:
            ewald_q_lim = [
                int(ewald_q_lim[0]) - max(self._margin, 1),
                int(ewald_q_lim[-1]) + max(self._margin, 1),
            ]

        else:
            return None

        # Choose most restrictive of Ewald and res limits.
        q_lim = sorted(res_q_lim + ewald_q_lim)
        q_lim = [q_lim[1], q_lim[2]]

        return q_lim

    def _r_limits(self, p, q, cq_beg, cq_end):
        """Calculate the values of r at which lines of constant p, q intersect
        the resolution limiting and the Ewald spheres, and return the
        appropriate overall limits"""

        # First the resolution limits. Set up the quadratic to solve
        a = self._cp_beg[0]
        b = cq_beg[0] + q * self._cp_beg[8]
        c = cq_beg[1] + q ** 2 * self._cp_beg[10] + q * cq_beg[2] - self._dstarmax2

        res_r_lim = solve_quad(a, b, c)
        res_r_lim = sorted([item for item in res_r_lim if item is not None])
        if len(res_r_lim) == 0:
            return None

        # Extend limits by the margin, ensuring there is a range even for
        # a single quadratic root
        res_r_lim = [
            int(res_r_lim[0]) - max(self._margin, 1),
            int(res_r_lim[-1]) + max(self._margin, 1),
        ]

        # Ewald sphere limits for the beginning setting
        b = cq_beg[0] + q * self._cp_beg[8] + self._cp_beg[12]
        c = (
            cq_beg[1]
            + q * (cq_beg[2] + self._cp_beg[13])
            + q ** 2 * self._cp_beg[10]
            + cq_beg[3]
        )

        ewald_r_lim_beg = solve_quad(a, b, c)
        ewald_r_lim_beg = [item for item in ewald_r_lim_beg if item is not None]

        # Ewald sphere limits for the end setting
        a = self._cp_end[0]
        b = cq_end[0] + q * self._cp_end[8] + self._cp_end[12]
        c = (
            cq_end[1]
            + q * (cq_end[2] + self._cp_end[13])
            + q ** 2 * self._cp_end[10]
            + cq_end[3]
        )

        ewald_r_lim_end = solve_quad(a, b, c)
        ewald_r_lim_end = [item for item in ewald_r_lim_end if item is not None]

        # if no intersections at all, return None

        if len(ewald_r_lim_beg) == 0 and len(ewald_r_lim_end) == 0:
            return None

        # if there are no intersections at the beginning setting, set up a
        # single loop covering the range between the intersections at the end
        # setting, and vice versa.
        if len(ewald_r_lim_beg) == 0:

            l1 = [
                int(min(ewald_r_lim_end)) - max(self._margin, 1),
                int(max(ewald_r_lim_end)) + max(self._margin, 1),
            ]
            l2 = [None]

        elif len(ewald_r_lim_end) == 0:

            l1 = [
                int(min(ewald_r_lim_beg)) - max(self._margin, 1),
                int(max(ewald_r_lim_beg)) + max(self._margin, 1),
            ]
            l2 = [None]

        # otherwise there is at least one intersection at both settings.
        # Set up two loops, one for each range swept out by a point of
        # intersection as it travels from the beginning to the end setting.
        else:

            l1 = sorted([min(ewald_r_lim_beg), min(ewald_r_lim_end)])
            l1 = [int(l1[0]) - max(self._margin, 1), int(l1[1]) + max(self._margin, 1)]
            l2 = sorted([max(ewald_r_lim_beg), max(ewald_r_lim_end)])
            l2 = [int(l2[0]) - max(self._margin, 1), int(l2[1]) + max(self._margin, 1)]

        # restrict loops according to the resolution limit
        l1[0] = max(res_r_lim[0], l1[0])
        l1[1] = min(res_r_lim[1], l1[1])
        if l1[0] >= l1[1]:
            l1 = [None]

        if l2 != [None]:
            l2[0] = max(res_r_lim[0], l2[0])
            l2[1] = min(res_r_lim[1], l2[1])
            if l2[0] >= l2[1]:
                l2 = [None]

        if l1 == [None] and l2 == [None]:
            return None

        return [tuple(l1), tuple(l2)]

    def generate_indices(self):
        """Determine looping limits for indices h, k and l using the Reeke
        algorithm. This is the top level method for this module. All other
        methods are (probably) called by this, and therefore may as well be
        private."""

        # The outer loop is between limits for the axis most closely parallel,
        # or antiparallel, to the X-ray beam, which is called 'p'.

        # Determine the limiting values of p
        p_lim = self._p_limits()

        # fill indices list by looping over p, q and r
        hkl = []

        for p in range(p_lim[0], p_lim[1] + 1):

            # quantities that vary with p but are constant with q, beginning setting
            cq_beg = [
                (p * self._cp_beg[7]),
                (p ** 2 * self._cp_beg[9]),
                (p * self._cp_beg[11]),
                (p * self._cp_beg[14]),
            ]

            # quantities that vary with p but are constant with q, end setting
            cq_end = [
                (p * self._cp_end[7]),
                (p ** 2 * self._cp_end[9]),
                (p * self._cp_end[11]),
                (p * self._cp_end[14]),
            ]

            # find the limiting values of q
            q_lim = self._q_limits(p)
            if q_lim is None:
                continue

            for q in range(q_lim[0], q_lim[1] + 1):

                # find the limiting values of r
                r_lim = self._r_limits(p, q, cq_beg, cq_end)
                if r_lim is None:
                    continue

                # make list of trials, removing any duplicates
                r_trials = []
                for item in r_lim:
                    if item[0] is None:
                        continue

                    r_seq = range(item[0], item[1] + 1)
                    r_trials += [e for e in r_seq if e not in r_trials]

                for r in r_trials:
                    hkl.append((self._permutation * (p, q, r)).elems)

        return hkl
