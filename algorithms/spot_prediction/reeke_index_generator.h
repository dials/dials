/*
 * reeke_index_generator.h
 *
 *  Copyright (C) 2013 CCP4, Diamond Light Source
 *
 *  Author: David Waterman (python code)
 *  Author: James Parkhurst (c++ port)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_REEKE_INDEX_GENERATOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_REEKE_INDEX_GENERATOR_H

#include <cmath>
#include <algorithm>
#include <boost/optional.hpp>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/small.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/spot_prediction/scan_varying_helpers.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using af::max_index;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::matrix::transpose_multiply;

  namespace reeke_detail {

    /** Sort a 2 element vector */
    template <typename T>
    void sort2(vec2<T>& a) {
      if (a[0] > a[1]) {
        std::swap(a[0], a[1]);
      }
    }

    /**
     * Struct to permute axes
     */
    struct permute_axes {
      /**
       * Find permutation of the columns of an orientation matrix so that column
       * p is closest to the source direction, column r is closest of q and r to
       * the spindle axis and column q is the remaining direction.
       */
      permute_axes(mat3<double> ub, vec3<double> axis, vec3<double> source)
          : index(0, 1, 2), permutation(0, 0, 0, 0, 0, 0, 0, 0, 0) {
        // Extract the reciprocal lattice directions from the columns of UB
        vec3<double> dir1(ub[0], ub[3], ub[6]);
        vec3<double> dir2(ub[1], ub[4], ub[7]);
        vec3<double> dir3(ub[2], ub[5], ub[8]);
        DIALS_ASSERT(dir1.length() > 0.0);
        DIALS_ASSERT(dir2.length() > 0.0);
        DIALS_ASSERT(dir3.length() > 0.0);
        // vec3.normalize will cause division by zero error if length is zero
        vec3<double> rl_dirs[3] = {
          dir1.normalize(), dir2.normalize(), dir3.normalize()};

        // Find reciprocal lattice axis closest to source direction and swap the
        // index order to put the 'p' axis first
        std::size_t p_index = max_index(vec3<double>(std::abs(rl_dirs[0] * source),
                                                     std::abs(rl_dirs[1] * source),
                                                     std::abs(rl_dirs[2] * source))
                                          .const_ref());
        std::swap(index[0], index[p_index]);

        // Now find which of the two remaining reciprocal lattice axes is closest
        // to the rotation axis and swap the index order to put r in last place
        std::size_t r_index = max_index(vec2<double>(std::abs(rl_dirs[index[1]] * axis),
                                                     std::abs(rl_dirs[index[2]] * axis))
                                          .const_ref());
        std::swap(index[2], index[r_index + 1]);

        // permutation matrix such that h, k, l = M * (p, q, r)
        permutation[3 * index[0] + 0] = 1;
        permutation[3 * index[1] + 1] = 1;
        permutation[3 * index[2] + 2] = 1;
      }

      vec3<std::size_t> index;
      mat3<std::size_t> permutation;
    };

    /**
     * Struct to permute the matrix
     */
    struct permute_matrix : public permute_axes {
      permute_matrix(mat3<double> ub_beg,
                     mat3<double> ub_end,
                     vec3<double> axis,
                     vec3<double> source)
          : permute_axes(ub_beg, axis, source) {
        // Set the reciprocal lattice axis vectors, in permuted order p, q and r
        // for both orientations
        rlv_beg = mat3<double>(ub_beg[index[0]],
                               ub_beg[index[0] + 3],
                               ub_beg[index[0] + 6],
                               ub_beg[index[1]],
                               ub_beg[index[1] + 3],
                               ub_beg[index[1] + 6],
                               ub_beg[index[2]],
                               ub_beg[index[2] + 3],
                               ub_beg[index[2] + 6]);
        rlv_end = mat3<double>(ub_end[index[0]],
                               ub_end[index[0] + 3],
                               ub_end[index[0] + 6],
                               ub_end[index[1]],
                               ub_end[index[1] + 3],
                               ub_end[index[1] + 6],
                               ub_end[index[2]],
                               ub_end[index[2] + 3],
                               ub_end[index[2] + 6]);
      }

      mat3<double> rlv_beg;
      mat3<double> rlv_end;
    };

    /**
     * Struct to compute variables constant with p
     */
    struct compute_constant_with_p {
      compute_constant_with_p(mat3<double> rlv_beg,
                              mat3<double> rlv_end,
                              vec3<double> axis,
                              vec3<double> source_beg,
                              vec3<double> source_end) {
        // Set permuted setting matrices
        mat3<double> p_beg = rlv_beg.transpose();
        mat3<double> p_end = rlv_end.transpose();

        /*
         * Define a new coordinate system concentric with the Ewald sphere.
         *
         * X' = X - source_x
         * Y' = Y - source_y
         * Z' = Z - source_z
         *
         * X = P' h'
         * -   =  -
         *                                    / p11 p12 p13 -source_X \
         * where h' = (p, q, r, 1)^T and P' = | p21 p22 p23 -source_y |
         *       -                       =    \ p31 p32 p33 -source_z /
         */

        // Calculate P' matrices for the beginning and end states
        double pp_beg[12] = {p_beg[0],
                             p_beg[1],
                             p_beg[2],
                             -source_beg[0],
                             p_beg[3],
                             p_beg[4],
                             p_beg[5],
                             -source_beg[1],
                             p_beg[6],
                             p_beg[7],
                             p_beg[8],
                             -source_beg[2]};
        double pp_end[12] = {p_end[0],
                             p_end[1],
                             p_end[2],
                             -source_end[0],
                             p_end[3],
                             p_end[4],
                             p_end[5],
                             -source_end[1],
                             p_end[6],
                             p_end[7],
                             p_end[8],
                             -source_end[2]};

        // Various quantities of interest are obtained from the reciprocal metric
        // tensor T of P'. These quantities are to be used (later) for solving the
        // intersection of a line of constant p, q index with the Ewald sphere. It
        // is efficient to calculate these before the outer loop. So, calculate T
        // for both beginning and end settings
        double t_beg[16], t_end[16];
        transpose_multiply(pp_beg, pp_beg, 3, 4, 4, t_beg);
        transpose_multiply(pp_end, pp_end, 3, 4, 4, t_end);

        // quantities that are constant with p, beginning setting
        cp_beg.resize(15);
        cp_beg[0] = t_beg[10];
        cp_beg[1] = t_beg[11] * t_beg[11];
        cp_beg[2] = t_beg[2] * t_beg[11] - t_beg[3] * t_beg[10];
        cp_beg[3] = t_beg[2] * t_beg[2] - t_beg[0] * t_beg[10];
        cp_beg[4] = t_beg[6] * t_beg[11] - t_beg[7] * t_beg[10];
        cp_beg[5] = t_beg[2] * t_beg[6] - t_beg[1] * t_beg[10];
        cp_beg[6] = t_beg[6] * t_beg[6] - t_beg[5] * t_beg[10];
        cp_beg[7] = 2.0 * t_beg[2];
        cp_beg[8] = 2.0 * t_beg[6];
        cp_beg[9] = t_beg[0];
        cp_beg[10] = t_beg[5];
        cp_beg[11] = 2.0 * t_beg[1];
        cp_beg[12] = 2.0 * t_beg[11];
        cp_beg[13] = 2.0 * t_beg[7];
        cp_beg[14] = 2.0 * t_beg[3];

        // quantities that are constant with p, end setting
        cp_end.resize(15);
        cp_end[0] = t_end[10];
        cp_end[1] = t_end[11] * t_end[11];
        cp_end[2] = t_end[2] * t_end[11] - t_end[3] * t_end[10];
        cp_end[3] = t_end[2] * t_end[2] - t_end[0] * t_end[10];
        cp_end[4] = t_end[6] * t_end[11] - t_end[7] * t_end[10];
        cp_end[5] = t_end[2] * t_end[6] - t_end[1] * t_end[10];
        cp_end[6] = t_end[6] * t_end[6] - t_end[5] * t_end[10];
        cp_end[7] = 2.0 * t_end[2];
        cp_end[8] = 2.0 * t_end[6];
        cp_end[9] = t_end[0];
        cp_end[10] = t_end[5];
        cp_end[11] = 2.0 * t_end[1];
        cp_end[12] = 2.0 * t_end[11];
        cp_end[13] = 2.0 * t_end[7];
        cp_end[14] = 2.0 * t_end[3];
      }

      af::small<double, 15> cp_beg;
      af::small<double, 15> cp_end;
    };

  }  // namespace reeke_detail

  /**
   * Class implementing the reeke model.
   */
  class ReekeModel {
  public:
    /**
     * Initialise and compute p limits.
     * @param ub_beg The beginning UB matrix
     * @param ub_end The end UB matrix
     * @param axis The rotation axis
     * @param source_beg The beginning source vector
     * @param source_end The end source vector
     * @param dmin The resolution
     * @param margin The margin to add around limits
     */
    ReekeModel(mat3<double> ub_beg,
               mat3<double> ub_end,
               vec3<double> axis,
               vec3<double> source_beg,
               vec3<double> source_end,
               double dmin,
               int margin) {
      // Save the source and axis
      source_beg_ = source_beg;
      source_end_ = source_end;
      axis_ = axis;
      margin_ = margin;

      // Set the wavelength
      wavelength_beg_ = 1.0 / source_beg_.length();
      wavelength_sq_beg_ = wavelength_beg_ * wavelength_beg_;
      wavelength_end_ = 1.0 / source_end_.length();
      wavelength_sq_end_ = wavelength_end_ * wavelength_end_;

      // the resolution limit
      dstarmax_ = 1.0 / dmin;
      dstarmax2_ = dstarmax_ * dstarmax_;

      // Determine the permutation order of columns of the setting matrix. Use
      // the setting ans source vector from the beginning for this. As a
      // side-effect set self._permutation.
      reeke_detail::permute_matrix perm(ub_beg, ub_end, axis, source_beg_);
      permutation_ = perm.permutation;

      // Compute the variables that are constant with p
      reeke_detail::compute_constant_with_p compute_cp(
        perm.rlv_beg, perm.rlv_end, axis, source_beg_, source_end_);
      cp_beg_ = compute_cp.cp_beg;
      cp_end_ = compute_cp.cp_end;

      // Compute and initialise the p limits
      compute_p_limits(perm.rlv_beg, perm.rlv_end);
    }

    /**
     * @returns the permutation
     */
    mat3<std::size_t> permutation() const {
      return permutation_;
    }

    /**
     * @returns The ewald sphere p limits
     */
    std::pair<vec2<double>, vec2<double> > ewald_sphere_p_limits() const {
      return std::make_pair(ewald_p_lim_beg_, ewald_p_lim_end_);
    }

    /**
     * @returns The resolution p limits.
     */
    std::pair<vec2<double>, vec2<double> > resolution_p_limits() const {
      return std::make_pair(res_p_lim_beg_, res_p_lim_end_);
    }

    /**
     * @returns The p limits
     */
    vec2<int> p_limits() const {
      return p_lim_;
    }

    /**
     * Calculate the values of q at which lines of constant p, q are tangential
     * to the circle intersecting the Ewald sphere at plane p, and values of q
     * at which lines of constant p, q are tangential to the circle intersecting
     * the resolution limiting sphere at plane p.i Return the appropriate
     * overall limits.
     * @returns The q limits
     */
    vec2<int> q_limits(int p) {
      // Get the resolution limits and the ewald sphere limits
      boost::optional<vec2<int> > res_q_lim = resolution_q_limits(p);
      boost::optional<vec2<int> > ewald_q_lim = ewald_sphere_q_limits(p);
      if (!res_q_lim || !ewald_q_lim) {
        return vec2<int>(0, 0);
      }

      // Choose most restrictive of Ewald and res limits.
      af::tiny<int, 4> limits(
        (*ewald_q_lim)[0], (*ewald_q_lim)[1], (*res_q_lim)[0], (*res_q_lim)[1]);
      std::sort(limits.begin(), limits.end());
      return vec2<int>(limits[1], limits[2] + 1);
    }

    /**
     * Calculate the values of r at which lines of constant p, q intersect
     * the resolution limiting and the Ewald spheres, and return the
     * appropriate overall limits
     * @returns The r limits
     */
    af::small<vec2<int>, 2> r_limits(double p, double q) {
      af::small<vec2<int>, 2> result;

      // quantities that vary with p but are constant with q, beginning setting
      af::tiny<double, 4> cq_beg(4);
      cq_beg[0] = p * cp_beg_[7];
      cq_beg[1] = p * p * cp_beg_[9];
      cq_beg[2] = p * cp_beg_[11];
      cq_beg[3] = p * cp_beg_[14];

      // quantities that vary with p but are constant with q, end setting
      af::tiny<double, 4> cq_end(4);
      cq_end[0] = p * cp_end_[7];
      cq_end[1] = p * p * cp_end_[9];
      cq_end[2] = p * cp_end_[11];
      cq_end[3] = p * cp_end_[14];

      // Get the resolution limits
      boost::optional<vec2<int> > res_r_lim = resolution_r_limits(p, q, cq_beg);
      if (!res_r_lim) {
        return result;
      }

      // Get the ewald sphere limits and restrict loops according to the
      // resolution limit
      af::small<vec2<int>, 2> ewald_r_lim = ewald_sphere_r_limits(p, q, cq_beg, cq_end);
      for (std::size_t i = 0; i < ewald_r_lim.size(); ++i) {
        if ((*res_r_lim)[0] > ewald_r_lim[i][0]) {
          ewald_r_lim[i][0] = (*res_r_lim)[0];
        }
        if ((*res_r_lim)[1] < ewald_r_lim[i][1]) {
          ewald_r_lim[i][1] = (*res_r_lim)[1];
        }
        if (ewald_r_lim[i][0] < ewald_r_lim[i][1]) {
          result.push_back(vec2<int>(ewald_r_lim[i][0], ewald_r_lim[i][1] + 1));
        }
      }

      // Ensure that if there are two ranges, they are in order and
      // non-overlapping.
      if (result.size() == 2) {
        if (result[1][0] < result[0][0]) {
          std::swap(result[0], result[1]);
        }
        if (result[1][0] < result[0][1]) {
          result[1][0] = result[0][1];
        }
      }
      return result;
    }

  private:
    /**
     * There are two planes of constant p that are tangential to the Ewald
     * sphere, on either side of the sphere. The smaller in magnitude of p is
     * the number of planes that fit in one radius of the Ewald sphere minus the
     * number of planes between the centre of the Ewald sphere and the p=0 plane
     * (a diagram helps!). The larger is the number of planes in one radius of
     * the Ewald sphere *plus* the the number of planes between the centre of
     * the Ewald sphere and p = 0.
     *
     * The correct sign is determined by whether the plane normal vector is more
     * closely parallel or antiparallel to the beam direction.
     */
    void compute_ewald_sphere_p_limits(vec3<double> v_beg,
                                       vec3<double> v_end,
                                       double dp_beg,
                                       double dp_end,
                                       double p_dist) {
      // Calculate beg limit
      int sign = v_beg * source_beg_ >= 0 ? 1 : -1;
      ewald_p_lim_beg_ = vec2<double>(-sign * (source_beg_.length() - dp_beg) / p_dist,
                                      sign * (source_beg_.length() + dp_beg) / p_dist);
      reeke_detail::sort2(ewald_p_lim_beg_);

      // Calculate end limit
      sign = v_end * source_end_ >= 0 ? 1 : -1;
      ewald_p_lim_end_ = vec2<double>(-sign * (source_end_.length() - dp_end) / p_dist,
                                      sign * (source_end_.length() + dp_end) / p_dist);
      reeke_detail::sort2(ewald_p_lim_end_);
    }

    /**
     * Compute the resolution p limits
     */
    void compute_resolution_p_limits(vec3<double> v_beg,
                                     vec3<double> v_end,
                                     double dp_beg,
                                     double dp_end,
                                     double p_dist) {
      // Determine limits for the planes of p that touch the circle of
      // intersection between the Ewald and resolution limiting spheres
      double sin_theta = 0.5 * wavelength_beg_ * dstarmax_;
      DIALS_ASSERT(sin_theta <= 1.0 && sin_theta >= -1.0);  // sanity check
      double sin_2theta = std::sin(2.0 * std::asin(sin_theta));

      // Calculate beg limit
      int sign = v_end * source_beg_ >= 0 ? 1 : -1;
      double e = 2.0 * sin_theta * sin_theta * dp_beg;
      double f = sin_2theta
                 * std::sqrt(std::max(1.0 / wavelength_sq_beg_ - dp_beg * dp_beg, 0.0));
      res_p_lim_beg_ = vec2<double>((sign * e - f) / p_dist, (sign * e + f) / p_dist);
      reeke_detail::sort2(res_p_lim_beg_);

      sin_theta = 0.5 * wavelength_end_ * dstarmax_;
      DIALS_ASSERT(sin_theta <= 1.0 && sin_theta >= -1.0);  // sanity check
      sin_2theta = std::sin(2.0 * std::asin(sin_theta));

      // Calculate end limit
      e = 2.0 * sin_theta * sin_theta * dp_end;
      f = sin_2theta
          * std::sqrt(std::max(1.0 / wavelength_sq_end_ - dp_end * dp_end, 0.0));
      res_p_lim_end_ = vec2<double>((sign * e - f) / p_dist, (sign * e + f) / p_dist);
      reeke_detail::sort2(res_p_lim_end_);
    }

    /**
     * Calculate the values of p at which planes of constant p are tangential
     * to the Ewald sphere, and values of p at which planes of constant p touch
     * the circle of intersection between the Ewald and resolution limiting
     * sphere.

     * Note p is the reciprocal cell axis given by the first column of the
     * permuted orientation matrix. Set the limits as attributes and return a
     * single set of overall limits.
     */
    void compute_p_limits(mat3<double> rlv_beg, mat3<double> rlv_end) {
      // Get the rows of the beginning and end matrices
      vec3<double> rlv_beg0 = rlv_beg.get_row(0);
      vec3<double> rlv_beg1 = rlv_beg.get_row(1);
      vec3<double> rlv_beg2 = rlv_beg.get_row(2);
      vec3<double> rlv_end0 = rlv_end.get_row(0);
      vec3<double> rlv_end1 = rlv_end.get_row(1);
      vec3<double> rlv_end2 = rlv_end.get_row(2);

      // Calculate unit vectors normal to planes of constant p, ensuring they
      // point in the direction of increasing p.
      vec3<double> v_beg = rlv_beg1.cross(rlv_beg2).normalize();
      vec3<double> v_end = rlv_end1.cross(rlv_end2).normalize();
      if (rlv_beg0 * v_beg < 0.0) v_beg = -v_beg;
      if (rlv_end0 * v_end < 0.0) v_end = -v_end;

      // Find distance between the planes of p and find distances between p = 0
      // and the plane passing through the centre of the Ewald sphere
      double p_dist = std::abs(rlv_beg0 * v_beg);
      double dp_beg = std::abs(v_beg * source_beg_);
      double dp_end = std::abs(v_end * source_end_);

      // Compute the ewald sphere and resolution limits
      compute_ewald_sphere_p_limits(v_beg, v_end, dp_beg, dp_end, p_dist);
      compute_resolution_p_limits(v_beg, v_end, dp_beg, dp_end, p_dist);

      // select between Ewald and resolution limits on the basis of sign
      af::tiny<double, 4> limits;
      int sign = v_end * source_end_ >= 0 ? 1 : -1;
      if (sign < 0) {
        // p axis aligned with beam, against source
        limits[0] = std::max(res_p_lim_beg_[0], ewald_p_lim_beg_[0]);
        limits[1] = std::max(res_p_lim_end_[0], ewald_p_lim_end_[0]);
        limits[2] = std::max(res_p_lim_beg_[1], ewald_p_lim_beg_[1]);
        limits[3] = std::max(res_p_lim_end_[1], ewald_p_lim_end_[1]);
      } else {
        // p axis aligned with source, against beam
        limits[0] = std::min(res_p_lim_beg_[0], ewald_p_lim_beg_[0]);
        limits[1] = std::min(res_p_lim_end_[0], ewald_p_lim_end_[0]);
        limits[2] = std::min(res_p_lim_beg_[1], ewald_p_lim_beg_[1]);
        limits[3] = std::min(res_p_lim_end_[1], ewald_p_lim_end_[1]);
      }

      // single set of limits covering overall range
      p_lim_ = vec2<int>((int)af::min(limits.const_ref()) - margin_,
                         (int)af::max(limits.const_ref()) + margin_ + 1);
    }

    /**
     * @returns The resolution q limits
     */
    boost::optional<vec2<int> > resolution_q_limits(int p) {
      // Find the resolution limits. Set up the quadratic to solve
      double a = cp_beg_[6];
      double b = 2.0 * p * cp_beg_[5];
      double c = p * p * cp_beg_[3] + cp_beg_[0] * dstarmax2_;
      af::small<double, 2> limits = reeke_detail::solve_quad(a, b, c);
      if (limits.size() == 0) {
        return boost::optional<vec2<int> >();
      }

      // Extend limits by the margin, ensuring there is a range even for
      // a single quadratic root
      return vec2<int>((int)af::min(limits.const_ref()) - margin_,
                       (int)af::max(limits.const_ref()) + margin_);
    }

    /**
     * @returns The ewald sphere q limits
     */
    boost::optional<vec2<int> > ewald_sphere_q_limits(int p) {
      // Ewald sphere limits for the beginning setting
      double a = cp_beg_[6];
      double b = 2.0 * (cp_beg_[4] + p * cp_beg_[5]);
      double c = cp_beg_[1] + p * (2 * cp_beg_[2] + p * cp_beg_[3]);
      af::small<double, 2> limits_beg = reeke_detail::solve_quad(a, b, c);

      // Ewald sphere limits for the end setting
      a = cp_end_[6];
      b = 2.0 * (cp_end_[4] + p * cp_end_[5]);
      c = cp_end_[1] + p * (2 * cp_end_[2] + p * cp_end_[3]);
      af::small<double, 2> limits_end = reeke_detail::solve_quad(a, b, c);

      // Determine the overall Ewald limits
      af::small<double, 4> limits;
      for (std::size_t i = 0; i < limits_beg.size(); ++i) {
        limits.push_back(limits_beg[i]);
      }
      for (std::size_t i = 0; i < limits_end.size(); ++i) {
        limits.push_back(limits_end[i]);
      }
      if (limits.size() == 0) {
        return boost::optional<vec2<int> >();
      }

      // Return the limits
      return vec2<int>((int)af::min(limits.const_ref()) - margin_,
                       (int)af::max(limits.const_ref()) + margin_);
    }

    /**
     * @returns The resolution r limits
     */
    boost::optional<vec2<int> > resolution_r_limits(double p,
                                                    double q,
                                                    af::tiny<double, 4> cq) {
      // First the resolution limits. Set up the quadratic to solve
      double a = cp_beg_[0];
      double b = cq[0] + q * cp_beg_[8];
      double c = cq[1] + q * q * cp_beg_[10] + q * cq[2] - dstarmax2_;
      af::small<double, 2> limits = reeke_detail::solve_quad(a, b, c);
      if (limits.size() == 0) {
        return boost::optional<vec2<int> >();
      }

      // Extend limits by the margin, ensuring there is a range even for
      // a single quadratic root
      return vec2<int>((int)af::min(limits.const_ref()) - margin_,
                       (int)af::max(limits.const_ref()) + margin_);
    }

    /**
     * @returns The ewald sphere r limits
     */
    af::small<vec2<int>, 2> ewald_sphere_r_limits(double p,
                                                  double q,
                                                  af::tiny<double, 4> cq_beg,
                                                  af::tiny<double, 4> cq_end) {
      // Ewald sphere limits for the beginning setting
      double a = cp_beg_[0];
      double b = cq_beg[0] + q * cp_beg_[8] + cp_beg_[12];
      double c =
        cq_beg[1] + q * (cq_beg[2] + cp_beg_[13]) + q * q * cp_beg_[10] + cq_beg[3];
      af::small<double, 2> limits_beg = reeke_detail::solve_quad(a, b, c);

      // Ewald sphere limits for the end setting
      a = cp_end_[0];
      b = cq_end[0] + q * cp_end_[8] + cp_end_[12];
      c = cq_end[1] + q * (cq_end[2] + cp_end_[13]) + q * q * cp_end_[10] + cq_end[3];
      af::small<double, 2> limits_end = reeke_detail::solve_quad(a, b, c);

      // Set up two loops, one for each range swept out by a point of
      // intersection as it travels from the beginning to the end setting.
      af::small<vec2<int>, 2> result;
      if (limits_beg.size() > 0 && limits_end.size() > 0) {
        double min_beg = af::min(limits_beg.const_ref());
        double max_beg = af::max(limits_beg.const_ref());
        double min_end = af::min(limits_end.const_ref());
        double max_end = af::max(limits_end.const_ref());
        result.push_back(vec2<int>((int)std::min(min_beg, min_end) - margin_,
                                   (int)std::max(min_beg, min_end) + margin_));
        result.push_back(vec2<int>((int)std::min(max_beg, max_end) - margin_,
                                   (int)std::max(max_beg, max_end) + margin_));
      } else {
        if (limits_beg.size() > 0) {
          result.push_back(vec2<int>((int)af::min(limits_beg.const_ref()) - margin_,
                                     (int)af::max(limits_beg.const_ref()) + margin_));
        }
        if (limits_end.size() > 0) {
          result.push_back(vec2<int>((int)af::min(limits_end.const_ref()) - margin_,
                                     (int)af::max(limits_end.const_ref()) + margin_));
        }
      }
      return result;
    }

    mat3<std::size_t> permutation_;
    vec3<double> source_beg_, source_end_, axis_;
    af::small<double, 15> cp_beg_;
    af::small<double, 15> cp_end_;
    vec2<int> p_lim_;
    double dstarmax_, dstarmax2_;
    double wavelength_beg_, wavelength_sq_beg_, wavelength_end_, wavelength_sq_end_;
    int margin_;
    vec2<double> ewald_p_lim_beg_;
    vec2<double> ewald_p_lim_end_;
    vec2<double> res_p_lim_beg_;
    vec2<double> res_p_lim_end_;
  };

  /**
   * A class to generate miller indices with the reeke algorithm
   */
  class ReekeIndexGenerator {
  public:
    /**
     * Initialise the reeke model with single s0 vector
     * @param ub_beg The starting UB matrix
     * @param ub_end The ending UB matrix
     * @param axis The rotation axis
     * @param s0 The incident beam vector
     * @param dmin The resolution limit
     * @param margin The additional margin to add around limits
     */
    ReekeIndexGenerator(mat3<double> ub_beg,
                        mat3<double> ub_end,
                        cctbx::sgtbx::space_group_type const& space_group_type,
                        vec3<double> axis,
                        vec3<double> s0,
                        double dmin,
                        int margin)
        : model_(ub_beg, ub_end, axis, -s0, -s0, dmin, margin),
          space_group_type_(space_group_type) {}

    /**
     * Initialise the reeke model with varying s0 vector
     * @param ub_beg The starting UB matrix
     * @param ub_end The ending UB matrix
     * @param axis The rotation axis
     * @param s0_beg The starting beam vector
     * @param s0_end The ending beam vector
     * @param dmin The resolution limit
     * @param margin The additional margin to add around limits
     */
    ReekeIndexGenerator(mat3<double> ub_beg,
                        mat3<double> ub_end,
                        cctbx::sgtbx::space_group_type const& space_group_type,
                        vec3<double> axis,
                        vec3<double> s0_beg,
                        vec3<double> s0_end,
                        double dmin,
                        int margin)
        : model_(ub_beg, ub_end, axis, -s0_beg, -s0_end, dmin, margin),
          space_group_type_(space_group_type) {}

    /**
     * @returns The next miller index to be generated
     */
    cctbx::miller::index<> next() {
      for (;;) {
        cctbx::miller::index<> h = model_.permutation() * next_pqr();
        if (h.is_zero()) {
          break;
        }
        if (!space_group_type_.group().is_sys_absent(h)) {
          return h;
        }
      }
      return cctbx::miller::index<>(0, 0, 0);
    }

    /**
     * Create an array of miller indices by calling next until (000) is reached.
     * @returns The array of valid miller indices
     */
    af::shared<cctbx::miller::index<> > to_array() {
      af::shared<cctbx::miller::index<> > result;
      for (;;) {
        cctbx::miller::index<> h = this->next();
        if (h.is_zero()) {
          break;
        }
        result.push_back(h);
      }
      return result;
    }

  private:
    /**
     * @returns The next pqr index
     */
    cctbx::miller::index<> next_pqr() {
      // Constants to make clearer
      const int enter = 0;
      const int yield = 1;

      // Static variables
      static vec2<int> p;
      static vec2<int> q;
      static af::small<vec2<int>, 2> r;
      static std::size_t ridx;
      static int state = enter;

      // This switch simulates a co-routine or python generator. The first time
      // the function is executed, control starts at the top (case 1). On
      // subsequent calls, control starts after the "yield" point (case 0), The
      // static variables ensure that the state is recovered on each subsequent
      // function call.
      cctbx::miller::index<> result(0, 0, 0);
      switch (state) {
      case enter:
        state = yield;
        p = model_.p_limits();
        for (; p[0] < p[1]; ++p[0]) {
          q = model_.q_limits(p[0]);
          for (; q[0] < q[1]; ++q[0]) {
            r = model_.r_limits(p[0], q[0]);
            ridx = 0;
            for (; ridx < r.size(); ++ridx) {
              for (; r[ridx][0] < r[ridx][1]; ++r[ridx][0]) {
                result = cctbx::miller::index<>(p[0], q[0], r[ridx][0]);
                if (!result.is_zero()) {
                  return result;
                case yield:;
                }
              }
            }
          }
        }
      }
      state = enter;
      return cctbx::miller::index<>(0, 0, 0);
    }

    ReekeModel model_;
    cctbx::sgtbx::space_group_type space_group_type_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_REEKE_INDEX_GENERATOR_H
