#ifndef DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_3D_H
#define DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_3D_H

#include <cmath>      // for exp, pow
#include <algorithm>  // for std::min, std::max
#include <scitbx/vec2.h>
#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <boost/math/special_functions/round.hpp>  // for iround
#include "gaussian_smoother.h"

namespace dials { namespace refinement {

  using scitbx::vec2;
  using scitbx::sparse::matrix;
  using scitbx::sparse::vector;

  // A Gaussian smoother, based largely on class SmoothedValue from Aimless.
  class GaussianSmoother3D {
  public:
    // Construct from range of raw unnormalised coordinate & number of sample
    // intervals. Set smoothing values to defaults, naverage = 3
    GaussianSmoother3D(vec2<double> x_range,
                       std::size_t num_x_intervals,
                       vec2<double> y_range,
                       std::size_t num_y_intervals,
                       vec2<double> z_range,
                       std::size_t num_z_intervals)
        : x0(x_range[0]),  // initialise members
          y0(y_range[0]),
          z0(z_range[0]),
          nxsample(num_x_intervals),
          nysample(num_y_intervals),
          nzsample(num_z_intervals) {
      DIALS_ASSERT(nxsample > 0);
      if (nxsample == 1) {
        nxvalues = 2;
      } else if (nxsample == 2) {
        nxvalues = 3;
      } else {
        nxvalues = nxsample + 2;
      }

      DIALS_ASSERT(nysample > 0);
      if (nysample == 1) {
        nyvalues = 2;
      } else if (nysample == 2) {
        nyvalues = 3;
      } else {
        nyvalues = nysample + 2;
      }

      DIALS_ASSERT(nzsample > 0);
      if (nzsample == 1) {
        nzvalues = 2;
      } else if (nzsample == 2) {
        nzvalues = 3;
      } else {
        nzvalues = nzsample + 2;
      }

      // smoothing spacing
      x_spacing_ = (x_range[1] - x_range[0]) / (double)nxsample;
      y_spacing_ = (y_range[1] - y_range[0]) / (double)nysample;
      z_spacing_ = (z_range[1] - z_range[0]) / (double)nzsample;

      // positions of the smoother parameters
      if (nxvalues < 4) {
        for (std::size_t i = 0; i < nxvalues; ++i) {
          x_positions_.push_back((double)i);
        }
      } else {
        for (std::size_t i = 0; i < nxvalues; ++i) {
          x_positions_.push_back((double)i - 0.5);
        }
      }

      if (nyvalues < 4) {
        for (std::size_t i = 0; i < nyvalues; ++i) {
          y_positions_.push_back((double)i);
        }
      } else {
        for (std::size_t i = 0; i < nyvalues; ++i) {
          y_positions_.push_back((double)i - 0.5);
        }
      }

      if (nzvalues < 4) {
        for (std::size_t i = 0; i < nzvalues; ++i) {
          z_positions_.push_back((double)i);
        }
      } else {
        for (std::size_t i = 0; i < nzvalues; ++i) {
          z_positions_.push_back((double)i - 0.5);
        }
      }

      // set default smoothing parameters - keep same as 1D for now.
      set_smoothing(3, -1.0);
    }

    /**
     * Set smoothing parameters. If sigma < 0, set to "optimum" (!) (or at
     * least suitable) value from num_average
     * @param num_average The number of points included in each calculation
     * @param sigma The width of the Gaussian used for smoothing
     */
    void set_smoothing(std::size_t num_average, double sigma) {
      // naverage cannot be greater than the number of values
      n_x_average = std::min(num_average, nxvalues);
      n_y_average = std::min(num_average, nyvalues);
      n_z_average = std::min(num_average, nzvalues);

      // In addition, like Aimless, limit it to the range [1,5]
      if (n_x_average < 1 || n_x_average > 5) {
        throw DIALS_ERROR("GaussianSmoother3D:: n_x_average must be between 1 & 5");
      }
      if (n_y_average < 1 || n_y_average > 5) {
        throw DIALS_ERROR("GaussianSmoother3D:: n_y_average must be between 1 & 5");
      }
      if (n_z_average < 1 || n_z_average > 5) {
        throw DIALS_ERROR("GaussianSmoother3D:: n_z_average must be between 1 & 5");
      }

      // sigma cannot be set to zero
      if (sigma == 0.0) {
        throw DIALS_ERROR("GaussianSmoother3D:: sigma cannot be set equal to zero");
      }

      half_nxaverage = (double)n_x_average / 2.0;
      half_nyaverage = (double)n_y_average / 2.0;
      half_nzaverage = (double)n_z_average / 2.0;
      std::size_t naverage = std::min(n_x_average, n_y_average);
      std::size_t min_naverage = std::min(naverage, n_z_average);
      sigma_ = sigma;

      if (sigma_ < 0.0) {
        //  Default values 0.65, 0.7, 0.75, 0.8 for nav = 2,3,4,5
        sigma_ = 0.65 + 0.05 * (min_naverage - 2);
      }
    }

    // Get number of values
    std::size_t num_x_values() {
      return nxvalues;
    }

    std::size_t num_y_values() {
      return nyvalues;
    }

    std::size_t num_z_values() {
      return nzvalues;
    }

    // Get number of sample intervals
    std::size_t num_x_samples() {
      return nxsample;
    }

    std::size_t num_y_samples() {
      return nysample;
    }

    std::size_t num_z_samples() {
      return nzsample;
    }

    // Get number of points averaged in each dimension
    std::size_t num_x_average() {
      return n_x_average;
    }

    std::size_t num_y_average() {
      return n_y_average;
    }

    std::size_t num_z_average() {
      return n_z_average;
    }

    // Get sigma smoothing factor
    double sigma() {
      return sigma_;
    }

    // Get spacing
    double x_spacing() {
      return x_spacing_;
    }

    double y_spacing() {
      return y_spacing_;
    }

    double z_spacing() {
      return z_spacing_;
    }

    // Get positions
    af::shared<double> x_positions() const {
      return x_positions_;
    }

    af::shared<double> y_positions() const {
      return y_positions_;
    }

    af::shared<double> z_positions() const {
      return z_positions_;
    }
    /**
     * Calculate a single interpolated value of param at a point using the
     * original unnormalised coordinate. Return this along with the weights at
     * each position, and their sum.
     * @param x The point to interpolate at
     * @param values The parameter values
     */
    SingleValueWeights value_weight(double x,
                                    double y,
                                    double z,
                                    const af::const_ref<double> values) {
      // use sparse storage as only naverage (default 3) values are non-zero
      vector<double> weight(nxvalues * nyvalues * nzvalues);

      // normalised coordinate
      double z1 = (x - x0) / x_spacing_;
      double z2 = (y - y0) / y_spacing_;
      double z3 = (z - z0) / z_spacing_;
      double sumwv = 0.0;
      double sumweight = 0.0;

      vec2<int> irange = idx_range(z1, nxvalues, half_nxaverage, n_x_average);
      vec2<int> jrange = idx_range(z2, nyvalues, half_nyaverage, n_y_average);
      vec2<int> krange = idx_range(z3, nzvalues, half_nzaverage, n_z_average);

      for (int i = irange[0]; i < irange[1]; ++i) {
        for (int j = jrange[0]; j < jrange[1]; ++j) {
          for (int k = krange[0]; k < krange[1]; ++k) {
            double ds = pow(pow(z1 - x_positions_[i], 2) + pow(z2 - y_positions_[j], 2)
                              + pow(z3 - z_positions_[k], 2),
                            0.5)
                        / sigma_;
            int idx = i + (j * nxvalues) + (k * nxvalues * nyvalues);
            weight[idx] = exp(-ds * ds);
            sumwv += weight[idx] * values[idx];
            sumweight += weight[idx];
          }
        }
      }

      double value;
      if (sumweight > 0.0) {
        value = sumwv / sumweight;
      } else {
        value = 0.0;
      }

      return SingleValueWeights(value, weight, sumweight);
    }

    /**
     * Calculate multiple interpolated values of the parameters at points using
     * the original unnormalised coordinate. Return this along with the matrix
     * of weights at each position, and their sums. The matrix of weights is
     * arranged such that each row contains the weights for a single point
     * taken from x.
     * @param x The array of points to interpolate at
     * @param param The parameter values
     */
    MultiValueWeights multi_value_weight(const af::const_ref<double> x,
                                         const af::const_ref<double> y,
                                         const af::const_ref<double> z,
                                         const af::const_ref<double> values) {
      // Use sparse storage as only naverage (default 3) values per row are
      // non-zero
      std::size_t nxpoints = x.size();
      std::size_t nypoints = y.size();
      std::size_t nzpoints = z.size();
      DIALS_ASSERT(nxpoints > 1);
      DIALS_ASSERT(nypoints == nxpoints);
      DIALS_ASSERT(nzpoints == nxpoints);
      matrix<double> weight(nxpoints, nxvalues * nyvalues * nzvalues);

      // Allocate space for the interpolated values and sumweights, with raw
      // refs for fastest access (See Michael Hohn's notes)
      af::shared<double> value(nxpoints, af::init_functor_null<double>());
      af::ref<double> value_ref = value.ref();
      af::shared<double> sumweight(nxpoints, af::init_functor_null<double>());
      af::ref<double> sumweight_ref = sumweight.ref();

      for (std::size_t irow = 0; irow < nxpoints; ++irow) {
        // normalised coordinate
        double z1 = (x[irow] - x0) / x_spacing_;
        double z2 = (y[irow] - y0) / y_spacing_;
        double z3 = (z[irow] - z0) / z_spacing_;
        double sumw = 0.0;
        double sumwv = 0.0;

        vec2<int> irange = idx_range(z1, nxvalues, half_nxaverage, n_x_average);
        vec2<int> jrange = idx_range(z2, nyvalues, half_nyaverage, n_y_average);
        vec2<int> krange = idx_range(z3, nzvalues, half_nzaverage, n_z_average);

        for (int icol = irange[0]; icol < irange[1]; ++icol) {
          for (int jcol = jrange[0]; jcol < jrange[1]; ++jcol) {
            for (int kcol = krange[0]; kcol < krange[1]; ++kcol) {
              double ds =
                pow(pow(z1 - x_positions_[icol], 2) + pow(z2 - y_positions_[jcol], 2)
                      + pow(z3 - z_positions_[kcol], 2),
                    0.5)
                / sigma_;
              double w = exp(-ds * ds);
              int idx = icol + (jcol * nxvalues) + (kcol * nxvalues * nyvalues);
              weight(irow, idx) = w;
              sumw += w;
              sumwv += w * values[idx];
            }
          }
        }
        sumweight_ref[irow] = sumw;

        if (sumw > 0.0) {
          value_ref[irow] = sumwv / sumw;
        } else {
          value_ref[irow] = 0.0;
        }
      }

      return MultiValueWeights(value, weight, sumweight);
    }

  private:
    vec2<int> idx_range(double z,
                        std::size_t nvalues,
                        double half_naverage,
                        std::size_t naverage) {
      int i1, i2;
      if (nvalues <= 3) {
        i1 = 0;
        i2 = nvalues;
      } else {  // in this case, 1st point in array (index 0) is at position -0.5
        i1 = boost::math::iround(z - half_naverage) + 1;
        i2 = i1 + naverage;
        if (i1 < 0) {  // set to beginning of range
          i1 = 0;
          i2 = std::max(2, i2);  // ensure a separation of at least 2
        }
        if (i2 > nvalues) {
          i2 = nvalues;
          i1 = std::min(i1, (int)nvalues - 2);  // ensure separation of >= 2
        }
      }
      return vec2<int>(i1, i2);
    }

    double x0;
    double y0;
    double z0;
    double x_spacing_;
    double y_spacing_;
    double z_spacing_;
    double half_nxaverage;
    double half_nyaverage;
    double half_nzaverage;
    double sigma_;  // keep same?
    std::size_t nxsample;
    std::size_t nysample;
    std::size_t nzsample;
    std::size_t nxvalues;
    std::size_t nyvalues;
    std::size_t nzvalues;
    std::size_t n_x_average;
    std::size_t n_y_average;
    std::size_t n_z_average;
    af::shared<double> x_positions_;
    af::shared<double> y_positions_;
    af::shared<double> z_positions_;
  };

}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_3D_H
