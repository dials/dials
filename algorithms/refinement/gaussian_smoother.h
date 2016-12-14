/*
 * gaussian_smoother.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_H
#define DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_H

#include <cmath> // for exp
#include <algorithm> // for std::min, std::max
#include <scitbx/vec2.h>
#include <scitbx/sparse/vector.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <boost/math/special_functions/round.hpp> // for iround

namespace dials { namespace refinement {

  using scitbx::vec2;
  using scitbx::sparse::vector;

  // A return type for GaussianSmoother::value_weight
  struct SingleValueWeights {
    double value;
    vector<double> weight;
    double sumweight;

    SingleValueWeights(double value_, vector<double> weight_, double sumweight_)
      : value(value_),
        weight(weight_),
        sumweight(sumweight_){}

    double get_value() const {
      return value;
    }

    vector<double> get_weight() const {
      return weight;
    }

    double get_sumweight() const {
      return sumweight;
    }
  };

  // A Gaussian smoother, based largely on class SmoothedValue from Aimless.
  class GaussianSmoother {
  public:

    // Construct from range of raw unnormalised coordinate & number of sample
    // intervals. Set smoothing values to defaults, naverage = 3
    GaussianSmoother(vec2<double> x_range,
                     std::size_t num_intervals)
      : x0(x_range[0]), // initialisation lists
        nsample(num_intervals){

      DIALS_ASSERT(nsample > 0);
      if (nsample == 1){
        nvalues = 2;
      } else if (nsample == 2) {
        nvalues = 3;
      } else {
        nvalues = nsample + 2;
      }

      // smoothing spacing
      spacing_ = (x_range[1] - x_range[0]) / (double)nsample;

      // positions of the smoother parameters
      if (nvalues < 4){
        for (std::size_t i = 0; i < nvalues; ++i) {
          positions_.push_back((double)i);
        }
      } else {
        for (std::size_t i = 0; i < nvalues; ++i) {
          positions_.push_back((double)i - 0.5);
        }
      }

      // set default smoothing parameters
      set_smoothing(3 , -1.0);
    }

    /**
     * Set smoothing parameters. If sigma < 0, set to "optimum" (!) (or at
     * least suitable) value from num_average
     * @param num_average The number of points included in each calculation
     * @param sigma The width of the Gaussian used for smoothing
     */
    void set_smoothing(std::size_t num_average, double sigma) {

      // naverage cannot be greater than the number of values
      naverage = std::min(num_average, nvalues);

      // In addition, like Aimless, limit it to the range [1,5]
      if (naverage < 1 || naverage > 5) {
        DIALS_ERROR("GaussianSmoother:: num_average must be between 1 & 5");
      }

      // sigma cannot be set to zero
      if (sigma == 0.0) {
        DIALS_ERROR("GaussianSmoother:: sigma cannot be set equal to zero");
      }

      half_naverage = (double)naverage / 2.0;
      sigma_ = sigma;

      if (sigma_ < 0.0) {
        //  Default values 0.65, 0.7, 0.75, 0.8 for nav = 2,3,4,5
          sigma_ = 0.65 + 0.05 * (naverage - 2);
      }
    }

    // Get number of values
    std::size_t num_values(){
      return nvalues;
    }

    // Get number of sample intervals
    std::size_t num_samples(){
      return nsample;
    }

    // Get number of points averaged
    std::size_t num_average(){
      return naverage;
    }

    // Get sigma smoothing factor
    double sigma(){
      return sigma_;
    }

    // Get spacing
    double spacing(){
      return spacing_;
    }

    // Get positions
    af::shared<double> positions() const {
      return positions_;
    }

    /**
     * Calculate a single interpolated value of param at a point using the
     * original unnormalised coordinate. Return this along with the weights at
     * each position, and their sum.
     * @param x The point to interpolate at
     * @param values The parameter values
     */
    SingleValueWeights value_weight(double x,
                                    const af::const_ref<double> values) {

      // use sparse storage as only naverage (default 3) values are non-zero
      vector<double> weight(nvalues);

      // normalised coordinate
      double z = (x - x0) / spacing_;
      double sumwv = 0.0;
      double sumweight = 0.0;

      int i1, i2;
      if (nvalues <= 3){
        i1 = 0;
        i2 = nvalues;
      } else { // in this case, 1st point in array (index 0) is at position -0.5
        i1 = boost::math::iround(z - half_naverage) + 1;
        i2 = i1 + naverage;
        if (i1 < 0){ // set to beginning of range
          i1 = 0;
          i2 = std::max(2, i2); //ensure a separation of at least 2
        }
        if (i2 > nvalues){
          i2 = nvalues;
          i1 = std::min(i1, (int)nvalues - 2); // ensure separation of >= 2
        }
      }

      for (int i = i1; i < i2; ++i) {
        double ds = (z - positions_[i]) / sigma_;
        weight[i] = exp(-ds*ds);
        sumwv += weight[i] * values[i];
        sumweight  += weight[i];
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
     * of weights at each position, and their sums.
     * @param x The points to interpolate at
     * @param param The parameter values
     */

     // FIXME TODO

  private:

    double x0;
    double spacing_;
    double half_naverage;
    double sigma_;
    std::size_t nsample;
    std::size_t nvalues;
    std::size_t naverage;
    af::shared<double> positions_;

  };

}} // namespace dials::refinement

#endif // DIALS_REFINEMENT_GAUSSIAN_SMOOTHER_H
