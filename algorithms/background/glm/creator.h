/*
 * creator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H

#include <scitbx/glmtbx/robust_glm.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::Shoebox;

  /**
   * A class to create the background model
   */
  class Creator {
  public:

    /**
     * Initialise the creator
     * @param tuning_constant The robust tuning constant
     * @param max_iter The maximum number of iterations
     */
    Creator(double tuning_constant, std::size_t max_iter)
      : tuning_constant_(tuning_constant),
        max_iter_(max_iter) {
      DIALS_ASSERT(tuning_constant > 0);
      DIALS_ASSERT(max_iter > 0);
    }

    /**
     * Compute the background values
     * @param sbox The shoeboxes
     * @returns Success True/False
     */
    af::shared<bool> operator()(af::ref< Shoebox<> > sbox) const {
      af::shared<bool> success(sbox.size(), true);
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        try {
          compute(sbox[i]);
        } catch(scitbx::error) {
          success[i] = false;
        } catch(dials::error) {
          success[i] = false;
        }
      }
      return success;
    }

  private:

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute(Shoebox<> &sbox) const {

      // Check shoebox is ok
      DIALS_ASSERT(sbox.is_consistent());

      // Compute number of background pixels
      std::size_t num_background = 0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & mask_code) == mask_code) {
          num_background++;
        }
      }
      DIALS_ASSERT(num_background > 0);

      // Allocate some arrays
      af::shared<double> Y(num_background, 0);
      af::shared<double> B(1, 0);
      af::versa< double, af::c_grid<2> > X(af::c_grid<2>(num_background, 1), 1);

      // Compute the median value for the starting value
      std::nth_element(Y.begin(), Y.begin() + Y.size() / 2, Y.end());
      double median = Y[Y.size() / 2];
      B[0] = median;

      // Compute the result
      scitbx::glmtbx::robust_glm<scitbx::glmtbx::poisson> result(
          X.const_ref(),
          Y.const_ref(),
          B.const_ref(),
          tuning_constant_,
          1e-3,
          max_iter_);
      DIALS_ASSERT(result.converged());

      // Compute the background
      double background = std::exp(result.parameters()[0]);

      // Fill in the background shoebox values
      for (std::size_t i = 0; i < sbox.background.size(); ++i) {
        sbox.background[i] = background;
      }
    }

    double tuning_constant_;
    std::size_t max_iter_;

  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
