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
     * Enumeration of model types
     */
    enum Model {
      Constant2d,
      Constant3d,
      LogLinear2d,
      LogLinear3d
    };

    /**
     * Initialise the creator
     * @param tuning_constant The robust tuning constant
     * @param max_iter The maximum number of iterations
     */
    Creator(Model model, double tuning_constant, std::size_t max_iter)
      : model_(model),
        tuning_constant_(tuning_constant),
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
      DIALS_ASSERT(sbox.is_consistent());
      switch (model_) {
      case Constant2d:
        compute_constant_2d(sbox);
        break;
      case Constant3d:
        compute_constant_3d(sbox);
        break;
      case LogLinear2d:
        compute_loglinear_2d(sbox);
        break;
      case LogLinear3d:
        compute_loglinear_3d(sbox);
        break;
      default:
        DIALS_ERROR("Unknown Model");
      };
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute_constant_2d(Shoebox<> &sbox) const {

      for (std::size_t k = 0; k < sbox.zsize(); ++k) {

        // Compute number of background pixels
        std::size_t num_background = 0;
        int mask_code = Valid | Background;
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              num_background++;
            }
          }
        }
        DIALS_ASSERT(num_background > 0);

        // Allocate some arrays
        af::shared<double> Y(num_background, 0);
        std::size_t l = 0;
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              DIALS_ASSERT(l < Y.size());
              DIALS_ASSERT(sbox.data(k,j,i) >= 0);
              Y[l++] = sbox.data(k,j,i);
            }
          }
        }
        DIALS_ASSERT(l == Y.size());

        // Compute the median value for the starting value
        std::nth_element(Y.begin(), Y.begin() + Y.size() / 2, Y.end());
        double median = Y[Y.size() / 2];
        if (median == 0) {
          median = 1.0;
        }

        // Compute the result
        RobustPoissonMean result(
            Y.const_ref(),
            median,
            tuning_constant_,
            1e-3,
            max_iter_);
        DIALS_ASSERT(result.converged());

        // Compute the background
        double background = result.mean();

        // Fill in the background shoebox values
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            sbox.background(k,j,i) = background;
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              sbox.mask(k,j,i) |= BackgroundUsed;
            }
          }
        }
      }
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute_constant_3d(Shoebox<> &sbox) const {

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
      std::size_t j = 0;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(j < Y.size());
          DIALS_ASSERT(sbox.data[i] >= 0);
          Y[j++] = sbox.data[i];
        }
      }
      DIALS_ASSERT(j == Y.size());

      // Compute the median value for the starting value
      std::nth_element(Y.begin(), Y.begin() + Y.size() / 2, Y.end());
      double median = Y[Y.size() / 2];
      if (median == 0) {
        median = 1.0;
      }

      // Compute the result
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant_,
          1e-3,
          max_iter_);
      DIALS_ASSERT(result.converged());

      // Compute the background
      double background = result.mean();

      // Fill in the background shoebox values
      for (std::size_t i = 0; i < sbox.background.size(); ++i) {
        sbox.background[i] = background;
        if ((sbox.mask[i] & mask_code) == mask_code) {
          sbox.mask[i] |= BackgroundUsed;
        }
      }
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute_loglinear_2d(Shoebox<> &sbox) const {

      for (std::size_t k = 0; k < sbox.zsize(); ++k) {

        // Compute number of background pixels
        std::size_t num_background = 0;
        int mask_code = Valid | Background;
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              num_background++;
            }
          }
        }
        DIALS_ASSERT(num_background > 0);

        // Allocate some arrays
        af::versa<double, af::c_grid<2> > X(af::c_grid<2>(num_background,3),0);
        af::shared<double> Y(num_background, 0);
        std::size_t l = 0;
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              DIALS_ASSERT(l < Y.size());
              DIALS_ASSERT(sbox.data(k,j,i) >= 0);
              Y[l] = sbox.data(k,j,i);
              X(l,0) = 1.0;
              X(l,1) = j;
              X(l,2) = i;
              l++;
            }
          }
        }
        DIALS_ASSERT(l == Y.size());

        // Compute the median value for the starting value
        std::nth_element(Y.begin(), Y.begin() + Y.size() / 2, Y.end());
        double median = Y[Y.size() / 2];
        if (median == 0) {
          median = 1.0;
        }

        // Setup the initial parameters
        af::shared<double> B(3);
        B[0] = std::log(median);
        B[1] = 0.0;
        B[2] = 0.0;

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
        B = result.parameters();
        DIALS_ASSERT(B.size() == 3);
        double b0 = B[0];
        double b1 = B[1];
        double b2 = B[2];

        // Fill in the background shoebox values
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            double background = std::exp(b0 + b1*j + b2*i);
            sbox.background(k,j,i) = background;
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              sbox.mask(k,j,i) |= BackgroundUsed;
            }
          }
        }
      }
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute_loglinear_3d(Shoebox<> &sbox) const {

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
      std::size_t j = 0;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(j < Y.size());
          DIALS_ASSERT(sbox.data[i] >= 0);
          Y[j++] = sbox.data[i];
        }
      }
      DIALS_ASSERT(j == Y.size());

      // Compute the median value for the starting value
      std::nth_element(Y.begin(), Y.begin() + Y.size() / 2, Y.end());
      double median = Y[Y.size() / 2];
      if (median == 0) {
        median = 1.0;
      }

      // Compute the result
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant_,
          1e-3,
          max_iter_);
      DIALS_ASSERT(result.converged());

      // Compute the background
      double background = result.mean();

      // Fill in the background shoebox values
      for (std::size_t i = 0; i < sbox.background.size(); ++i) {
        sbox.background[i] = background;
        if ((sbox.mask[i] & mask_code) == mask_code) {
          sbox.mask[i] |= BackgroundUsed;
        }
      }
    }

    Model model_;
    double tuning_constant_;
    std::size_t max_iter_;

  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
