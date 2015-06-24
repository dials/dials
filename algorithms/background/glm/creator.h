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

  namespace detail {

    template <typename T>
    T median(const af::const_ref<T> &x) {
      af::shared<T> temp(x.begin(), x.end());
      std::nth_element(temp.begin(), temp.begin() + temp.size() / 2, temp.end());
      return temp[temp.size() / 2];
    }

  }

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
        double median = detail::median(Y.const_ref());
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
      double median = detail::median(Y.const_ref());
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
              X(l,1) = j+0.5;
              X(l,2) = i+0.5;
              l++;
            }
          }
        }
        DIALS_ASSERT(l == Y.size());

        // Check that we have a spread of X/Y
        std::size_t countx = 0;
        std::size_t county = 0;
        for (std::size_t l = 1; l < X.accessor()[0]; ++l) {
          if (X(l,1) != X(l-1,1)) {
            countx++;
          }
          if (X(l,2) != X(l-1,2)) {
            county++;
          }
        }
        DIALS_ASSERT(countx > 0 && county > 0);

        // Compute the median value for the starting value
        double median = detail::median(Y.const_ref());
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
        DIALS_ASSERT(b0 > -300 && b0 < 300);
        DIALS_ASSERT(b1 > -300 && b1 < 300);
        DIALS_ASSERT(b2 > -300 && b2 < 300);

        // Fill in the background shoebox values
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            double eta = b0 + b1 * (j+0.5) + b2*(i+0.5);
            double background = std::exp(eta);
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
      for (std::size_t k = 0; k < sbox.zsize(); ++k) {
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              num_background++;
            }
          }
        }
      }
      DIALS_ASSERT(num_background > 0);

      // Allocate some arrays
      af::versa<double, af::c_grid<2> > X(af::c_grid<2>(num_background,4),0);
      af::shared<double> Y(num_background, 0);
      std::size_t l = 0;
      for (std::size_t k = 0; k < sbox.zsize(); ++k) {
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              DIALS_ASSERT(l < Y.size());
              DIALS_ASSERT(sbox.data(k,j,i) >= 0);
              Y[l] = sbox.data(k,j,i);
              X(l,0) = 1.0;
              X(l,1) = k;
              X(l,2) = j;
              X(l,3) = i;
              l++;
            }
          }
        }
      }
      DIALS_ASSERT(l == Y.size());

      // Check that we have a spread of X/Y
      std::size_t countx = 0;
      std::size_t county = 0;
      std::size_t countz = 0;
      for (std::size_t l = 1; l < X.accessor()[0]; ++l) {
        if (X(l,1) != X(l-1,1)) {
          countx++;
        }
        if (X(l,2) != X(l-1,2)) {
          county++;
        }
        if (X(l,3) != X(l-1,3)) {
          county++;
        }
      }
      DIALS_ASSERT(countx > 0 && county > 0 && countz > 0);

      // Compute the median value for the starting value
      double median = detail::median(Y.const_ref());
      if (median == 0) {
        median = 1.0;
      }

      // Setup the initial parameters
      af::shared<double> B(4);
      B[0] = std::log(median);
      B[1] = 0.0;
      B[2] = 0.0;
      B[3] = 0.0;

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
      DIALS_ASSERT(B.size() == 4);
      double b0 = B[0];
      double b1 = B[1];
      double b2 = B[2];
      double b3 = B[3];
      DIALS_ASSERT(b0 > -300 && b0 < 300);
      DIALS_ASSERT(b1 > -300 && b1 < 300);
      DIALS_ASSERT(b2 > -300 && b2 < 300);
      DIALS_ASSERT(b3 > -300 && b3 < 300);

      // Fill in the background shoebox values
      for (std::size_t k = 0; k < sbox.zsize(); ++k) {
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            double background = std::exp(b0 + b1*k + b2*j + b3*i);
            sbox.background(k,j,i) = background;
            if ((sbox.mask(k,j,i) & mask_code) == mask_code) {
              sbox.mask(k,j,i) |= BackgroundUsed;
            }
          }
        }
      }
    }

    Model model_;
    double tuning_constant_;
    std::size_t max_iter_;

  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
