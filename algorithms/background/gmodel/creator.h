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

#include <dials/array_family/reflection_table.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/image_volume.h>
#include <dials/algorithms/background/gmodel/model.h>
#include <dials/algorithms/background/gmodel/robust_estimator.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::Shoebox;
  using model::ImageVolume;
  using model::MultiPanelImageVolume;

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
     * Initialise the creator
     * @param tuning_constant The robust tuning constant
     * @param max_iter The maximum number of iterations
     */
    Creator(
          boost::shared_ptr<BackgroundModel> model,
          bool robust,
          double tuning_constant,
          std::size_t max_iter)
      : model_(model),
        robust_(robust),
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
    af::shared<bool> shoebox(af::reflection_table reflections) const {
      DIALS_ASSERT(reflections.contains("shoebox"));
      af::ref< Shoebox<> > sbox = reflections["shoebox"];
      af::shared<bool> success(sbox.size(), true);
      af::shared<double> scale = reflections["background.scale"];
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        try {
          DIALS_ASSERT(sbox[i].is_consistent());
          scale[i] = compute(
              sbox[i].panel,
              sbox[i].bbox,
              sbox[i].data.const_ref(),
              sbox[i].background.ref(),
              sbox[i].mask.ref());
        } catch(scitbx::error) {
          success[i] = false;
        } catch(dials::error) {
          success[i] = false;
        }
      }
      return success;
    }

    /**
     * Compute the background values
     * @param reflections The reflection table
     * @param volume The image volume
     * @returns Success True/False
     */
    af::shared<bool> volume(
        af::reflection_table reflections,
        MultiPanelImageVolume<> volume) const {
      typedef MultiPanelImageVolume<>::float_type FloatType;

      DIALS_ASSERT(reflections.contains("bbox"));
      DIALS_ASSERT(reflections.contains("panel"));
      af::const_ref< int6 > bbox = reflections["bbox"];
      af::const_ref< std::size_t > panel = reflections["panel"];
      af::shared<bool> success(bbox.size(), true);
      af::shared<double> scale = reflections["background.scale"];
      for (std::size_t i = 0; i < bbox.size(); ++i) {

        // Get the image volume
        ImageVolume<> v = volume.get(panel[i]);

        // Trim the bbox
        int6 b = v.trim_bbox(bbox[i]);
        std::size_t p = panel[i];

        // Extract from image volume
        af::versa< FloatType, af::c_grid<3> > data = v.extract_data(b);
        af::versa< FloatType, af::c_grid<3> > bgrd = v.extract_background(b);
        af::versa< int,    af::c_grid<3> > mask = v.extract_mask(b, i);

        // Compute the background
        try {
          scale[i] = compute(
              p,
              b,
              data.const_ref(),
              bgrd.ref(),
              mask.ref());

          // Need to set the background in volume
          v.set_background(b, bgrd.const_ref());
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
    template <typename T>
    double compute(
        std::size_t panel,
        int6 bbox,
        const af::const_ref< T, af::c_grid<3> > &data,
        af::ref< T, af::c_grid<3> > background,
        af::ref< int, af::c_grid<3> > mask) const {
      if (robust_) {
        return compute_robust(panel, bbox, data, background, mask);
      } else {
        return compute_non_robust(panel, bbox, data, background, mask);
      }
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    template <typename T>
    double compute_non_robust(
        std::size_t panel,
        int6 bbox,
        const af::const_ref< T, af::c_grid<3> > &data,
        af::ref< T, af::c_grid<3> > background,
        af::ref< int, af::c_grid<3> > mask) const {
      af::versa< double, af::c_grid<3> > model = model_->extract(panel, bbox);
      double sum1 = 0;
      double sum2 = 0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          sum1 += data[i];
          sum2 += model[i];
        }
      }
      DIALS_ASSERT(sum1 >= 0);
      DIALS_ASSERT(sum2 > 0);
      double scale = sum1 / sum2;
      for (std::size_t i = 0; i < data.size(); ++i) {
        background[i] = scale * model[i];
        if ((mask[i] & mask_code) == mask_code) {
          mask[i] |= BackgroundUsed;
        }
      }
      return scale;
    }

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    template <typename T>
    double compute_robust(
        std::size_t panel,
        int6 bbox,
        const af::const_ref< T, af::c_grid<3> > &data,
        af::ref< T, af::c_grid<3> > background,
        af::ref< int, af::c_grid<3> > mask) const {
      af::versa< double, af::c_grid<3> > model = model_->extract(panel, bbox);

      // Compute number of background pixels
      std::size_t num_background = 0;
      double sum_model = 0.0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          num_background++;
          sum_model += model[i];
        }
      }
      DIALS_ASSERT(sum_model > 0);
      DIALS_ASSERT(num_background > 0);

      // Allocate some arrays
      af::shared<double> X(num_background, 0);
      af::shared<double> Y(num_background, 0);
      std::size_t l = 0;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(l < Y.size());
          DIALS_ASSERT(data[i] >= 0);
          X[l] = model[i];
          Y[l] = data[i];
          l++;
        }
      }
      DIALS_ASSERT(l == Y.size());

      // Compute the median value for the starting value
      double median = detail::median(Y.const_ref());
      if (median == 0) {
        median = 1.0;
      }

      // Compute robust Poisson mean as in the GLM algorithm
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant_,
          1e-3,
          max_iter_);
      DIALS_ASSERT(result.converged());

      // Get the mean value
      double mean = result.mean();
      double svar = std::sqrt(mean);

      // Compute am initial scale
      double sum1 = 0.0;
      double sum2 = 0.0;
      for (std::size_t i = 0; i < Y.size(); ++i) {
        double res = (Y[i] - mean) / svar;
        double y = Y[i];
        if (res < -tuning_constant_) {
          y = mean - tuning_constant_ * svar;
        } else if (res > tuning_constant_) {
          y = mean + tuning_constant_ * svar;
        }
        sum1 += y;
        sum2 += X[i];
      }
      DIALS_ASSERT(sum1 >= 0);
      DIALS_ASSERT(sum2 > 0);
      double scale = sum1 / sum2;

      // Fill in the background shoebox values
      for (std::size_t i = 0; i < data.size(); ++i) {
        background[i] = scale * model[i];
        if ((mask[i] & mask_code) == mask_code) {
          mask[i] |= BackgroundUsed;
        }
      }

      return scale;
    }

    boost::shared_ptr<BackgroundModel> model_;
    bool robust_;
    double tuning_constant_;
    std::size_t max_iter_;

  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
