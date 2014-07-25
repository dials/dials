/*
 * fast_integrator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H

#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/reflection_table.h>
#include <dials/model/data/image.h>

namespace dials { namespace algorithms {

  using scitbx::af::int4;
  using scitbx::af::int6;
  using model::Image;

  /**
   * A class to contain results from the workers
   */
  class FastIntegratorResult {
  public:

    /**
     * Initialise with the result data.
     */
    FastIntegratorResult(
          std::size_t index,
          const af::const_ref<double> &intensity,
          const af::const_ref<double> &variance)
      : index_(index),
        intensity_(intensity.begin(), intensity.end()),
        variance_(variance.begin(), variance.end()) {
      DIALS_ASSERT(intensity_.size() == variance_.size());
    }

    /**
     * @returns The thread index
     */
    std::size_t index() const {
      return index_;
    }

    /**
     * @returns The intensity values.
     */
    af::shared<double> intensity() const {
      return intensity_;
    }

    /**
     * @returns The variance values.
     */
    af::shared<double> variance() const {
      return variance_;
    }

  private:

    std::size_t index_;
    af::shared<double> intensity_;
    af::shared<double> variance_;
  };

  class FastIntegratorWorker {
  public:

    FastIntegratorWorker(
          std::size_t index, 
          std::size_t first, 
          std::size_t last,
          const af::const_ref<std::size_t> &panel,
          const af::const_ref<std::size_t> &frame,
          const af::const_ref<int4> &bbox) 
      : index_(index),
        first_(first),
        last_(last),
        current_(first),
        intensity_(panel.size()),
        variance_(panel.size()) {
      DIALS_ASSERT(first_ < last_);
      DIALS_ASSERT(panel.size() == frame.size());
      DIALS_ASSERT(panel.size() == bbox.size());
    }

    std::size_t first() const {
      return first_;
    }

    std::size_t last() const {
      return last_;
    }

    std::size_t current() const {
      return current_;
    }

    void next(const Image &image) {
      DIALS_ASSERT(!finished());
      std::cout << "Processing Image: " << current_ << std::endl;
     
      

      current_++;
    }

    bool finished() const {
      return current_ >= last_;
    }

    FastIntegratorResult result() {
      DIALS_ASSERT(finished());
      return FastIntegratorResult(
          index_,
          intensity_.const_ref(),
          variance_.const_ref());
    }

  private:

    std::size_t index_;
    std::size_t first_, last_, current_;
    af::shared<std::size_t> panel_;
    af::shared<std::size_t> frame_;
    af::shared<int4> bbox_;
    af::shared<double> intensity_;
    af::shared<double> variance_;
  };

  /**
   * A class to manage a fast integration of the data. Each image is integrated
   * separately so each slice of a partial reflection is integrated by summation
   * and the results are accumulated at the end of the processing. Since each
   * image is done separately, they can all be done in parallel. To facilitate
   * this, the class returns a worker object which is called for each thread
   * with a sequence of images.
   */
  class FastIntegrator {
  public:

    struct sort_by_frame {
      const std::size_t *frame_;
      sort_by_frame(const std::size_t *frame) : frame_(frame) {}
      bool operator()(std::size_t a, std::size_t b) const {
        return frame_[a] < frame_[b];
      }
    };

    FastIntegrator(
          af::reflection_table predicted, 
          int image0,
          int image1,
          std::size_t nproc)
      : result_(predicted),
        image0_(image0),
        image1_(image1),
        nframes_(image1 - image0),
        nproc_(nproc),
        accumulated_(nproc, false) {
      
      // Check input is ok    
      DIALS_ASSERT(predicted.size() > 0);
      DIALS_ASSERT(predicted.contains("panel"));
      DIALS_ASSERT(predicted.contains("bbox"));
      DIALS_ASSERT(image1 > image0);
      DIALS_ASSERT(nproc > 0);
      DIALS_ASSERT(nframes_ > 0);
      DIALS_ASSERT(nproc < nframes_);

      // Get relavant columns from predicted
      af::shared<std::size_t> panel = predicted["panel"];
      af::shared<int6> bbox = predicted["bbox"];

      // Compute the number of partials
      std::size_t num_partial = 0;
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        int z0 = std::max(image0_, bbox[i][4]);
        int z1 = std::min(image1_, bbox[i][5]);
        DIALS_ASSERT(z1 > z0);
        num_partial += (z1 - z0);
      }

      // initialise the arrays
      bbox2_.resize(num_partial);
      panel_.resize(num_partial);
      frame_.resize(num_partial);
      indices_.resize(num_partial);
      partial_.resize(num_partial);
      for (std::size_t i = 0; i < indices_.size(); ++i) {
        indices_[i] = i;
      }

      // Fill the partial arrays
      std::size_t j1 = 0;
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        int z0 = std::max(image0_, bbox[i][4]);
        int z1 = std::min(image1_, bbox[i][5]);
        for (int z = z0; z < z1; ++z) {
          bbox2_[j1][0] = bbox[i][0];
          bbox2_[j1][1] = bbox[i][1];
          bbox2_[j1][2] = bbox[i][2];
          bbox2_[j1][3] = bbox[i][3];
          panel_[j1] = panel[i];
          frame_[j1] = z - image0_;
          partial_[j1] = i;
          j1++;
        }
      }
      DIALS_ASSERT(j1 == num_partial);

      // Sort the indices by frame 
      std::sort(indices_.begin(), indices_.end(), sort_by_frame(&frame_[0]));

      // Create the frame ranges
      frame_range_.resize(nproc_+1);
      double div = (double)nframes_ / (double)nproc_;
      for (std::size_t i = 0; i < nproc_+1; ++i) {
        frame_range_[i] = (std::size_t)std::floor(div * i + 0.5);
      }
      DIALS_ASSERT(frame_range_.back() == nframes_);

      // Compute the offsets into the index array
      offset_.resize(nproc_+1);
      offset_[0] = 0;
      std::size_t j2 = 0;
      for (std::size_t i = 0; i < nproc_; ++i) {
        std::size_t f = frame_range_[i+1];
        while (j2 < indices_.size() && frame_[indices_[j2]] < f) ++j2;
        offset_[i+1] = j2;
      }
      DIALS_ASSERT(offset_.back() == indices_.size());
    }

    /**
     * @returns The number of workers.
     */
    std::size_t size() const {
      return nproc_;
    }
   
    /**
     * @returns is the processing finished.
     */
    bool finished() const {
      bool f = true;
      for (std::size_t i = 0; i < accumulated_.size(); ++i) {
        if (!accumulated_[i]) {
          f = false;
          break;
        }
      }
      return f;
    }

    /**
     * Construct the worker for a particular thread.
     * @param index The thread index.
     * @returns The worker object
     */
    FastIntegratorWorker worker(std::size_t index) const {
      DIALS_ASSERT(index < size());
  
      // Get the frame range
      std::size_t first = frame_range_[index];
      std::size_t last = frame_range_[index+1];

      // Get the offset and number
      std::size_t off = offset_[index];
      std::size_t num = offset_[index+1] - off;

      // Arrays to pass into worker object
      af::shared<std::size_t> panel(num);
      af::shared<std::size_t> frame(num);
      af::shared<int4> bbox2(num);
      for (std::size_t i = 0; i < num; ++i) {
        std::size_t j = indices_[i+off];
        panel[i] = panel_[j];
        frame[i] = frame_[j];
        bbox2[i] = bbox2_[j];
      }

      // Create the fast integrator worker
      return FastIntegratorWorker(
          index, first, last, 
          panel.const_ref(), 
          frame.const_ref(), 
          bbox2.const_ref());
    }

    /**
     * Accumulate the results from the partial sums.
     * @param result The results from a single thread
     */
    void accumulate(const FastIntegratorResult &result) {

      // Ensure the worker index is valid
      std::size_t index = result.index();
      DIALS_ASSERT(index < size());
      DIALS_ASSERT(accumulated_[index] == false);
      
      // Get the offset and number
      std::size_t off = offset_[index];
      std::size_t num = offset_[index+1] - off;

      // Get the columns from the table
      af::shared<double> IR = result_["intensity.sum.value"];
      af::shared<double> VR = result_["intensity.sum.variance"];
      DIALS_ASSERT(IR.size() == VR.size());

      // Get the arrays from the worker result
      af::shared<double> IA = result.intensity();
      af::shared<double> VA = result.variance();
      DIALS_ASSERT(IA.size() == VA.size());
      DIALS_ASSERT(IA.size() == num);

      // Accumulate the results
      for (std::size_t i = 0; i < num; ++i) {
        std::size_t j = partial_[indices_[i+off]];
        DIALS_ASSERT(j < IR.size());
        IR[j] += IA[i];
        VR[j] += VA[i];
      }

      // Set the accumulated flag
      accumulated_[index] = true;
    }

    /**
     * @returns The results
     */
    af::reflection_table result() const {
      DIALS_ASSERT(finished());
      return result_;
    }

  private:

    af::reflection_table result_;
    int image0_;
    int image1_;
    std::size_t nframes_;
    std::size_t nproc_;
    af::shared<bool> accumulated_;
    af::shared<std::size_t> panel_;
    af::shared<std::size_t> frame_;
    af::shared<int4> bbox2_;
    af::shared<std::size_t> indices_;
    af::shared<std::size_t> partial_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> frame_range_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
