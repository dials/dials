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
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/integration/summation.h>

namespace dials { namespace algorithms {

  using scitbx::af::int4;
  using scitbx::af::int6;
  using model::Image;
  using model::Valid;

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

    /** Arbitrary constant to limit size of measurement box. */
    static const std::size_t MAX_BOX_SIZE = 100;

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
        panel_(panel.begin(), panel.end()),
        frame_(frame.begin(), frame.end()),
        bbox_(bbox.begin(), bbox.end()),
        offset_(last-first+1),
        intensity_(panel.size()),
        variance_(panel.size()) {

      // Check input is valid
      DIALS_ASSERT(first_ < last_);
      DIALS_ASSERT(panel.size() > 0);
      DIALS_ASSERT(panel.size() == frame.size());
      DIALS_ASSERT(panel.size() == bbox.size());

      // Check the frame numbers are in order and within the bounds.  Also
      // create an array of offsets into the data for each frame.  Also compute
      // the maximum buffer size needed to store the shoeboxes for each image.
      std::size_t max_buffer_size = 0;
      std::size_t i = 0;
      offset_[0] = 0;
      DIALS_ASSERT(frame.front() >= first_);
      DIALS_ASSERT(frame.back() < last_);
      DIALS_ASSERT(frame.front() <= frame.back());
      for (std::size_t f = first_; f < last_; ++f) {
        while (i < frame.size() && frame[i] == f) {
          DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
          DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
          std::size_t xs = bbox[i][1] - bbox[i][0];
          std::size_t ys = bbox[i][3] - bbox[i][2];
          DIALS_ASSERT(xs < MAX_BOX_SIZE);
          DIALS_ASSERT(ys < MAX_BOX_SIZE);
          max_buffer_size = std::max(max_buffer_size, xs * ys);
          i++;
        }
        offset_[f - first_ + 1] = i;
      }
      DIALS_ASSERT(i == frame.size());
      DIALS_ASSERT(offset_.back() == frame.size());

      // Allocate the buffers
      data_.resize(max_buffer_size);
      mask_.resize(max_buffer_size);
      bgrd_.resize(max_buffer_size);
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

      // Loop through all the reflections
      std::size_t index = current_ - first_;
      for (std::size_t ind = offset_[index]; ind < offset_[index+1]; ++ind) {
        DIALS_ASSERT(frame_[ind] == current_);

        // Get the panel and bounding box
        std::size_t panel = panel_[ind];
        int4 bbox = bbox_[ind];

        // Compute the grid size
        af::c_grid<2> grid(bbox[3] - bbox[2], bbox[1] - bbox[0]);
        DIALS_ASSERT(grid.size_1d() <= data_.size());

        // Create the data, background and mask arrays
        af::ref< double, af::c_grid<2> > data(&data_[0], grid);
        af::ref< double, af::c_grid<2> > bgrd(&bgrd_[0], grid);
        af::ref< int, af::c_grid<2> > mask(&mask_[0], grid);

        // Get the data for the panel
        Image::int_const_ref_type imdata = image.data(panel);
        Image::bool_const_ref_type immask = image.mask(panel);

        // Copy the image and data
        int x0 = std::max(bbox[0], 0);
        int y0 = std::max(bbox[2], 0);
        int x1 = std::min(bbox[1], (int)imdata.accessor()[1]);
        int y1 = std::min(bbox[3], (int)imdata.accessor()[0]);
        std::fill(data.begin(), data.end(), 0);
        std::fill(bgrd.begin(), bgrd.end(), 0);
        std::fill(mask.begin(), mask.end(), 0);
        for (int  y = y0; y < y1; ++y) {
          for (int x = x0; x < x1; ++x) {
            int j = y - bbox[2];
            int i = x - bbox[0];
            data(j,i) = imdata(y,x);
            mask(j,i) = immask(y,x) ? Valid : 0;
          }
        }

        // Do the summation integration
        Summation<> summation(data, bgrd, mask);
        intensity_[ind] = summation.intensity();
        variance_[ind] = summation.variance();
      }

      // Increment the frame number
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
    af::shared<std::size_t> offset_;
    af::shared<double> intensity_;
    af::shared<double> variance_;
    af::shared<double> data_;
    af::shared<double> bgrd_;
    af::shared<int> mask_;
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

    /**
     * Initialise the book-keeping for the integration.
     * @param predicted The predicted reflections.
     * @param image0 The first image
     * @param image1 The last image
     * @param nproc The number of processors
     */
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
