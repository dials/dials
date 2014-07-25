

#ifndef DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H

#include <dials/array_family/reflection_table.h>

namespace dials { namespace algorithms {

  class FastIntegratorResult {
  public:

    FastIntegratorResult(std::size_t index) 
      : index_(index) {

    }

    std::size_t index() const {
      return index_;
    }

    af::shared<double> intensity() const {
      return af::shared<double>();
    }

    af::shared<double> variance() const {
      return af::shared<double>();
    }

    af::shared<std::size_t> indices() const {
      return af::shared<std::size_t>();
    }

  private:

    std::size_t index_;
  };

  class FastIntegratorWorker {
  public:

    FastIntegratorWorker(
          std::size_t index, 
          std::size_t first, 
          std::size_t last) 
      : first_(first),
        last_(last),
        current_(first),
        result_(index) {
      DIALS_ASSERT(first_ < last_);
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

    void next() {
      DIALS_ASSERT(!finished());

      std::cout << "Processing Image: " << current_ << std::endl;
      sleep(1);
      current_++;
    }

    bool finished() const {
      return current_ >= last_;
    }

    FastIntegratorResult result() {
      DIALS_ASSERT(finished());
      return result_;
    }

  private:

    std::size_t first_, last_, current_;
    FastIntegratorResult result_;
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

    FastIntegrator(
          af::reflection_table predicted, 
          std::size_t nframes,
          std::size_t nproc)
      : result_(predicted),
        nframes_(nframes),
        nproc_(nproc),
        accumulated_(nproc, false) {
      DIALS_ASSERT(nproc > 0);
      DIALS_ASSERT(nframes > 0);
      DIALS_ASSERT(nproc < nframes);
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
      std::size_t first = 0, last = 0;
      compute_frame_range(index, first, last);
      return FastIntegratorWorker(index, first, last);
    }

    /**
     * Accumulate the results from the partial sums.
     * @param result The results from a single thread
     */
    void accumulate(const FastIntegratorResult &result) {

      // Ensure the worker index is valid
      DIALS_ASSERT(result.index() < accumulated_.size());
      DIALS_ASSERT(accumulated_[result.index()] == false);

      // Get the columns from the table
      af::shared<double> IR = result_["intensity.sum.value"];
      af::shared<double> VR = result_["intensity.sum.variance"];
      DIALS_ASSERT(IR.size() == VR.size());

      // Get the arrays from the worker result
      af::shared<double> IA = result.intensity();
      af::shared<double> VA = result.variance();
      af::shared<std::size_t> indices = result.indices();
      DIALS_ASSERT(IA.size() == VA.size());
      DIALS_ASSERT(IA.size() == indices.size());

      // Accumulate the results
      for (std::size_t i = 0; i < indices.size(); ++i) {
        std::size_t j = indices[i];
        DIALS_ASSERT(j < IR.size());
        IR[j] += IA[i];
        VR[j] += VA[i];
      }

      // Set the accumulated flag
      accumulated_[result.index()] = true;
    }


    /**
     * @returns The results
     */
    af::reflection_table result() const {
      DIALS_ASSERT(finished());
      return result_;
    }

  private:

    /**
     * Compute the frame range.
     * @param index The index to compute for
     * @param first The beginning of the range
     * @param last The end of the range
     */
    void compute_frame_range(std::size_t index, 
        std::size_t &first, std::size_t &last) const {
      double div = (double)nframes_ / (double)nproc_;
      first = (std::size_t)ceil(div * index);
      last = (std::size_t)ceil(div * (index + 1));
      DIALS_ASSERT(first < last);
      if (last < nframes_) {
        DIALS_ERROR("Error computing frames");
      } else if (last == nframes_) {
        // Do nothing
      } else if (last == nframes_ + 1) {
        last -= 1; 
        DIALS_ASSERT(first < last);
      } else {
        DIALS_ERROR("Error computing frames");
      }
    }

    af::reflection_table result_;
    std::size_t nframes_;
    std::size_t nproc_;
    af::shared<bool> accumulated_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
