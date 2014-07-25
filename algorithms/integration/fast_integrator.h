

#ifndef DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H

#include <dials/array_family/reflection_table.h>

namespace dials { namespace algorithms {

  class FastIntegratorResult {
  public:

    FastIntegratorResult() {

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

    FastIntegratorWorker() {

    }

    void next() {
    }

    bool finished() const {
      return true;
    }

    FastIntegratorResult result() {
      return result_;
    }

  private:

    FastIntegratorResult result_;
  };

  class FastIntegrator {
  public:

    FastIntegrator(
          af::reflection_table predicted, 
          std::size_t nproc)
      : result_(predicted),
        accumulated_(nproc, false) {
      DIALS_ASSERT(nproc > 0);
    }

    std::size_t size() const {
      return accumulated_.size();
    }
    
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

    FastIntegratorWorker worker(std::size_t index) const {
      return FastIntegratorWorker();
    }

    void accumulate(const FastIntegratorResult &result) {

      // Ensure the worker index is valid
      DIALS_ASSERT(result.index() < accumulated_.size());

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


    af::reflection_table result() const {
      DIALS_ASSERT(finished());
      return result_;
    }

  private:

    af::reflection_table result_;
    af::shared<bool> accumulated_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FAST_INTEGRATOR_H
