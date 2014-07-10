


#ifndef DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  class Flattened2DProfileFitting {
  public:

    Flattened2DProfileFitting() {

    }

    af::shared<double> intensity() const {
      return intensity_;
    }

    af::shared<double> variance() const {
      return variance_;
    }

  private:

    af::shared<double> intensity_;
    af::shared<double> variance_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
