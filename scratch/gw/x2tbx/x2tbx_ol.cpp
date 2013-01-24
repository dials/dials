#include <x2tbx.h>

namespace x2tbx {

  ObservationList::ObservationList(void)
  {
    imean = 0.0;
    sigimean = 0.0;
  }

  ObservationList::~ObservationList(void) { }

  void
  ObservationList::add(i_sig_type o)
  {
    observations.push_back(o);
    imean = 0.0;
    sigimean = 0.0;
  }

  void
  ObservationList::merge(void)
  {
    CCTBX_ASSERT(observations.size() > 0);
    float sum_wi = 0.0;
    float sum_w = 0.0;

    for(size_t j = 0; j < observations.size(); j ++) {
      float i = observations[j][0];
      float w = 1.0 / (observations[j][1] * observations[j][1]);
      sum_w += w;
      sum_wi += w * i;
    }

    imean = sum_wi / sum_w;
    sigimean = 1.0 / sqrt(sum_w);
  }

  i_sig_type
  ObservationList::get_i_sigma(void)
  {
    return i_sig_type(imean, sigimean);
  }

  size_t
  ObservationList::multiplicity(void)
  {
    return observations.size();
  }

  float
  ObservationList::rmerge(void)
  {
    CCTBX_ASSERT(observations.size() > 0);
    CCTBX_ASSERT(sigimean > 0.0);
    float sum_di = 0.0;

    for(size_t j = 0; j < observations.size(); j ++) {
      sum_di += fabs(observations[j][0] - imean);
    }

    return sum_di;
  }
}
