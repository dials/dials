#include <x2tbx.h>

namespace x2tbx {

  ObservationList::ObservationList()
  {
    imean_ = 0.0;
    sigimean_ = 0.0;
  }

  /**
   * Add one observation
   * @param o I and sigI of observation
   */

  void
  ObservationList::add(i_sig_type o)
  {
    observations_.push_back(o);
  }

  /**
   * Merge observations - necessary after all reflections loaded.
   */

  void
  ObservationList::merge()
  {
    CCTBX_ASSERT(observations_.size() > 0);
    float sum_wi = 0.0;
    float sum_w = 0.0;
    total_i_sigi_ = 0.0;

    for(size_t j = 0; j < observations_.size(); j ++) {
      float i = observations_[j][0];
      float w = 1.0 / (observations_[j][1] * observations_[j][1]);
      sum_w += w;
      sum_wi += w * i;

      total_i_sigi_ += observations_[j][0] / observations_[j][1];
    }

    imean_ = sum_wi / sum_w;
    sigimean_ = 1.0 / sqrt(sum_w);
  }

  /**
   * Get average I/sigma for observations
   * @returns i_sig_type containing Imean, sigImean
   */

  i_sig_type
  ObservationList::i_sigma()
  {
    return i_sig_type(imean_, sigimean_);
  }

  /**
   * Get total (unweighted) I/sigma for observations
   * @returns total I/sigma - useful for computing average unmerged I/sigma
   */

  float
  ObservationList::total_i_sigma()
  {
    return total_i_sigi_;
  }

  /**
   * Get multiplicity of observations of this reflection.
   * @returns multiplicity
   */

  size_t
  ObservationList::multiplicity()
  {
    return observations_.size();
  }

  /**
   * Get rmerge contribution from observations of this reflection, equal to
   * sum(|I - <I>|) for all observations I and weighted mean <I>.
   * @returns rmerge contribution
   */

  float
  ObservationList::rmerge()
  {
    CCTBX_ASSERT(observations_.size() > 0);
    CCTBX_ASSERT(sigimean_ > 0.0);
    float sum_di = 0.0;

    for(size_t j = 0; j < observations_.size(); j ++) {
      sum_di += fabs(observations_[j][0] - imean_);
    }

    return sum_di;
  }
}
