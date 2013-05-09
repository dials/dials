/*
 * gaussian_corrected_mean_subtractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_GAUSSIAN_CORRECTED_MEAN_SUBTRACTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_GAUSSIAN_CORRECTED_MEAN_SUBTRACTOR_H

#include <dials/algorithms/background/mean_subtractor.h>

namespace dials { namespace algorithms {

  /**
   * Class for background subtraction by mean background value with
   * a correction for gaussian spots.
   */
  class GaussianCorrectedMeanSubtractor : public MeanSubtractor {
  public:

    /** Initialise the class. */
    GaussianCorrectedMeanSubtractor() {}

    /**
     * Calculate the background intensity
     * @params reflection The reflection
     * @return The background intensity
     */
    virtual int background_intensity(Reflection &reflection) const {
      return MeanSubtractor::background_intensity(reflection);
    }
 };

}}

#endif /* DIALS_ALGORITHMS_BACKGROUND_GAUSSIAN_CORRECTED_MEAN_SUBTRACTOR_H */
