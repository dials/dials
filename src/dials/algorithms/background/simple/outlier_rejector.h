/*
 * outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_OUTLIER_REJECTOR_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  class OutlierRejector {
  public:
    virtual ~OutlierRejector() {}

    virtual void mark(const af::const_ref<double, af::c_grid<3> > &data,
                      af::ref<int, af::c_grid<3> > mask) const = 0;
  };

}}}  // namespace dials::algorithms::background

#endif  // DIALS_ALGORITHMS_BACKGROUND_OUTLIER_REJECTOR_H
