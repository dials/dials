/*
 * interfaces.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_INTERFACES_H
#define DIALS_ALGORITHMS_INTEGRATION_INTERFACES_H

#include <dials/array_family/reflection.h>

namespace dials { namespace algorithms {

  /**
   * Interface class for computing the reflection mask
   */
  class MaskCalculatorIface {
  public:
    virtual ~MaskCalculatorIface() = 0;

    virtual void operator()(af::Reflection &reflection,
                            bool adjacent = false) const = 0;
  };

  // Implementation for pure virtual destructor
  MaskCalculatorIface::~MaskCalculatorIface() {}

  /**
   * Interface class for computing the reflection background
   */
  class BackgroundCalculatorIface {
  public:
    virtual ~BackgroundCalculatorIface() = 0;

    virtual void operator()(af::Reflection &reflection) const = 0;
  };

  // Implementation for pure virtual destructor
  BackgroundCalculatorIface::~BackgroundCalculatorIface() {}

  /**
   * Interface class for computing the reflection intensity
   */
  class IntensityCalculatorIface {
  public:
    virtual ~IntensityCalculatorIface() = 0;

    virtual void operator()(
      af::Reflection &reflection,
      const std::vector<af::Reflection> &adjacent_reflections) const = 0;
  };

  // Implementation for pure virtual destructor
  IntensityCalculatorIface::~IntensityCalculatorIface() {}

  /**
   * Interface class for computing the reference profiles
   */
  class ReferenceCalculatorIface {
  public:
    virtual ~ReferenceCalculatorIface() = 0;

    virtual void operator()(af::Reflection &reflection) = 0;
  };

  // Implementation for pure virtual destructor
  ReferenceCalculatorIface::~ReferenceCalculatorIface() {}

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_INTERFACES_H
