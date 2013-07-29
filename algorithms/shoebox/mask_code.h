/*
 * mask_code.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_MASK_CODE_H
#define DIALS_ALGORITHMS_SHOEBOX_MASK_CODE_H

namespace dials { namespace algorithms { namespace shoebox {

  /**
   * An enumeration of shoebox mask codes. If a pixel is labelled as:
   *  a) Invalid. This means that the pixel either belongs to another
   *     reflection or that the pixel is bad or outside the detector.
   *  b) Background. This means that the pixel is to be used for background
   *     determination.
   *  c) Foreground. This means that the pixel is to be used for integration.
   */
  enum MaskCode {
    Invalid = (1 << 1),     ///< Pixel is invalid for this shoebox
    Valid = (1 << 2),       ///< Pixel is valid for this shoebox
    Background = (1 << 3),  ///< Pixel is in the background
    Foreground = (1 << 4),  ///< Pixel is in the foreground
  };

}}}; // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_CODE_H */
