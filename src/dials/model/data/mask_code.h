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
#ifndef DIALS_MODEL_DATA_MASK_CODE_H
#define DIALS_MODEL_DATA_MASK_CODE_H

namespace dials { namespace model {

  /**
   * An enumeration of shoebox mask codes. If a pixel is labelled as:
   *  a) Valid. This means that the pixel belongs to this reflection and
   *     is not a bad pixel etc.
   *  b) Background. This means that the pixel is to be used for background
   *     determination.
   *  c) Foreground. This means that the pixel is to be used for integration.
   *  d) Strong. This means that the pixel is defined as strong
   */
  enum MaskCode {
    Valid = (1 << 0),           ///< Pixel is valid for this shoebox
    Background = (1 << 1),      ///< Pixel is in the background
    Foreground = (1 << 2),      ///< Pixel is in the foreground
    Strong = (1 << 3),          ///< Pixel is a strong pixel
    BackgroundUsed = (1 << 4),  ///< Pixel was used in background calculation
    Overlapped = (1 << 5),      ///< Pixel overlaps another reflection foreground
  };

}}  // namespace dials::model

#endif
