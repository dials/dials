/*
 * goniometer.h
 *
 *   Copyright (C) 2013 Diamond Light Source, James Parkhurst
 *
 *   This code is distributed under the BSD license, a copy of which is
 *   included in the root directory of this package.
 */
#ifndef DIALS_MODEL_GONIOMETER_H
#define DIALS_MODEL_GONIOMETER_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace dials { namespace model {

  using scitbx::vec3;
  using scitbx::mat3;

  /**
   * A class to represent the rotation axis for a standard rotation
   * geometry diffraction data set.
   *
   * The rotation axis assumed to have it's origin at the origin of the
   * laboratory coordinate system. The rotation axis vector is normalized
   * when set using either the constructor or the rotation axis 'setter' 
   *
   * The fixed rotation matrix, F, represents the additional fixed rotation of 
   * the sample attached to the rotating axis - for example the effects of 
   * kappa and phi for an omega scan on a kappa goniometer. The fixed rotation 
   * matrix which should be applied as:
   *
   *  A = [R][F][U][B]
   *
   * in a standard orientation matrix.
   *
   * In the document detailing the conventions used:
   *    rotation_axis -> m2
   *    fixed_rotation_matrix -> *unspecified* 
   */
  class Goniometer {

  public:

    /** 
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix and the rotation axis is set to the x axis
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer() 
      : rotation_axis_(1.0, 0.0, 0.0),
        fixed_rotation_(
          1.0, 0.0, 0.0, 
          0.0, 1.0, 0.0, 
          0.0, 0.0, 1.0) {}

    /** 
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix.
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer(vec3 <double> rotation_axis) 
      : rotation_axis_(rotation_axis.normalize()),
        fixed_rotation_(
          1.0, 0.0, 0.0, 
          0.0, 1.0, 0.0, 
          0.0, 0.0, 1.0) {}

    /** 
     * Initialise the goniometer. 
     * @param rotation_axis The goniometer rotation axis
     * @param fixed_rotation The additional fixed rotation of the sample
     */
    Goniometer(vec3 <double> rotation_axis, 
               mat3 <double> fixed_rotation) 
      : rotation_axis_(rotation_axis.normalize()),
        fixed_rotation_(fixed_rotation) {}

    /** Get the rotation axis */
    vec3 <double> get_rotation_axis() const {
      return rotation_axis_;
    }
    
    /** Get the fixed rotation matrix */
    mat3 <double> get_fixed_rotation_matrix() const {
      return fixed_rotation_;
    }
    
    /** Set the rotation axis */
    void set_rotation_axis(vec3 <double> rotation_axis) {
      rotation_axis_ = rotation_axis.normalize();
    }

    /** Set the fixed rotation matrix */
    void set_fixed_rotation_matrix(mat3 <double> fixed_rotation) {
      fixed_rotation_ = fixed_rotation;
    }
        
  private:

    vec3 <double> rotation_axis_;
    mat3 <double> fixed_rotation_;
  };

}} // namespace dials::model

#endif // DIALS_MODEL_GONIOMETER_H
