#ifndef DIALS_REFLECTION_REFLECTION_H
#define DIALS_REFLECTION_REFLECTION_H

#include <map>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/miller.h>

namespace dials {

class ReflectionData {
public:

    ReflectionData() 
        : miller_index_(0, 0, 0),
          rotation_angle_(0.0),
          beam_vector_(0.0, 0.0, 0.0),
          image_coord_(0.0, 0.0, 0.0) {}
          
    ReflectionData(cctbx::miller::index <> miller_index,
                   double rotation_angle,
                   scitbx::vec3 <double> beam_vector,
                   scitbx::vec3 <double> image_coord)
        : miller_index_(miller_index),
          rotation_angle_(rotation_angle),
          beam_vector_(beam_vector),
          image_coord_(image_coord) {}

public:

    cctbx::miller::index <> get_miller_index() const {
        return miller_index_;
    }
    
    double get_rotation_angle() const {
        return rotation_angle_;
    }
    
    scitbx::vec3 <double> get_beam_vector() const {
        return beam_vector_;
    }
    
    scitbx::vec3 <double> get_image_coord() const {
        return image_coord_;
    }

    void set_miller_index(cctbx::miller::index <> miller_index) {
        miller_index_ = miller_index;
    }
    
    void set_rotation_angle(double rotation_angle) {
        rotation_angle_ = rotation_angle;
    }
    
    void set_beam_vector(scitbx::vec3 <double> beam_vector) {
        beam_vector_ = beam_vector;
    }
    
    void set_image_coord(scitbx::vec3 <double>  image_coord) {
        image_coord_ = image_coord;
    }
    
    bool is_zero() {
        return miller_index_.is_zero();
    }

private:

    cctbx::miller::index <> miller_index_;
    double                  rotation_angle_;
    scitbx::vec3 <double>   beam_vector_;
    scitbx::vec3 <double>   image_coord_;
};

class Reflection : public ReflectionData {
public:
    Reflection() : ReflectionData() {}
    
    Reflection(cctbx::miller::index <> miller_index,
               double rotation_angle,
               scitbx::vec3 <double> beam_vector,
               scitbx::vec3 <double> image_coord)
        : ReflectionData(
            miller_index, 
            rotation_angle, 
            beam_vector, 
            image_coord) {}    
};

typedef scitbx::af::shared <Reflection> ReflectionList;
//typedef std::map <cctbx::miller::index <>, Reflection> ReflectionList;
    
} // namespace dials

#endif /* DIALS_REFLECTION_REFLECTION_H */

