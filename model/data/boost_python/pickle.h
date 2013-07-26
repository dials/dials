
#ifndef DIALS_MODEL_DATA_BOOST_PYTHON_PICKLE_H
#define DIALS_MODEL_DATA_BOOST_PYTHON_PICKLE_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/miller.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <dials/error.h>

namespace dials { namespace model { namespace boost_python {
    namespace reflection {

  using scitbx::af::flex_grid;

  struct to_string : scitbx::af::boost_python::pickle_double_buffered::to_string
  {
    using scitbx::af::boost_python::pickle_double_buffered::to_string::operator<<;

    to_string()
    {
      unsigned int version = 1;
      *this << version;
    }

    template <typename ProfileType>
    void profile_to_string(const ProfileType &p) {
      *this << p.accessor().all().size();
      for (std::size_t i = 0; i < p.accessor().all().size(); ++i) {
        *this << p.accessor().all()[i];
      }
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this << p[i];
      }
    }

    to_string& operator<<(const Reflection &val) {
      *this << val.get_miller_index()[0]
            << val.get_miller_index()[1]
            << val.get_miller_index()[2]
            << val.get_status()
            << val.get_entering()
            << val.get_rotation_angle()
            << val.get_beam_vector()[0]
            << val.get_beam_vector()[1]
            << val.get_beam_vector()[2]
            << val.get_image_coord_px()[0]
            << val.get_image_coord_px()[1]
            << val.get_image_coord_mm()[0]
            << val.get_image_coord_mm()[1]
            << val.get_frame_number()
            << val.get_panel_number()
            << val.get_bounding_box()[0]
            << val.get_bounding_box()[1]
            << val.get_bounding_box()[2]
            << val.get_bounding_box()[3]
            << val.get_bounding_box()[4]
            << val.get_bounding_box()[5]
            << val.get_centroid_position()[0]
            << val.get_centroid_position()[1]
            << val.get_centroid_position()[2]
            << val.get_centroid_variance()[0]
            << val.get_centroid_variance()[1]
            << val.get_centroid_variance()[2]
            << val.get_centroid_sq_width()[0]
            << val.get_centroid_sq_width()[1]
            << val.get_centroid_sq_width()[2]
            << val.get_intensity()
            << val.get_intensity_variance();

      profile_to_string(val.get_shoebox());
      profile_to_string(val.get_shoebox_mask());
      profile_to_string(val.get_shoebox_background());
      profile_to_string(val.get_transformed_shoebox());

      return *this;
    }
  };

  struct from_string : scitbx::af::boost_python::pickle_double_buffered::from_string
  {
    using scitbx::af::boost_python::pickle_double_buffered::from_string::operator>>;

    from_string(const char* str_ptr)
    : scitbx::af::boost_python::pickle_double_buffered::from_string(str_ptr)
    {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    template <typename ProfileType>
    void profile_from_string(ProfileType &p) {
      typename ProfileType::index_type shape;
      typename ProfileType::size_type n_dim;
      *this >> n_dim;
      shape.resize(n_dim);
      std::size_t capacity = 0;
      for (std::size_t i = 0; i < n_dim; ++i) {
        *this >> shape[i];
        capacity *= shape[i];
      }
      p.resize(capacity);
      p.resize(flex_grid<>(shape));
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this >> p[i];
      }
      std::cout << p.size() << std::endl;
    }

    from_string& operator>>(Reflection& val)
    {
      *this >> val.miller_index_[0]
            >> val.miller_index_[1]
            >> val.miller_index_[2]
            >> val.status_
            >> val.entering_
            >> val.rotation_angle_
            >> val.beam_vector_[0]
            >> val.beam_vector_[1]
            >> val.beam_vector_[2]
            >> val.image_coord_px_[0]
            >> val.image_coord_px_[1]
            >> val.image_coord_mm_[0]
            >> val.image_coord_mm_[1]
            >> val.frame_number_
            >> val.panel_number_
            >> val.bounding_box_[0]
            >> val.bounding_box_[1]
            >> val.bounding_box_[2]
            >> val.bounding_box_[3]
            >> val.bounding_box_[4]
            >> val.bounding_box_[5]
            >> val.centroid_position_[0]
            >> val.centroid_position_[1]
            >> val.centroid_position_[2]
            >> val.centroid_variance_[0]
            >> val.centroid_variance_[1]
            >> val.centroid_variance_[2]
            >> val.centroid_sq_width_[0]
            >> val.centroid_sq_width_[1]
            >> val.centroid_sq_width_[2]
            >> val.intensity_
            >> val.intensity_variance_;

      profile_from_string(val.shoebox_);
      profile_from_string(val.shoebox_mask_);
      profile_from_string(val.shoebox_background_);
      profile_from_string(val.transformed_shoebox_);
      return *this;
    }

    unsigned int version;
  };

}}}} // namespace dials::model::boost_python::reflection

#endif /* DIALS_MODEL_DATA_BOOST_PYTHON_PICKLE_H */
