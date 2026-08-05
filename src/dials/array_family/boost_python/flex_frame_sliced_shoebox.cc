/*
 * flex_frame_sliced_shoebox.cc
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/shoebox.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using dials::model::FrameSlicedShoebox;
  using dials::model::Shoebox;

  /**
   * Construct an array of frame sliced shoeboxes from an array of shoeboxes.
   *
   * @param phi The scan rotation angle in radians at the centre of each frame
   *            of the scan, as held by FrameOrientations.phi
   * @param first_frame The image number that phi[0] refers to, as held by
   *                    FrameOrientations.images[0]
   */
  template <typename FloatType>
  typename af::flex<FrameSlicedShoebox<FloatType> >::type* from_shoeboxes(
    const af::const_ref<Shoebox<FloatType> >& shoeboxes,
    const af::const_ref<double>& phi,
    int first_frame) {
    af::shared<FrameSlicedShoebox<FloatType> > result(shoeboxes.size());
    for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
      result[i] = FrameSlicedShoebox<FloatType>(shoeboxes[i], phi, first_frame);
    }
    return new typename af::flex<FrameSlicedShoebox<FloatType> >::type(
      result, af::flex_grid<>(result.size()));
  }

  /**
   * Get the number of frames of each frame sliced shoebox
   */
  template <typename FloatType>
  af::shared<std::size_t> num_frames(
    const af::const_ref<FrameSlicedShoebox<FloatType> >& self) {
    af::shared<std::size_t> result(self.size());
    for (std::size_t i = 0; i < self.size(); ++i) {
      result[i] = self[i].size();
    }
    return result;
  }

  /**
   * Construct an array of frame sliced shoeboxes from the flattened per-frame
   * arrays, as returned by get_frame_sliced_shoebox_data_arrays.
   *
   * @param num_frames The number of frames of each frame sliced shoebox, which
   *                   gives the lengths of the slices taken from the other
   *                   arrays
   */
  template <typename FloatType>
  typename af::flex<FrameSlicedShoebox<FloatType> >::type* from_data_arrays(
    const af::const_ref<std::size_t>& num_frames,
    const af::const_ref<int>& frames,
    const af::const_ref<double>& phi,
    const af::const_ref<int>& foreground_pixel_count,
    const af::const_ref<int>& valid_pixel_count,
    const af::const_ref<double>& foreground_sum_raw,
    const af::const_ref<double>& foreground_sum_minus_background,
    const af::const_ref<double>& summation_intensity,
    const af::const_ref<double>& summation_intensity_variance,
    const af::const_ref<bool>& summation_intensity_valid) {
    std::size_t total = 0;
    for (std::size_t i = 0; i < num_frames.size(); ++i) {
      total += num_frames[i];
    }
    DIALS_ASSERT(frames.size() == total);
    DIALS_ASSERT(phi.size() == total);
    DIALS_ASSERT(foreground_pixel_count.size() == total);
    DIALS_ASSERT(valid_pixel_count.size() == total);
    DIALS_ASSERT(foreground_sum_raw.size() == total);
    DIALS_ASSERT(foreground_sum_minus_background.size() == total);
    DIALS_ASSERT(summation_intensity.size() == total);
    DIALS_ASSERT(summation_intensity_variance.size() == total);
    DIALS_ASSERT(summation_intensity_valid.size() == total);

    af::shared<FrameSlicedShoebox<FloatType> > result(num_frames.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < num_frames.size(); ++i) {
      std::size_t nz = num_frames[i];
      if (nz == 0) {
        continue;
      }
      result[i] = FrameSlicedShoebox<FloatType>(
        af::shared<int>(&frames[k], &frames[k] + nz),
        af::shared<double>(&phi[k], &phi[k] + nz),
        af::shared<int>(&foreground_pixel_count[k], &foreground_pixel_count[k] + nz),
        af::shared<int>(&valid_pixel_count[k], &valid_pixel_count[k] + nz),
        af::shared<double>(&foreground_sum_raw[k], &foreground_sum_raw[k] + nz),
        af::shared<double>(&foreground_sum_minus_background[k],
                           &foreground_sum_minus_background[k] + nz),
        af::shared<double>(&summation_intensity[k], &summation_intensity[k] + nz),
        af::shared<double>(&summation_intensity_variance[k],
                           &summation_intensity_variance[k] + nz),
        af::shared<bool>(&summation_intensity_valid[k],
                         &summation_intensity_valid[k] + nz));
      k += nz;
    }
    return new typename af::flex<FrameSlicedShoebox<FloatType> >::type(
      result, af::flex_grid<>(result.size()));
  }

  /**
   * Append a per-frame array to the flat array collecting it. The source is
   * taken by value, as the getters of FrameSlicedShoebox return a copy.
   */
  template <typename ElementType>
  void extend_data_array(af::shared<ElementType>& dst,
                         const af::shared<ElementType> src) {
    dst.extend(src.begin(), src.end());
  }

  /**
   * Get the per-frame arrays of every frame sliced shoebox, each concatenated
   * into a single flat array. Use num_frames to split them up again.
   *
   * @returns A tuple of the frames, rotation angles, foreground pixel counts,
   *          valid pixel counts, raw foreground sums, background-subtracted
   *          foreground sums, summation intensities, their variances and their
   *          validity
   */
  template <typename FloatType>
  boost::python::tuple get_frame_sliced_shoebox_data_arrays(
    const af::const_ref<FrameSlicedShoebox<FloatType> >& self) {
    af::shared<int> frames;
    af::shared<double> phi;
    af::shared<int> foreground_pixel_count;
    af::shared<int> valid_pixel_count;
    af::shared<double> foreground_sum_raw;
    af::shared<double> foreground_sum_minus_background;
    af::shared<double> summation_intensity;
    af::shared<double> summation_intensity_variance;
    af::shared<bool> summation_intensity_valid;
    for (std::size_t i = 0; i < self.size(); ++i) {
      extend_data_array(frames, self[i].frames());
      extend_data_array(phi, self[i].phi());
      extend_data_array(foreground_pixel_count, self[i].foreground_pixel_count());
      extend_data_array(valid_pixel_count, self[i].valid_pixel_count());
      extend_data_array(foreground_sum_raw, self[i].foreground_sum_raw());
      extend_data_array(foreground_sum_minus_background,
                        self[i].foreground_sum_minus_background());
      extend_data_array(summation_intensity, self[i].summation_intensity());
      extend_data_array(summation_intensity_variance,
                        self[i].summation_intensity_variance());
      extend_data_array(summation_intensity_valid, self[i].summation_intensity_valid());
    }
    return boost::python::make_tuple(frames,
                                     phi,
                                     foreground_pixel_count,
                                     valid_pixel_count,
                                     foreground_sum_raw,
                                     foreground_sum_minus_background,
                                     summation_intensity,
                                     summation_intensity_variance,
                                     summation_intensity_valid);
  }

  /**
   * A class to convert the frame sliced shoebox class to a string for pickling
   */
  template <typename FloatType>
  struct frame_sliced_shoebox_to_string : pickle_double_buffered::to_string {
    using pickle_double_buffered::to_string::operator<<;

    typedef FrameSlicedShoebox<FloatType> frame_sliced_shoebox_type;

    /** Initialise with the version for checking */
    frame_sliced_shoebox_to_string() {
      unsigned int version = 3;
      *this << version;
    }

    /** Convert a single frame sliced shoebox instance to string */
    frame_sliced_shoebox_to_string& operator<<(const frame_sliced_shoebox_type& val) {
      // The arrays are all the same length, so only write it once
      *this << val.size();

      array_to_string(val.frames());
      array_to_string(val.phi());
      array_to_string(val.foreground_pixel_count());
      array_to_string(val.valid_pixel_count());
      array_to_string(val.foreground_sum_raw());
      array_to_string(val.foreground_sum_minus_background());
      array_to_string(val.summation_intensity());
      array_to_string(val.summation_intensity_variance());
      array_to_string(val.summation_intensity_valid());

      return *this;
    }

    /** Convert a per-frame array to string */
    template <typename ElementType>
    void array_to_string(const af::shared<ElementType>& a) {
      for (std::size_t i = 0; i < a.size(); ++i) {
        *this << a[i];
      }
    }
  };

  /**
   * A class to convert a string to a frame sliced shoebox for unpickling
   */
  template <typename FloatType>
  struct frame_sliced_shoebox_from_string : pickle_double_buffered::from_string {
    using pickle_double_buffered::from_string::operator>>;

    typedef FrameSlicedShoebox<FloatType> frame_sliced_shoebox_type;

    /** Initialise the class with the string. Get the version and check */
    frame_sliced_shoebox_from_string(const char* str_ptr)
        : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 3);
    }

    /** Get a single frame sliced shoebox instance from a string */
    frame_sliced_shoebox_from_string& operator>>(frame_sliced_shoebox_type& val) {
      std::size_t nz;
      *this >> nz;

      // NB the arrays are read into locals first, as the order in which function
      // arguments are evaluated is not specified
      af::shared<int> frames = array_from_string<int>(nz);
      af::shared<double> phi = array_from_string<double>(nz);
      af::shared<int> foreground_pixel_count = array_from_string<int>(nz);
      af::shared<int> valid_pixel_count = array_from_string<int>(nz);
      af::shared<double> foreground_sum_raw = array_from_string<double>(nz);
      af::shared<double> foreground_sum_minus_background =
        array_from_string<double>(nz);
      af::shared<double> summation_intensity = array_from_string<double>(nz);
      af::shared<double> summation_intensity_variance = array_from_string<double>(nz);
      af::shared<bool> summation_intensity_valid = array_from_string<bool>(nz);

      val = frame_sliced_shoebox_type(frames,
                                      phi,
                                      foreground_pixel_count,
                                      valid_pixel_count,
                                      foreground_sum_raw,
                                      foreground_sum_minus_background,
                                      summation_intensity,
                                      summation_intensity_variance,
                                      summation_intensity_valid);

      return *this;
    }

    /** Get a per-frame array from a string */
    template <typename ElementType>
    af::shared<ElementType> array_from_string(std::size_t n) {
      af::shared<ElementType> a(n);
      for (std::size_t i = 0; i < n; ++i) {
        *this >> a[i];
      }
      return a;
    }

    unsigned int version;
  };

  template <typename FloatType>
  typename scitbx::af::boost_python::
    flex_wrapper<FrameSlicedShoebox<FloatType>, return_internal_reference<> >::class_f_t
    flex_frame_sliced_shoebox_wrapper(const char* name) {
    typedef FrameSlicedShoebox<FloatType> frame_sliced_shoebox_type;

    return scitbx::af::boost_python::
      flex_wrapper<frame_sliced_shoebox_type, return_internal_reference<> >::plain(name)
        .def("__init__",
             make_constructor(from_shoeboxes<FloatType>,
                              default_call_policies(),
                              (boost::python::arg("shoeboxes"),
                               boost::python::arg("phi"),
                               boost::python::arg("first_frame"))))
        .def("__init__",
             make_constructor(from_data_arrays<FloatType>,
                              default_call_policies(),
                              (boost::python::arg("num_frames"),
                               boost::python::arg("frames"),
                               boost::python::arg("phi"),
                               boost::python::arg("foreground_pixel_count"),
                               boost::python::arg("valid_pixel_count"),
                               boost::python::arg("foreground_sum_raw"),
                               boost::python::arg("foreground_sum_minus_background"),
                               boost::python::arg("summation_intensity"),
                               boost::python::arg("summation_intensity_variance"),
                               boost::python::arg("summation_intensity_valid"))))
        .def("num_frames", &num_frames<FloatType>)
        .def("get_frame_sliced_shoebox_data_arrays",
             &get_frame_sliced_shoebox_data_arrays<FloatType>)
        .def_pickle(
          flex_pickle_double_buffered<frame_sliced_shoebox_type,
                                      frame_sliced_shoebox_to_string<FloatType>,
                                      frame_sliced_shoebox_from_string<FloatType> >());
  }

  void export_flex_frame_sliced_shoebox() {
    flex_frame_sliced_shoebox_wrapper<ProfileFloatType>("frame_sliced_shoebox");
  }

}}}  // namespace dials::af::boost_python
