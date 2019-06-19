/*
 * ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/spot_finding/helen/blobfinder.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  class BlobFinderWrapper {
  public:
    BlobFinderWrapper(int pixels_per_row,
                      int row_count,
                      int exp_spot_dimension,
                      double global_threshold,
                      double min_blob_score,
                      int num_passes)
        : alg_(pixels_per_row,
               row_count,
               exp_spot_dimension,
               global_threshold,
               min_blob_score,
               num_passes) {}

    template <typename T>
    void threshold(const af::const_ref<T, af::c_grid<2> > &image,
                   const af::const_ref<bool, af::c_grid<2> > &mask,
                   af::ref<bool, af::c_grid<2> > result) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(result.accessor()));

      // FIXME Make a copy (blob finder not marked as const)
      af::shared<T> image_tmp(image.size());
      af::shared<bool> mask_tmp(mask.size());
      std::copy(image.begin(), image.end(), image_tmp.begin());
      std::copy(mask.begin(), mask.end(), mask_tmp.begin());

      alg_.threshold(&image_tmp[0], &mask_tmp[0], &result[0]);
    }

    template <typename T>
    af::versa<bool, af::c_grid<2> > threshold_w_return(
      const af::const_ref<T, af::c_grid<2> > &image,
      const af::const_ref<bool, af::c_grid<2> > &mask) {
      af::versa<bool, af::c_grid<2> > result(image.accessor());
      threshold(image, mask, result.ref());
      return result;
    }

    template <typename T>
    af::versa<double, af::c_grid<2> > correlation_image(
      const af::const_ref<T, af::c_grid<2> > &image,
      const af::const_ref<bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      // FIXME Make a copy (blob finder not marked as const)
      af::shared<T> image_tmp(image.size());
      af::shared<bool> mask_tmp(mask.size());
      std::copy(image.begin(), image.end(), image_tmp.begin());
      std::copy(mask.begin(), mask.end(), mask_tmp.begin());
      af::versa<double, af::c_grid<2> > result(image.accessor());

      alg_.correlation_image(&image_tmp[0], &mask_tmp[0], &result[0]);

      return result;
    }

  private:
    BlobThresholdAlgorithm alg_;
  };

  BOOST_PYTHON_MODULE(dials_algorithms_spot_finding_helen_ext) {
    class_<BlobFinderWrapper>("BlobThresholdAlgorithm", no_init)
      .def(init<int, int, int, double, double, int>((arg("pixels_per_row"),
                                                     arg("row_count"),
                                                     arg("exp_spot_dimension") = 3,
                                                     arg("global_threshold") = 100,
                                                     arg("min_blob_score") = 0.7,
                                                     arg("num_passes") = 1)))
      .def("correlation", &BlobFinderWrapper::correlation_image<double>)
      .def("correlation", &BlobFinderWrapper::correlation_image<int>)
      .def("threshold", &BlobFinderWrapper::threshold<double>)
      .def("threshold", &BlobFinderWrapper::threshold<int>)
      .def("threshold", &BlobFinderWrapper::threshold_w_return<int>)
      .def("threshold", &BlobFinderWrapper::threshold_w_return<double>);
  }

}}}  // namespace dials::algorithms::boost_python
