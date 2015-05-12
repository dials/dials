/*
 * creator.cc
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
#include <dials/algorithms/background/simple/creator.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  af::shared<bool> call_1(
      const Creator<FloatType> &self,
      const af::const_ref< Shoebox<FloatType> > &sbox) {
    af::shared<double> mse(sbox.size());
    af::shared<double> dispersion(sbox.size());
    return self(sbox, mse.ref(), dispersion.ref());
  }

  template <typename FloatType>
  af::shared<bool> call_2(
      const Creator<FloatType> &self,
      const af::const_ref< Shoebox<FloatType> > &sbox,
      af::ref<double> mse,
      af::ref<double> dispersion) {
    return self(sbox, mse, dispersion);
  }

  template <typename FloatType>
  af::tiny<FloatType,2> call_3(
      const Creator<FloatType> &self,
      Shoebox<FloatType> shoebox) {
    return self(shoebox);
  }

  template <typename FloatType>
  af::tiny<FloatType,2> call_4(
      const Creator<FloatType> &self,
      const af::const_ref< FloatType, af::c_grid<3> > &data,
      af::ref< int, af::c_grid<3> > mask,
      af::ref< FloatType, af::c_grid<3> > background) {
    return self(data, mask, background);
  }


  template <typename FloatType>
  void creator_wrapper() {

    class_< Creator<FloatType> >("Creator", no_init)
      .def(init<
          boost::shared_ptr<Modeller>
        >((arg("modeller"))))
      .def(init<
          boost::shared_ptr<Modeller>,
          boost::shared_ptr<OutlierRejector>
        >((arg("modeller"),
           arg("rejector"))))
      .def("__call__", &call_1<FloatType>)
      .def("__call__", &call_2<FloatType>)
      .def("__call__", &call_3<FloatType>)
      .def("__call__", &call_4<FloatType>);
  }


  void export_creator()
  {
    creator_wrapper<float>();
  }

}}}} // namespace = dials::algorithms::background::boost_python
