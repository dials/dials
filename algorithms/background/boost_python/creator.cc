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
#include <dials/algorithms/background/creator.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  af::shared<bool> call_for_shoeboxes(
      const Creator &self,
      const af::const_ref< Shoebox<> > &sbox) {
    af::shared<double> mse(sbox.size());
    return self(sbox, mse.ref());
  }

  void export_creator()
  {
    typedef af::shared<bool>(Creator::*call_for_shoeboxes_w_err)(
        const af::const_ref< Shoebox<> >&,
        af::ref<double>)const;

    typedef double(Creator::*call_for_shoebox)(Shoebox<>)const;

    typedef double(Creator::*call_for_data)(
        const af::const_ref< double, af::c_grid<3> >&,
        af::ref< int, af::c_grid<3> >,
        af::ref< double, af::c_grid<3> >)const;

    class_<Creator>("Creator", no_init)
      .def(init<
          boost::shared_ptr<Modeller>
        >((arg("modeller"))))
      .def(init<
          boost::shared_ptr<Modeller>,
          boost::shared_ptr<OutlierRejector>
        >((arg("modeller"),
           arg("rejector"))))
      .def("__call__", &call_for_shoeboxes)
      .def("__call__", (call_for_shoeboxes_w_err)&Creator::operator())
      .def("__call__", (call_for_shoebox)&Creator::operator())
      .def("__call__", (call_for_data)&Creator::operator());
  }

}}}} // namespace = dials::algorithms::background::boost_python
