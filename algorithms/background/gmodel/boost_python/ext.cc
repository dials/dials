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
#include <dials/algorithms/background/gmodel/creator.h>
#include <dials/algorithms/background/gmodel/model.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  struct StaticBackgroundModelPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getstate(const StaticBackgroundModel &obj)
    {
      boost::python::list data;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        data.append(obj.data(i));
      }
      return boost::python::make_tuple(data);
    }

    static
    void setstate(StaticBackgroundModel &obj, boost::python::tuple state)
    {
      DIALS_ASSERT(boost::python::len(state) == 1);
      boost::python::list data = boost::python::extract<boost::python::list>(state[0]);
      for (std::size_t i = 0; i < boost::python::len(data); ++i) {
        af::const_ref< double, af::c_grid<2> > arr = boost::python::extract<
          af::const_ref< double, af::c_grid<2> > >(data[i]);
        obj.add(arr);
      }
    }
  };

  BOOST_PYTHON_MODULE(dials_algorithms_background_gmodel_ext)
  {
    class_<BackgroundModel, boost::noncopyable, boost::shared_ptr<BackgroundModel> >("BackgroundModel", no_init)
      .def("extract", pure_virtual(&BackgroundModel::extract))
      ;

    class_< StaticBackgroundModel, bases<BackgroundModel> >("StaticBackgroundModel")
      .def("add", &StaticBackgroundModel::add)
      .def("__len__", &StaticBackgroundModel::size)
      .def("data", &StaticBackgroundModel::data)
      .def_pickle(StaticBackgroundModelPickleSuite())
      ;

    class_<Creator> creator("Creator", no_init);
    creator
      .def(init<
          boost::shared_ptr<BackgroundModel>,
          bool,
          double,
          std::size_t>((
              arg("model"),
              arg("robust"),
              arg("tuning_constant"),
              arg("max_iter"))))
      .def("__call__", &Creator::shoebox)
      .def("__call__", &Creator::volume)
      ;

    class_<DispersionThreshold>("DispersionThreshold", no_init)
      .def(init< std::size_t,
                 std::size_t,
                 double,
                 double,
                 double,
                 int >())
      .def("__call__", &DispersionThreshold::threshold<int>)
      .def("__call__", &DispersionThreshold::threshold<double>)
      ;
  }

}}}} // namespace = dials::algorithms::background::boost_python
