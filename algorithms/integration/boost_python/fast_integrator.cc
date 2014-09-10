/*
 * fast_integrator.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
/* #include <boost/python.hpp> */
/* #include <boost/python/def.hpp> */
/* #include <boost/python/iterator.hpp> */
/* #include <dials/algorithms/integration/fast_integrator.h> */

/* namespace dials { namespace algorithms { namespace boost_python { */

/*   using namespace boost::python; */

/*   class WorkerIterator { */
/*   public: */

/*     typedef std::size_t base_iterator; */
/*     typedef ptrdiff_t difference_type; */
/*     typedef std::forward_iterator_tag iterator_category; */
/*     typedef FastIntegratorWorker value_type; */
/*     typedef const value_type *pointer; */
/*     typedef const value_type reference; */

/*     WorkerIterator(const FastIntegrator &self, base_iterator it) */
/*       : self_(self), */
/*         it_(it) {} */

/*     reference operator*() { */
/*       return self_.worker(it_); */
/*     } */

/*     WorkerIterator& operator++() { */
/*       ++it_; */
/*       return *this; */
/*     } */

/*     WorkerIterator operator++(int) { */
/*       WorkerIterator result(*this); */
/*       ++(*this); */
/*       return result; */
/*     } */

/*     bool operator==(const WorkerIterator& rhs) const { */
/*       return it_ == rhs.it_; */
/*     } */

/*     bool operator!=(const WorkerIterator& rhs) const { */
/*       return !(*this == rhs); */
/*     } */

/*   private: */
/*     const FastIntegrator &self_; */
/*     base_iterator it_; */
/*   }; */

/*   struct make_worker_iterator { */
/*     static */
/*     WorkerIterator begin(const FastIntegrator &self) { */
/*       return WorkerIterator(self, 0); */
/*     } */

/*     static */
/*     WorkerIterator end(const FastIntegrator &self) { */
/*       return WorkerIterator(self, self.size()); */
/*     } */

/*     static */
/*     object range() { */
/*       return boost::python::range( */
/*         &make_worker_iterator::begin, */
/*         &make_worker_iterator::end); */
/*     } */
/*   }; */

/*   struct FastIntegratorResultPickleSuite : */
/*       public boost::python::pickle_suite { */

/*     static */
/*     boost::python::tuple */
/*     getinitargs(const FastIntegratorResult &self) { */
/*       return boost::python::make_tuple( */
/*           self.index(), */
/*           self.intensity(), */
/*           self.variance() */
/*           ); */
/*     } */
/*   }; */

/*   void export_fast_integrator() */
/*   { */
/*     class_<FastIntegratorResult>("FastIntegratorResult", no_init) */
/*       .def(init<std::size_t, */
/*                 const af::const_ref<double>&, */
/*                 const af::const_ref<double>& >()) */
/*       .def_pickle(FastIntegratorResultPickleSuite()) */
/*       ; */

/*     class_<FastIntegratorWorker>("FastIntegratorWorker", no_init) */
/*       .def("first", &FastIntegratorWorker::first) */
/*       .def("last", &FastIntegratorWorker::last) */
/*       .def("current", &FastIntegratorWorker::current) */
/*       .def("next", &FastIntegratorWorker::next) */
/*       .def("finished", &FastIntegratorWorker::finished) */
/*       .def("result", &FastIntegratorWorker::result) */
/*       ; */

/*     class_<FastIntegrator>("FastIntegratorInternal", no_init) */
/*       .def(init<af::reflection_table, */
/*                 int, */
/*                 int, */
/*                 std::size_t>()) */
/*       .def("__len__", &FastIntegrator::size) */
/*       .def("finished", &FastIntegrator::finished) */
/*       .def("accumulate", &FastIntegrator::accumulate) */
/*       .def("worker", &FastIntegrator::worker) */
/*       .def("result", &FastIntegrator::result) */
/*       .def("workers", make_worker_iterator::range()) */
/*       ; */
/*   } */

/* }}} // namespace = dials::algorithms::boost_python */
