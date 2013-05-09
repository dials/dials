/*
 * poisson_discriminator.cc
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
#include <dials/algorithms/background/poisson_discriminator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_poisson_discriminator()
  {
    class_<PoissonDiscriminator, bases<DiscriminatorStrategy> >(
        "PoissonDiscriminator")
      .def(init<>());
  }

}}} // namespace = dials::algorithms::boost_python
