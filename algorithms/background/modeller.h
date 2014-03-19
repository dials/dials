/*
 * modeller.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_MODELLER_H
#define DIALS_ALGORITHMS_BACKGROUND_MODELLER_H

#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <scitbx/matrix/inversion.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  using scitbx::matrix::inversion_in_place;

  /**
   * An abtract class for the background model
   */
  class Model {
  public:
    virtual
    double value(double z, double y, double x) const = 0;
  };

  /**
   * An abstract class for the background modeller
   */
  class Modeller {
  public:
    virtual
    boost::shared_ptr<Model> create(
        const af::const_ref< double, af::c_grid<3> > &data,
        const af::const_ref< bool, af::c_grid<3> > &mask) const = 0;
  };

  /**
   * A background model that is constant per image.
   */
  class Constant2dModel : public Model {
  public:

    Constant2dModel(af::shared<double> a_)
      : a(a_) {}

    virtual
    double value(double z, double y, double x) const {
      int index = (int)std::floor(z);
      DIALS_ASSERT(index >= 0 && index < a.size());
      return a[index];
    }

  private:
    af::shared<double> a;
  };

  /**
   * Create a background that is constant per image.
   */
  class Constant2dModeller : public Modeller {
  public:

    virtual
    boost::shared_ptr<Model> create(
        const af::const_ref< double, af::c_grid<3> > &data,
        const af::const_ref< bool, af::c_grid<3> > &mask) const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      af::shared<double> mean(data.accessor()[0], 0.0);
      for (std::size_t k = 0; k < data.accessor()[0]; ++k) {
        std::size_t count = 0;
        for (std::size_t j = 0; j < data.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < data.accessor()[2]; ++i) {
            if (mask(k, j, i)) {
              mean[k] += data(k, j, i);
              count++;
            }
          }
        }
        if (count > 0) {
          mean[k] /= count;
        }
      }
      return boost::make_shared<Constant2dModel>(mean);
    }
  };

  /**
   * A background model that is constant over the reflection.
   */
  class Constant3dModel : public Model {
  public:

    Constant3dModel(double a_)
      : a(a_) {}

    virtual
    double value(double z, double y, double x) const {
      return a;
    }

  private:
    double a;
  };

  /**
   * Create a background model that is constant over the reflection.
   */
  class Constant3dModeller : public Modeller {
  public:

    virtual
    boost::shared_ptr<Model> create(
        const af::const_ref< double, af::c_grid<3> > &data,
        const af::const_ref< bool, af::c_grid<3> > &mask) const {
      DIALS_ASSERT(data.size() == mask.size());
      double mean = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if (mask[i]) {
          mean += data[i];
          count++;
        }
      }
      if (count > 0) {
        mean /= count;
      }
      return boost::make_shared<Constant3dModel>(mean);
    }
  };

  /**
   * A background model that is a plane per image.
   */
  class Linear2dModel : public Model {
  public:

    Linear2dModel(
        af::shared<double> a_,
        af::shared<double> b_,
        af::shared<double> c_)
      : a(a_), b(b_), c(c_) {
      DIALS_ASSERT(a.size() == b.size());
      DIALS_ASSERT(a.size() == c.size());
    }

    virtual
    double value(double z, double y, double x) const {
      int index = (int)std::floor(z);
      DIALS_ASSERT(index >= 0 && index < c.size());
      return a[index] + b[index]*x + c[index]*y;
    }

  private:
    af::shared<double> a, b, c;
  };

  /**
   * Create a background model that is a plane per image.
   */
  class Linear2dModeller : public Modeller {
  public:

    virtual
    boost::shared_ptr<Model> create(
      const af::const_ref< double, af::c_grid<3> > &data,
      const af::const_ref< bool, af::c_grid<3> > &mask) const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      af::shared<double> a(mask.accessor()[0], 0);
      af::shared<double> b(mask.accessor()[0], 0);
      af::shared<double> c(mask.accessor()[0], 0);
      for (std::size_t k = 0; k < mask.accessor()[0]; ++k) {
        af::shared<double> A(3*3, 0);
        af::shared<double> B(3, 0);
        for (std::size_t j = 0; j < mask.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < mask.accessor()[2]; ++i) {
            if (mask(k,j,i)) {
              double x = (i + 0.5);
              double y = (j + 0.5);
              double p = data(k,j,i);
              A[0] += 1;
              A[1] += x;
              A[2] += y;
              A[4] += x * x;
              A[5] += x * y;
              A[8] += y * y;
              B[0] += p;
              B[1] += x * p;
              B[2] += y * p;
            }
          }
        }
        A[3] = A[1];
        A[6] = A[2];
        A[7] = A[5];
        inversion_in_place(A.begin(), 3, B.begin(), 1);
        a[k] = B[0];
        b[k] = B[1];
        c[k] = B[2];
      }
      return boost::make_shared<Linear2dModel>(a, b, c);
    }
  };

  /**
   * A background model that is a 3d hyper-plane.
   */
  class Linear3dModel : public Model {
  public:

    Linear3dModel(double a_, double b_, double c_, double d_)
      : a(a_), b(b_), c(c_), d(d_) {}

    virtual
    double value(double z, double y, double x) const {
      return a + b*x + c*y + d*z;
    }

  private:
    double a, b, c, d;
  };

  /**
   * Create a background model that is a 3d hyper-plane.
   */
  class Linear3dModeller : public Modeller {
  public:

    virtual
    boost::shared_ptr<Model> create(
      const af::const_ref< double, af::c_grid<3> > &data,
      const af::const_ref< bool, af::c_grid<3> > &mask) const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      af::shared<double> A(4*4, 0);
      af::shared<double> B(4, 0);
      for (std::size_t k = 0; k < mask.accessor()[0]; ++k) {
        for (std::size_t j = 0; j < mask.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < mask.accessor()[2]; ++i) {
            if (mask(k,j,i)) {
              double x = (i + 0.5);
              double y = (j + 0.5);
              double z = (k + 0.5);
              double p = data(k,j,i);
              A[0] += 1;
              A[1] += x;
              A[2] += y;
              A[3] += z;
              A[5] += x * x;
              A[6] += x * y;
              A[7] += x * z;
              A[10] += y * y;
              A[11] += y * z;
              A[15] += z * z;
              B[0] += p;
              B[1] += x * p;
              B[2] += y * p;
              B[3] += z * p;
            }
          }
        }
      }
      A[4] = A[1];
      A[8] = A[2];
      A[9] = A[6];
      A[12] = A[3];
      A[13] = A[7];
      A[14] = A[11];
      inversion_in_place(A.begin(), 4, B.begin(), 1);
      return boost::make_shared<Linear3dModel>(
          B[0], B[1], B[2], B[3]);
    }

  };

}}} // namespace dials::algorithms::background

#endif // DIALS_ALGORITHMS_BACKGROUND_MODELLER_H
