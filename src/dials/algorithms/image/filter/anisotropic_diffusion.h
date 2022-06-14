/*
 * anisotropic_diffusion.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_ANISOTROPIC_DIFFUSION_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_ANISOTROPIC_DIFFUSION_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * Do anisotropic filtering on an image
   * @param data The image
   * @param niter The number of iterations
   * @param kappa The diffusion parameter, small values stop diffusion across edges
   * @param gamma The step for each iteration (0 < gamma < 1.0)
   * @return The filtered image
   */
  inline af::versa<double, af::c_grid<2> > anisotropic_diffusion(
    const af::const_ref<double, af::c_grid<2> > &data,
    std::size_t niter,
    double kappa,
    double gamma) {
    // Check input
    DIALS_ASSERT(niter > 0);
    DIALS_ASSERT(kappa > 0);
    DIALS_ASSERT(gamma > 0);

    // Initialise stuff
    std::size_t height = data.accessor()[0];
    std::size_t width = data.accessor()[1];
    af::versa<double, af::c_grid<2> > AA(data.accessor());
    af::versa<double, af::c_grid<2> > BB(data.accessor(), 0);
    af::ref<double, af::c_grid<2> > A = AA.ref();
    af::ref<double, af::c_grid<2> > B = BB.ref();
    std::copy(data.begin(), data.end(), A.begin());

    // Compute inv of kappa
    double kappa_inv2 = 1.0 / (kappa * kappa);

    // Iterate
    for (std::size_t iter = 0; iter < niter; ++iter) {
      for (std::size_t j = 1; j < height - 1; ++j) {
        for (std::size_t i = 1; i < width - 1; ++i) {
          // Gradients
          double DN = A(j, i) - A(j - 1, i);
          double DE = A(j, i) - A(j, i - 1);
          double DS = A(j + 1, i) - A(j, i);
          double DW = A(j, i + 1) - A(j, i);

          // Diffusion stuff
          // double CN = std::exp(-(DN*DN*kappa_inv2));
          // double CE = std::exp(-(DE*DE*kappa_inv2));
          // double CS = std::exp(-(DS*DS*kappa_inv2));
          // double CW = std::exp(-(DW*DW*kappa_inv2));
          double CN = 1.0 / (1.0 + (DN * DN * kappa_inv2));
          double CE = 1.0 / (1.0 + (DE * DE * kappa_inv2));
          double CS = 1.0 / (1.0 + (DS * DS * kappa_inv2));
          double CW = 1.0 / (1.0 + (DW * DW * kappa_inv2));

          // Components
          double N = CN * DN;
          double E = CE * DE;
          double S = CS * DS;
          double W = CW * DW;

          // Update image
          B(j, i) = gamma * (S - N + W - E);
        }
      }
      for (std::size_t i = 0; i < B.size(); ++i) {
        A[i] += B[i];
      }
    }

    // Return filtered image
    return AA;
  }

  /**
   * Do anisotropic filtering on an image
   * @param data The image
   * @param mask The mask
   * @param niter The number of iterations
   * @param kappa The diffusion parameter, small values stop diffusion across edges
   * @param gamma The step for each iteration (0 < gamma < 1.0)
   * @return The filtered image
   */
  inline af::versa<double, af::c_grid<2> > masked_anisotropic_diffusion(
    const af::const_ref<double, af::c_grid<2> > &data,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    std::size_t niter,
    double kappa,
    double gamma) {
    // Check input
    DIALS_ASSERT(niter > 0);
    DIALS_ASSERT(kappa > 0);
    DIALS_ASSERT(gamma > 0);

    // Initialise stuff
    std::size_t height = data.accessor()[0];
    std::size_t width = data.accessor()[1];
    af::versa<double, af::c_grid<2> > AA(data.accessor());
    af::versa<double, af::c_grid<2> > BB(data.accessor(), 0);
    af::ref<double, af::c_grid<2> > A = AA.ref();
    af::ref<double, af::c_grid<2> > B = BB.ref();
    std::copy(data.begin(), data.end(), A.begin());

    // Compute inv of kappa
    double kappa_inv2 = 1.0 / (kappa * kappa);

    // Iterate
    for (std::size_t iter = 0; iter < niter; ++iter) {
      for (std::size_t j = 1; j < height - 1; ++j) {
        for (std::size_t i = 1; i < width - 1; ++i) {
          if (mask(j, i)) {
            // Points
            double AP = A(j, i);
            double AN = mask(j - 1, i) ? A(j - 1, i) : AP;
            double AS = mask(j + 1, i) ? A(j + 1, i) : AP;
            double AE = mask(j, i - 1) ? A(j, i - 1) : AP;
            double AW = mask(j, i + 1) ? A(j, i + 1) : AP;

            // Gradients
            double DN = AP - AN;
            double DE = AP - AE;
            double DS = AS - AP;
            double DW = AW - AP;

            // Diffusion stuff
            // double CN = std::exp(-(DN*DN*kappa_inv2));
            // double CE = std::exp(-(DE*DE*kappa_inv2));
            // double CS = std::exp(-(DS*DS*kappa_inv2));
            // double CW = std::exp(-(DW*DW*kappa_inv2));
            double CN = 1.0 / (1.0 + (DN * DN * kappa_inv2));
            double CE = 1.0 / (1.0 + (DE * DE * kappa_inv2));
            double CS = 1.0 / (1.0 + (DS * DS * kappa_inv2));
            double CW = 1.0 / (1.0 + (DW * DW * kappa_inv2));

            // Components
            double N = CN * DN;
            double E = CE * DE;
            double S = CS * DS;
            double W = CW * DW;

            // Update image
            B(j, i) = gamma * (S - N + W - E);
          } else {
            B(j, i) = 0;
          }
        }
      }
      for (std::size_t i = 0; i < B.size(); ++i) {
        A[i] += B[i];
      }
    }

    // Return filtered image
    return AA;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_FILTER_ANISOTROPIC_DIFFUSION_H
