#ifndef DIALS_SCALING_SCALING_HELPER_H
#define DIALS_SCALING_SCALING_HELPER_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/math/zernike.h>
#include <dials/error.h>

namespace dials_scaling {

  /**
   * Elementwise squaring of a matrix
   */
  scitbx::sparse::matrix<double> elementwise_square(scitbx::sparse::matrix<double> m){

    scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

    // outer loop iterate over the columns
    for (std::size_t j=0; j < m.n_cols(); j++) {

      // inner loop iterate over the non-zero elements of the column
      for (scitbx::sparse::matrix<double>::row_iterator p=m.col(j).begin(); p != m.col(j).end(); ++p)
      {
        std::size_t i = p.index();
        result(i, j) = *p * *p;
      }
    }
    return result;
  }

  scitbx::sparse::matrix<double> row_multiply(scitbx::sparse::matrix<double> m,
                                              scitbx::af::const_ref<double> v){

    DIALS_ASSERT(m.n_rows() == v.size());

    // call compact to ensure that each elt of the matrix is only defined once
    m.compact();

    scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

    // outer loop iterate over the columns
    for (std::size_t j=0; j < m.n_cols(); j++) {

      // inner loop iterate over the non-zero elements of the column
      for (scitbx::sparse::matrix<double>::row_iterator p=m.col(j).begin(); p != m.col(j).end(); ++p)
      {
        std::size_t i = p.index();
        result(i, j) = *p * v[i];
      }
    }
    return result;
  }

  scitbx::af::shared<scitbx::vec2<double> > calc_theta_phi(
    scitbx::af::shared<scitbx::vec3<double> > xyz){
      int n_obs = xyz.size();
      scitbx::af::shared<scitbx::vec2<double> > theta_phi(n_obs);
      for (int i=0; i< n_obs; i++){
        theta_phi[i] = scitbx::vec2<double>(
          std::atan2(pow(pow(xyz[i][1], 2) + pow(xyz[i][0], 2), 0.5), xyz[i][2]),
          std::atan2(xyz[i][1], xyz[i][0]));
      }
      return theta_phi;
    }

  scitbx::af::shared<double> calc_sigmasq(
    scitbx::sparse::matrix<double> jacobian_transpose,
    scitbx::sparse::matrix<double> var_cov_matrix)
    {
      int n_cols = jacobian_transpose.n_cols();
      scitbx::af::shared<double> sigmasq(n_cols);
      for (int i=0; i < n_cols; i++){
        //scitbx::af::shared<double> result(n_refl);
        //for (int j=0; j < n_refl; j++){
        //  result[j] += jacobian_transpose.col(i) * var_cov_matrix.col(j);
        //  }

        for (scitbx::sparse::matrix<double>::row_iterator
          p = jacobian_transpose.col(i).begin();
          p != jacobian_transpose.col(i).end(); ++p){
            int k = p.index();
            sigmasq[i] += *p * (jacobian_transpose.col(i) * var_cov_matrix.col(k));// *p * result[k];
          }
      }
      return sigmasq;
    }


  scitbx::af::shared<scitbx::vec3<double> > rotate_vectors_about_axis(
    scitbx::af::shared<scitbx::vec3<double> > const rot_axis,
    scitbx::af::shared<scitbx::vec3<double> > const vectors,
    scitbx::af::shared<double> const angles)
    {
      // first normalise rotation axis
      double modulus = pow((pow(rot_axis[0][0],2) + pow(rot_axis[0][1], 2) +
        pow(rot_axis[0][2], 2)), 0.5);
      double ux = rot_axis[0][0]/modulus;
      double uy = rot_axis[0][1]/modulus;
      double uz = rot_axis[0][2]/modulus;
      int n_obs = angles.size();
      scitbx::af::shared<scitbx::vec3<double> > rotated_vectors(n_obs);

      for (int i=0; i < n_obs; i++){
        double cos_angle = std::cos(angles[i]);
        double sin_angle = std::sin(angles[i]);
        rotated_vectors[i] = scitbx::vec3<double>(
          (((cos_angle + ((pow(ux, 2)) * (1.0 - cos_angle))) * vectors[i][0])
            + (((ux * uy * (1.0 - cos_angle)) - (uz * sin_angle)) * vectors[i][1])
            + (((uz * ux * (1.0 - cos_angle)) + (uy * sin_angle)) * vectors[i][2])),
          ((((ux * uy * (1.0 - cos_angle)) + (uz * sin_angle)) * vectors[i][0])
            + ((cos_angle + ((pow(uy, 2)) * (1.0 - cos_angle))) * vectors[i][1])
            + (((uz * uy * (1.0 - cos_angle)) - (ux * sin_angle)) * vectors[i][2])),
          ((((ux * uz * (1.0 - cos_angle)) - (uy * sin_angle)) * vectors[i][0])
            + (((uy * uz * (1.0 - cos_angle)) + (ux * sin_angle)) * vectors[i][1])
            + ((cos_angle + ((pow(uz, 2)) * (1.0 - cos_angle))) * vectors[i][2]))
        );
        }
        return rotated_vectors;
      }

  /**
   * Spherical harmonic table
   */

  using scitbx::sparse::matrix;
  using scitbx::math::zernike::log_factorial_generator;
  using scitbx::math::zernike::nss_spherical_harmonics;
  using scitbx::sparse::vector;

  matrix<double> create_sph_harm_table(
      scitbx::af::shared<scitbx::vec2<double> > const s0_theta_phi,
      scitbx::af::shared<scitbx::vec2<double> > const s1_theta_phi,
      int lmax)
      {
        nss_spherical_harmonics<double> nsssphe(
          lmax, 50000, log_factorial_generator<double>((2 * lmax) + 1));
        int n_abs_param = (2 * lmax) + (pow(double(lmax), 2));
        int n_obs = s1_theta_phi.size();
        matrix<double> sph_harm_terms_(n_abs_param, n_obs);
        double sqrt2 = 1.414213562;
        int counter = 0;
        for (int l=1; l < lmax+1; l++) {
          for (int m=-1*l; m < l+1; m++) {
            if (m < 0) {
              double prefactor = sqrt2 * pow(-1.0, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(counter, i) = prefactor * (
                  nsssphe.spherical_harmonic_direct(l, -1*m, s0_theta_phi[i][0], s0_theta_phi[i][1]).imag()
                  + nsssphe.spherical_harmonic_direct(l, -1*m, s1_theta_phi[i][0], s1_theta_phi[i][1]).imag());
                }
              }
            else if (m == 0) {
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(counter, i) = (0.5 * (
                  nsssphe.spherical_harmonic_direct(l, 0, s0_theta_phi[i][0], s0_theta_phi[i][1]).real()
                  + nsssphe.spherical_harmonic_direct(l, 0, s1_theta_phi[i][0], s1_theta_phi[i][1]).real()));
              }
            }
            else {
              double prefactor = sqrt2 * pow(-1.0, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                double val = prefactor * (
                  nsssphe.spherical_harmonic_direct(l, m, s0_theta_phi[i][0], s0_theta_phi[i][1]).real()
                  + nsssphe.spherical_harmonic_direct(l, m, s1_theta_phi[i][0], s1_theta_phi[i][1]).real());
                sph_harm_terms_(counter, i) = val;
              }
            }
          counter += 1;
          }
        }
      return sph_harm_terms_;
      }

} // dials_scaling

#endif // DIALS_SCALING_SCALING_HELPER_H
