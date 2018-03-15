#ifndef DIALS_SCRATCH_SCALING_SCALING_HELPER_H
#define DIALS_SCRATCH_SCALING_SCALING_HELPER_H

#include <scitbx/sparse/matrix.h>
#include <scitbx/math/zernike.h>

namespace dials_scratch { namespace scaling {
  
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

  /**
   * Spherical harmonic table
   */

  using scitbx::sparse::matrix;
  using scitbx::math::zernike::log_factorial_generator;
  using scitbx::math::zernike::nss_spherical_harmonics;
  using scitbx::sparse::vector;

  matrix<double> create_sph_harm_table(scitbx::af::shared<double> const s0_theta,
      scitbx::af::shared<double> const s0_phi,
      scitbx::af::shared<double> const s1_theta,
      scitbx::af::shared<double> const s1_phi,
      int lmax)
      {
        nss_spherical_harmonics<double> nsssphe(
          lmax, 50000, log_factorial_generator<double>((2 * lmax) + 1));
        int n_abs_param = (2 * lmax) + (pow(lmax, 2));
        int n_obs = s1_theta.size();
        matrix<double> sph_harm_terms_(n_obs, n_abs_param);
        double sqrt2 = 1.414213562;
        int counter = 0;
        for (int l=1; l < lmax+1; l++) {
          for (int m=-1*l; m < l+1; m++) {
            if (m < 0) {
              double prefactor = sqrt2 * pow(-1, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(i, counter) = prefactor * (
                  nsssphe.spherical_harmonic_direct(l, -1*m, s0_theta[i], s0_phi[i]).imag()
                  + nsssphe.spherical_harmonic_direct(l, -1*m, s1_theta[i], s1_phi[i]).imag());
                }
              }
            else if (m == 0) {
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(i, counter) = (0.5 * (
                  nsssphe.spherical_harmonic_direct(l, 0, s0_theta[i], s0_phi[i]).real()
                  + nsssphe.spherical_harmonic_direct(l, 0, s1_theta[i], s1_phi[i]).real()));
              }
            }
            else {
              double prefactor = sqrt2 * pow(-1, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                double val = prefactor * (
                  nsssphe.spherical_harmonic_direct(l, m, s0_theta[i], s0_phi[i]).real()
                  + nsssphe.spherical_harmonic_direct(l, m, s1_theta[i], s1_phi[i]).real());
                sph_harm_terms_(i, counter) = val;
              }
            }
          counter += 1;
          }
        }
      return sph_harm_terms_;
      }

  /*class sph_harm_table {
  public:
    
    sph_harm_table(scitbx::af::shared<double> s0_phi_,
      scitbx::af::shared<double> s0_theta_,
      scitbx::af::shared<double> s1_phi_,
      scitbx::af::shared<double> s1_theta_,
      int lmax_)
      : lmax(lmax_), s0_phi(s0_phi_), s0_theta(s0_theta_), s1_phi(s1_phi_),
        s1_theta(s1_theta_),
        nsssphe(lmax_, 50000, log_factorial_generator<double>(2 * lmax_ + 1)) {
        
        //std::cout<< "here \n";
        int n_abs_param = (2 * lmax) + (pow(lmax, 2));
        int n_obs = s1_theta.size();
        //std::cout << n_abs_param << "\n";
        //std::cout << n_obs << "\n";
        matrix<double> sph_harm_terms_(n_obs, n_abs_param);
        std::cout << sph_harm_terms_.n_cols() << std::endl;
        std::cout << sph_harm_terms_.n_rows() << std::endl;
        double sqrt2 = 1.414213562;
        //std::cout<<"here \n";
        // nss_spherical_harmonics<double> nsssphe(lmax, 50000, lgf);
        nsssphe = nss_spherical_harmonics(<int> lmax, <int> 50000,
          <log_factorial_generator<double>> lgf);

        int counter = 0;

        for (int l=1; l < lmax+1; l++) {
          for (int m= -1*l ; m < l+1; m++) {
            std::cout << l << std::endl;
            std::cout << m << std::endl;
            if (m < 0) {
              double prefactor = sqrt2 * pow(-1, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(i, counter) = prefactor * (
                  nsssphe.spherical_harmonic_pc(l, -1*m, s0_theta[i], s0_phi[i]).imag()
                  + nsssphe.spherical_harmonic_pc(l, -1*m, s1_theta[i], s1_phi[i]).imag());
                }
              }
            else if (m == 0) {
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(i, counter) = (0.5 * (
                  nsssphe.spherical_harmonic_pc(l, 0, s0_theta[i], s0_phi[i]).real()
                  + nsssphe.spherical_harmonic_pc(l, 0, s1_theta[i], s1_phi[i]).real()));
              }
            }
            else {
              double prefactor = sqrt2 * pow(-1, m) / 2.0;
              for (int i=0; i < n_obs; i++){
                sph_harm_terms_(i, counter) = prefactor * (
                  nsssphe.spherical_harmonic_pc(l, m, s0_theta[i], s0_phi[i]).real()
                  + nsssphe.spherical_harmonic_pc(l, m, s1_theta[i], s1_phi[i]).real());
              }
            counter += 1;
            }
          }
        }
      }

    scitbx::sparse::matrix<double> sph_harm_terms(){
      return sph_harm_terms_;
    }

  private:

    int lmax;
    scitbx::af::shared<double> s0_phi;
    scitbx::af::shared<double> s0_theta;
    scitbx::af::shared<double> s1_phi;
    scitbx::af::shared<double> s1_theta;
    scitbx::sparse::matrix<double> sph_harm_terms_;
    //scitbx::math::zernike::log_factorial_generator<double> lgf;
    scitbx::math::zernike::nss_spherical_harmonics<double> nsssphe;
  };*/
}} // dials_scratch::scaling

#endif // DIALS_SCRATCH_SCALING_SCALING_HELPER_H
