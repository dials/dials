#ifndef DIALS_SCALING_SCALING_HELPER_H
#define DIALS_SCALING_SCALING_HELPER_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/math/zernike.h>
#include <dials/error.h>
#include <math.h>

typedef scitbx::sparse::matrix<double>::column_type col_type;

namespace dials_scaling {

using namespace boost::python;

/**
 * Elementwise squaring of a matrix
 */
scitbx::sparse::matrix<double> elementwise_square(scitbx::sparse::matrix<double> m) {
  scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

  // outer loop iterate over the columns
  for (std::size_t j = 0; j < m.n_cols(); j++) {
    // inner loop iterate over the non-zero elements of the column
    for (scitbx::sparse::matrix<double>::row_iterator p = m.col(j).begin();
         p != m.col(j).end();
         ++p) {
      std::size_t i = p.index();
      result(i, j) = *p * *p;
    }
  }
  return result;
}

scitbx::sparse::matrix<double> calculate_dIh_by_dpi(
  scitbx::af::shared<double> dIh,
  scitbx::af::shared<double> sumgsq,
  scitbx::sparse::matrix<double> h_index_mat,
  scitbx::sparse::matrix<double> derivatives) {
  // derivatives is a matrix where rows are params and cols are reflections
  int n_params = derivatives.n_rows();
  int n_groups = h_index_mat.n_cols();
  scitbx::sparse::matrix<double> dIh_by_dpi(n_groups, n_params);

  for (int i = 0; i < h_index_mat.n_cols(); ++i) {
    const col_type column = h_index_mat.col(i);
    // loop over reflection groups
    column.compact();
    for (col_type::const_iterator it = column.begin(); it != column.end(); ++it) {
      // so it.index gives index of a refl from group i
      int refl_idx = it.index();
      const col_type dgidx_by_dpi = derivatives.col(refl_idx);
      // deriv of one refl wrt all params
      for (col_type::const_iterator dgit = dgidx_by_dpi.begin();
           dgit != dgidx_by_dpi.end();
           ++dgit) {
        // dgit.index indicates which params have nonzero derivs
        dIh_by_dpi(i, dgit.index()) += (dIh[refl_idx] * *dgit / sumgsq[i]);
      }
    }
  }
  dIh_by_dpi.compact();
  return dIh_by_dpi;
}

scitbx::sparse::matrix<double> calculate_dIh_by_dpi_transpose(
  scitbx::af::shared<double> dIh,
  scitbx::af::shared<double> sumgsq,
  scitbx::sparse::matrix<double> h_index_mat,
  scitbx::sparse::matrix<double> derivatives) {
  // derivatives is a matrix where rows are params and cols are reflections
  int n_params = derivatives.n_rows();
  int n_groups = h_index_mat.n_cols();
  scitbx::sparse::matrix<double> dIh_by_dpi(n_params, n_groups);

  for (int i = 0; i < h_index_mat.n_cols(); ++i) {
    const col_type column = h_index_mat.col(i);
    // first loop over h_idx to get indices for a reflection group
    column.compact();
    scitbx::sparse::vector<double> deriv_of_group_by_params(n_params);
    for (col_type::const_iterator it = column.begin(); it != column.end(); ++it) {
      // it.index gives index of a refl from group i
      int refl_idx = it.index();
      const col_type dgidx_by_dpi = derivatives.col(refl_idx);
      // deriv of one refl wrt all params
      for (col_type::const_iterator dgit = dgidx_by_dpi.begin();
           dgit != dgidx_by_dpi.end();
           ++dgit) {
        // dgit.index indicates which params have nonzero derivs
        dIh_by_dpi(dgit.index(), i) += (dIh[refl_idx] * *dgit / sumgsq[i]);
      }
    }
  }
  dIh_by_dpi.compact();
  return dIh_by_dpi;
}

scitbx::sparse::matrix<double> calc_jacobian(scitbx::sparse::matrix<double> derivatives,
                                             scitbx::sparse::matrix<double> h_index_mat,
                                             scitbx::af::shared<double> Ih,
                                             scitbx::af::shared<double> g,
                                             scitbx::af::shared<double> dIh,
                                             scitbx::af::shared<double> sumgsq) {
  // derivatives is a matrix where rows are params and cols are reflections
  int n_params = derivatives.n_rows();
  int n_refl = derivatives.n_cols();

  scitbx::sparse::matrix<double> dIhbydpiT =
    calculate_dIh_by_dpi_transpose(dIh, sumgsq, h_index_mat, derivatives);
  scitbx::sparse::matrix<double> Jacobian(n_refl, n_params);

  for (int i = 0; i < h_index_mat.n_cols(); ++i) {
    const col_type column = h_index_mat.col(i);
    // first loop over h_idx to get indices for a reflection group
    column.compact();
    for (col_type::const_iterator it = column.begin(); it != column.end(); ++it) {
      // it.index gives index of a refl from group i
      int refl_idx = it.index();
      const col_type dgidx_by_dpi = derivatives.col(refl_idx);
      // deriv of one refl wrt all params
      dgidx_by_dpi.compact();
      // loop over nonzero elements of dgidx by dpi
      for (col_type::const_iterator dgit = dgidx_by_dpi.begin();
           dgit != dgidx_by_dpi.end();
           ++dgit) {
        // dgit.index indicates which params have nonzero derivs
        Jacobian(refl_idx, dgit.index()) -= *dgit * Ih[refl_idx];
      }
      // now loop over nonzero elements of dIhbydpi
      // get col corresponding to group
      const col_type dIh_col = dIhbydpiT.col(i);
      // now loop over nonzero params
      for (col_type::const_iterator dIit = dIh_col.begin(); dIit != dIh_col.end();
           ++dIit) {
        Jacobian(refl_idx, dIit.index()) -= g[refl_idx] * *dIit;
      }
    }
  }
  Jacobian.compact();
  return Jacobian;
}

scitbx::sparse::matrix<double> row_multiply(scitbx::sparse::matrix<double> m,
                                            scitbx::af::const_ref<double> v) {
  DIALS_ASSERT(m.n_rows() == v.size());

  // call compact to ensure that each elt of the matrix is only defined once
  m.compact();

  scitbx::sparse::matrix<double> result(m.n_rows(), m.n_cols());

  // outer loop iterate over the columns
  for (std::size_t j = 0; j < m.n_cols(); j++) {
    // inner loop iterate over the non-zero elements of the column
    for (scitbx::sparse::matrix<double>::row_iterator p = m.col(j).begin();
         p != m.col(j).end();
         ++p) {
      std::size_t i = p.index();
      result(i, j) = *p * v[i];
    }
  }
  return result;
}

scitbx::af::shared<scitbx::vec2<double> > calc_theta_phi(
  scitbx::af::shared<scitbx::vec3<double> > xyz) {
  // theta from -pi to pi. phi from 0 to pi
  int n_obs = xyz.size();
  scitbx::af::shared<scitbx::vec2<double> > theta_phi(n_obs);
  for (int i = 0; i < n_obs; i++) {
    theta_phi[i] = scitbx::vec2<double>(
      std::atan2(pow(pow(xyz[i][1], 2) + pow(xyz[i][0], 2), 0.5), xyz[i][2]),
      std::atan2(xyz[i][1], xyz[i][0]));
  }
  return theta_phi;
}

scitbx::af::shared<std::size_t> calc_lookup_index(
  scitbx::af::shared<scitbx::vec2<double> > theta_phi,
  double points_per_degree) {
  // theta index 0 to 2pi (shift by pi from input theta),
  scitbx::af::shared<std::size_t> lookup_index(theta_phi.size());
  for (int i = 0; i < theta_phi.size(); ++i) {
    lookup_index[i] =
      int(360.0 * points_per_degree
            * floor(theta_phi[i][0] * points_per_degree * 180.0 / scitbx::constants::pi)
          + floor((theta_phi[i][1] + scitbx::constants::pi) * 180.0 * points_per_degree
                  / scitbx::constants::pi));
  }
  return lookup_index;
}

boost::python::tuple determine_outlier_indices(
  scitbx::sparse::matrix<double> h_index_mat,
  scitbx::af::shared<double> z_scores,
  double zmax) {
  scitbx::af::shared<std::size_t> outlier_indices;
  scitbx::af::shared<std::size_t> other_potential_outlier_indices;
  for (int i = 0; i < h_index_mat.n_cols(); ++i) {
    const col_type column = h_index_mat.col(i);
    column.compact();
    double max_z = zmax;  // copy value//
    int n_elem = 0;
    int index_of_max = 0;
    for (col_type::const_iterator it = column.begin(); it != column.end(); ++it) {
      double val = z_scores[it.index()];
      if (val > max_z) {
        max_z = val;
        index_of_max = it.index();
      }
      ++n_elem;
    }
    if (n_elem > 2 && max_z > zmax) {
      // want to get indices of other potential outliers too
      outlier_indices.push_back(index_of_max);
      for (col_type::const_iterator it = column.begin(); it != column.end(); ++it) {
        if (it.index() != index_of_max) {
          other_potential_outlier_indices.push_back(it.index());
        }
      }
    }
  }
  return boost::python::make_tuple(outlier_indices, other_potential_outlier_indices);
}

scitbx::af::shared<double> calc_sigmasq(
  scitbx::sparse::matrix<double> jacobian_transpose,
  scitbx::sparse::matrix<double> var_cov_matrix) {
  int n_cols = jacobian_transpose.n_cols();
  scitbx::af::shared<double> sigmasq(n_cols);
  for (int i = 0; i < n_cols; i++) {
    // scitbx::af::shared<double> result(n_refl);
    // for (int j=0; j < n_refl; j++){
    //  result[j] += jacobian_transpose.col(i) * var_cov_matrix.col(j);
    //  }

    for (scitbx::sparse::matrix<double>::row_iterator p =
           jacobian_transpose.col(i).begin();
         p != jacobian_transpose.col(i).end();
         ++p) {
      int k = p.index();
      sigmasq[i] +=
        *p * (jacobian_transpose.col(i) * var_cov_matrix.col(k));  // *p * result[k];
    }
  }
  return sigmasq;
}

scitbx::af::shared<scitbx::vec3<double> > rotate_vectors_about_axis(
  scitbx::af::shared<scitbx::vec3<double> > const rot_axis,
  scitbx::af::shared<scitbx::vec3<double> > const vectors,
  scitbx::af::shared<double> const angles) {
  // first normalise rotation axis
  double modulus = pow(
    (pow(rot_axis[0][0], 2) + pow(rot_axis[0][1], 2) + pow(rot_axis[0][2], 2)), 0.5);
  double ux = rot_axis[0][0] / modulus;
  double uy = rot_axis[0][1] / modulus;
  double uz = rot_axis[0][2] / modulus;
  int n_obs = angles.size();
  scitbx::af::shared<scitbx::vec3<double> > rotated_vectors(n_obs);

  for (int i = 0; i < n_obs; i++) {
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
       + ((cos_angle + ((pow(uz, 2)) * (1.0 - cos_angle))) * vectors[i][2])));
  }
  return rotated_vectors;
}

/**
 * Spherical harmonic table
 */

using scitbx::math::zernike::log_factorial_generator;
using scitbx::math::zernike::nss_spherical_harmonics;
using scitbx::sparse::matrix;
using scitbx::sparse::vector;

boost::python::tuple calculate_harmonic_tables_from_selections(
  scitbx::af::shared<std::size_t> s0_selection,
  scitbx::af::shared<std::size_t> s1_selection,
  boost::python::list coefficients_list) {
  int n_refl = s0_selection.size();
  int n_param = boost::python::len(coefficients_list);
  boost::python::list output_coefficients_list;
  matrix<double> coefficients_matrix(n_refl, n_param);
  // loop though each param first, then loop over selection
  for (int i = 0; i < n_param; ++i) {
    scitbx::af::shared<double> coefs =
      boost::python::extract<scitbx::af::shared<double> >(coefficients_list[i]);
    scitbx::af::shared<double> coef_for_output(n_refl);
    // loop though each reflection
    for (int j = 0; j < n_refl; ++j) {
      double val0 = coefs[s0_selection[j]];
      double val1 = coefs[s1_selection[j]];
      double value = (val0 + val1) / 2.0;
      coefficients_matrix(j, i) = value;
      coef_for_output[j] = value;
    }
    output_coefficients_list.append(coef_for_output);
  }
  return boost::python::make_tuple(output_coefficients_list, coefficients_matrix);
}

matrix<double> create_sph_harm_table(
  scitbx::af::shared<scitbx::vec2<double> > const s0_theta_phi,
  scitbx::af::shared<scitbx::vec2<double> > const s1_theta_phi,
  int lmax) {
  nss_spherical_harmonics<double> nsssphe(
    lmax, 50000, log_factorial_generator<double>((2 * lmax) + 1));
  int n_abs_param = (2 * lmax) + (pow(double(lmax), 2));
  int n_obs = s1_theta_phi.size();
  matrix<double> sph_harm_terms_(n_abs_param, n_obs);
  double sqrt2 = 1.414213562;
  int counter = 0;
  for (int l = 1; l < lmax + 1; l++) {
    for (int m = -1 * l; m < l + 1; m++) {
      if (m < 0) {
        double prefactor = sqrt2 * pow(-1.0, m) / 2.0;
        for (int i = 0; i < n_obs; i++) {
          sph_harm_terms_(counter, i) =
            prefactor
            * (nsssphe
                 .spherical_harmonic_direct(
                   l, -1 * m, s0_theta_phi[i][0], s0_theta_phi[i][1])
                 .imag()
               + nsssphe
                   .spherical_harmonic_direct(
                     l, -1 * m, s1_theta_phi[i][0], s1_theta_phi[i][1])
                   .imag());
        }
      } else if (m == 0) {
        for (int i = 0; i < n_obs; i++) {
          sph_harm_terms_(counter, i) =
            (0.5
             * (nsssphe
                  .spherical_harmonic_direct(
                    l, 0, s0_theta_phi[i][0], s0_theta_phi[i][1])
                  .real()
                + nsssphe
                    .spherical_harmonic_direct(
                      l, 0, s1_theta_phi[i][0], s1_theta_phi[i][1])
                    .real()));
        }
      } else {
        double prefactor = sqrt2 * pow(-1.0, m) / 2.0;
        for (int i = 0; i < n_obs; i++) {
          double val = prefactor
                       * (nsssphe
                            .spherical_harmonic_direct(
                              l, m, s0_theta_phi[i][0], s0_theta_phi[i][1])
                            .real()
                          + nsssphe
                              .spherical_harmonic_direct(
                                l, m, s1_theta_phi[i][0], s1_theta_phi[i][1])
                              .real());
          sph_harm_terms_(counter, i) = val;
        }
      }
      counter += 1;
    }
  }
  return sph_harm_terms_;
}

boost::python::list create_sph_harm_lookup_table(int lmax, int points_per_degree) {
  nss_spherical_harmonics<double> nsssphe(
    lmax, 50000, log_factorial_generator<double>((2 * lmax) + 1));
  boost::python::list coefficients_list;
  double sqrt2 = 1.414213562;
  int n_items = 360 * 180 * points_per_degree * points_per_degree;
  double delta = scitbx::constants::pi / (360 * points_per_degree);
  for (int l = 1; l < lmax + 1; l++) {
    for (int m = -1 * l; m < l + 1; m++) {
      scitbx::af::shared<double> coefficients(n_items);
      if (m < 0) {
        double prefactor = sqrt2 * pow(-1.0, m);
        for (int i = 0; i < n_items; i++) {
          double theta = (floor(i / (360.0 * points_per_degree)) * scitbx::constants::pi
                          / (180.0 * points_per_degree))
                         + delta;
          double phi = (((i % (360 * points_per_degree)) * scitbx::constants::pi
                         / (180.0 * points_per_degree))
                        - scitbx::constants::pi)
                       + delta;
          coefficients[i] =
            prefactor
            * (nsssphe.spherical_harmonic_direct(l, -1 * m, theta, phi).imag());
        }
      } else if (m == 0) {
        for (int i = 0; i < n_items; i++) {
          double theta = (floor(i / (360.0 * points_per_degree)) * scitbx::constants::pi
                          / (180.0 * points_per_degree))
                         + delta;
          double phi = (((i % (360 * points_per_degree)) * scitbx::constants::pi
                         / (180.0 * points_per_degree))
                        - scitbx::constants::pi)
                       + delta;
          coefficients[i] = nsssphe.spherical_harmonic_direct(l, 0, theta, phi).real();
        }
      } else {
        double prefactor = sqrt2 * pow(-1.0, m);
        for (int i = 0; i < n_items; i++) {
          double theta = (floor(i / (360.0 * points_per_degree)) * scitbx::constants::pi
                          / (180.0 * points_per_degree))
                         + delta;
          double phi = (((i % (360 * points_per_degree)) * scitbx::constants::pi
                         / (180.0 * points_per_degree))
                        - scitbx::constants::pi)
                       + delta;
          coefficients[i] =
            prefactor * (nsssphe.spherical_harmonic_direct(l, m, theta, phi).real());
        }
      }
      coefficients_list.append(coefficients);
    }
  }
  return coefficients_list;
}

}  // namespace dials_scaling

#endif  // DIALS_SCALING_SCALING_HELPER_H
