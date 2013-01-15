
#ifndef DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H
#define DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H

#include <exception>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace geometry { namespace algorithm {

/** Class representing the XDS transform grid */
class XdsTransformGrid {

public:

    /** Default constructor */
    XdsTransformGrid() {}

    /** 
     * Initialise grid
     * @param n_ref The number of reflections
     * @param origin The origin of the grid (size = 2 * origin + 1)
     * @param sigma_divergence The standard deviation of the beam divergence
     * @param sigma_mosaicity The standard deviation of the mosaicity
     * @param n_sigma The number of standard deviations to cover
     */
    XdsTransformGrid(std::size_t n_ref,
                     scitbx::vec3 <std::size_t> origin,
                     double sigma_divergence,
                     double sigma_mosaicity,
                     double n_sigma = 10.0)
        : n_ref_(n_ref),
          origin_(origin),
          sigma_divergence_(sigma_divergence),
          sigma_mosaicity_(sigma_mosaicity),
          delta_divergence_(sigma_divergence * n_sigma), 
          delta_mosaicity_(sigma_mosaicity  * n_sigma)
    {
        // Check input to ensure grid is valid
        if (n_ref <= 0) {
            throw std::runtime_error("n_reflections must be > 0");
        }
        if (origin[0] <= 0 || origin[1] <= 0 || origin[2] <= 0) {
            throw std::runtime_error("origin must be > 0");
        }
        if (sigma_divergence <= 0.0 || sigma_mosaicity <= 0.0) {
            throw std::runtime_error("sigma_d and sigma_m must be > 0");
        }
        if (n_sigma <= 0.0) {
            throw std::runtime_error("n_sigma must be > 0");
        } 

        // Calculate the size of the grid
        size_[0] = 2 * origin_[0] + 1;
        size_[1] = 2 * origin_[1] + 1;
        size_[2] = 2 * origin_[2] + 1;
    
        // Calculate the step size of the grid
        step_size_[0] = delta_divergence_ / size_[0];
        step_size_[1] = delta_divergence_ / size_[1];
        step_size_[2] = delta_mosaicity_  / size_[2];
        
        // Allocate memory for the grid array
        data_ = scitbx::af::flex_double(scitbx::af::flex_grid <> (
            n_ref, size_[2], size_[1], size_[0]));
    }

public:

    /** Get the number of reflections */
    std::size_t get_n_reflections() {
        return n_ref_;
    }
    
    /** Get the size of the individual grids */
    scitbx::vec3 <std::size_t> get_size() {
        return size_;
    }

    /** Get the grid origin */
    scitbx::vec3 <std::size_t> get_origin() {
        return origin_;
    }

    /** Get the grid step sizes */
    scitbx::vec3 <double> get_step_size() {
        return step_size_;
    }
    
    /** Get the standard deviation of the beam divergence */
    double get_sigma_divergence() {
        return sigma_divergence_;
    }
    
    /** Get the standard deviation of the mosaicity */
    double get_sigma_mosaicity() {
        return sigma_mosaicity_;
    }
    
    /** Get the multiple of the standard deviation of the beam divergence */
    double get_delta_divergence() {
        return delta_divergence_;
    }
    
    /** Get the multiple of the standard deviation of the mosaicity */
    double get_delta_mosaicity() {
        return delta_mosaicity_;
    }

    /** Get the grid data */
    scitbx::af::flex_double get_data() {
        return data_;
    }

private:

    scitbx::af::flex_double data_;
    std::size_t n_ref_;
    scitbx::vec3 <std::size_t> origin_;
    scitbx::vec3 <std::size_t> size_;
    scitbx::vec3 <double> step_size_;
    double sigma_divergence_;
    double sigma_mosaicity_;
    double delta_divergence_;
    double delta_mosaicity_;
};

}}} // namespace = dials::geometry::algorithm

#endif // DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H
