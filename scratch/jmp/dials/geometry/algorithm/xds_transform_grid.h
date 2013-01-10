
#ifndef DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H
#define DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H

#include <exception>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace geometry { namespace algorithm {

/** Class representing the XDS transform grid */
class xds_transform_grid {

public:

    /** Default constructor */
    xds_transform_grid() {}

    /** 
     * Initialise grid 
     */
    xds_transform_grid(std::size_t n_ref,
                       scitbx::vec3 <std::size_t> origin,
                       double sigma_divergence,
                       double sigma_mosaicity,
                       double n_sigma = 10.0)
        : _n_ref(n_ref),
          _origin(origin),
          _sigma_divergence(sigma_divergence),
          _sigma_mosaicity(sigma_mosaicity),
          _delta_divergence(sigma_divergence * n_sigma), 
          _delta_mosaicity(sigma_mosaicity  * n_sigma)
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
        _size[0] = 2 * _origin[0] + 1;
        _size[1] = 2 * _origin[1] + 1;
        _size[2] = 2 * _origin[2] + 1;
    
        // Calculate the step size of the grid
        _step_size[0] = _delta_divergence / _size[0];
        _step_size[1] = _delta_divergence / _size[1];
        _step_size[2] = _delta_mosaicity  / _size[2];
        
        // Allocate memory for the grid array
        _data = scitbx::af::flex_double(scitbx::af::flex_grid <> (
            n_ref, _size[0], _size[1], _size[2]));
    }

public:

    std::size_t get_n_reflections() {
        return _n_ref;
    }
    
    scitbx::vec3 <std::size_t> get_size() {
        return _size;
    }

    scitbx::vec3 <std::size_t> get_origin() {
        return _origin;
    }

    scitbx::vec3 <double> get_step_size() {
        return _step_size;
    }
    
    double get_sigma_divergence() {
        return _sigma_divergence;
    }
    
    double get_sigma_mosaicity() {
        return _sigma_mosaicity;
    }
    
    double get_delta_divergence() {
        return _delta_divergence;
    }
    
    double get_delta_mosaicity() {
        return _delta_mosaicity;
    }

    scitbx::af::flex_double get_data() {
        return _data;
    }

private:

    scitbx::af::flex_double _data;
    std::size_t _n_ref;
    scitbx::vec3 <std::size_t> _origin;
    scitbx::vec3 <std::size_t> _size;
    scitbx::vec3 <double> _step_size;
    double _sigma_divergence;
    double _sigma_mosaicity;
    double _delta_divergence;
    double _delta_mosaicity;
};

}}}

#endif // DIALS_GEOMETRY_XDS_TRANSFORM_GRID_H
