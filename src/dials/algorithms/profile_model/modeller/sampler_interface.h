/*
 * sampler_interface.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SAMPLER_INTERFACE_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SAMPLER_INTERFACE_H

#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms {

  using scitbx::af::double3;
  using scitbx::af::int3;

  /**
   * Class for sampler interface
   */
  class SamplerIface {
  public:
    virtual ~SamplerIface() {}

    /**
     * @returns The total number of grid points
     */
    virtual std::size_t size() const = 0;

    /**
     * Find the nearest reference profile to the given point.
     * @param panel The panel
     * @param xyz The coordinate
     * @returns The index of the reference profile
     */
    virtual std::size_t nearest(std::size_t panel, double3 xyz) const = 0;

    /**
     * Find the nearest n reference profiles to the given point.
     * This is mostly used during learning to get the neighbouring profiles.
     * @param panel The panel
     * @param xyz The coordinate
     * @returns A list of reference profile indices
     */
    virtual af::shared<std::size_t> nearest_n(std::size_t panel, double3 xyz) const = 0;

    /**
     * Get the weight for the given profile at the given coordinate.
     * @param index The profile index
     * @param panel The panel
     * @param xyz The coordinate
     * @returns The weight (between 1.0 and 0.0)
     */
    virtual double weight(std::size_t index, std::size_t panel, double3 xyz) const = 0;

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    virtual double3 coord(std::size_t index) const = 0;

    /**
     * Return the neighbouring grid points.
     */
    virtual af::shared<std::size_t> neighbours(std::size_t index) const = 0;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SAMPLER_INTERFACE_H */
