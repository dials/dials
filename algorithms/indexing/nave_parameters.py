#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

from __future__ import absolute_import, division
from scitbx.matrix import sqr,col
from dials.array_family import flex
import math

"""
Class to determine mosaicity and effective domain size for a crystal given a set of indexed reflections
"""

class nave_parameters(object):
  def __init__(self, params, experiments, reflections, refinery, graph_verbose=True):
    self.params = params
    self.experiments = experiments
    self.reflections = reflections
    self.refinery = refinery
    self.graph_verbose = graph_verbose

  def __call__(self):
    """Determine optimal mosaicity and domain size model (monochromatic)"""
    if self.refinery is None:
      RR = self.reflections
    else:
      RR = self.refinery.predict_for_reflection_table(self.reflections)
    excursion_rad = RR["delpsical.rad"]
    delta_psi_deg = excursion_rad * 180./math.pi
    print
    print flex.max(delta_psi_deg), flex.min(delta_psi_deg)
    mean_excursion = flex.mean(delta_psi_deg)
    print "The mean excursion is %7.3f degrees, r.m.s.d %7.3f"%(mean_excursion, math.sqrt(flex.mean(RR["delpsical2"])))

    crystal = self.experiments[0].crystal
    beam = self.experiments[0].beam
    miller_indices = self.reflections["miller_index"]

    # FIXME XXX revise this formula so as to use a different wavelength potentially for each reflection
    two_thetas = crystal.get_unit_cell().two_theta(miller_indices,beam.get_wavelength(),deg=True)
    dspacings = crystal.get_unit_cell().d(miller_indices)
    dspace_sq = dspacings * dspacings

    #  First -- try to get a reasonable envelope for the observed excursions.
        ## minimum of three regions; maximum of 50 measurements in each bin
    print "fitting parameters on %d spots"%len(excursion_rad)
    n_bins = min(max(3, len(excursion_rad)//25),50)
    bin_sz = len(excursion_rad)//n_bins
    print "nbins",n_bins,"bin_sz",bin_sz
    order = flex.sort_permutation(two_thetas)
    two_thetas_env = flex.double()
    dspacings_env = flex.double()
    excursion_rads_env = flex.double()
    for x in xrange(0,n_bins):
      subset = order[x*bin_sz:(x+1)*bin_sz]
      two_thetas_env.append(flex.mean(two_thetas.select(subset)))
      dspacings_env.append(flex.mean(dspacings.select(subset)))
      excursion_rads_env.append(flex.max(flex.abs(excursion_rad.select(subset))))

    #  Second -- parameter fit
        ## solve the normal equations
    sum_inv_u_sq = flex.sum(dspacings_env * dspacings_env)
    sum_inv_u    = flex.sum(dspacings_env)
    sum_te_u     = flex.sum(dspacings_env * excursion_rads_env)
    sum_te       = flex.sum(excursion_rads_env)
    Normal_Mat   = sqr((sum_inv_u_sq, sum_inv_u, sum_inv_u, len(dspacings_env)))
    Vector       = col((sum_te_u, sum_te))
    solution     = Normal_Mat.inverse() * Vector
    s_ang = 1./(2*solution[0])
    print "Best LSQ fit Scheerer domain size is %9.2f ang"%(
      s_ang)

    tan_phi_rad = dspacings / (2. * s_ang)
    tan_phi_deg = tan_phi_rad * 180./math.pi
    k_degrees = solution[1]* 180./math.pi
    print "The LSQ full mosaicity is %8.5f deg; half-mosaicity %9.5f"%(2*k_degrees, k_degrees)
    tan_outer_deg = tan_phi_deg + k_degrees

    from xfel.mono_simulation.max_like import minimizer
    # coerce the estimates to be positive for max-likelihood
    lower_limit_domain_size = math.pow(crystal.get_unit_cell().volume(),
     1./3.)*3 # params.refinement.domain_size_lower_limit

    d_estimate = max(s_ang, lower_limit_domain_size)
    M = minimizer(d_i = dspacings, psi_i = excursion_rad, eta_rad = abs(2. * solution[1]),
                  Deff = d_estimate)
    print "ML: mosaicity FW=%4.2f deg, Dsize=%5.0fA on %d spots"%(M.x[1]*180./math.pi, 2./M.x[0], len(two_thetas))
    tan_phi_rad_ML = dspacings / (2. / M.x[0])
    tan_phi_deg_ML = tan_phi_rad_ML * 180./math.pi
    tan_outer_deg_ML = tan_phi_deg_ML + 0.5*M.x[1]*180./math.pi

    self.nv_acceptance_flags = flex.abs(delta_psi_deg) < tan_outer_deg_ML

    if self.graph_verbose: #params.refinement.mosaic.enable_AD14F7B: # Excursion vs resolution fit
      AD1TF7B_MAX2T = 30.
      AD1TF7B_MAXDP = 1.
      from matplotlib import pyplot as plt
      plt.plot(two_thetas, delta_psi_deg, "bo")
      minplot = flex.min(two_thetas)
      plt.plot([0,minplot],[mean_excursion,mean_excursion],"k-")
      LR = flex.linear_regression(two_thetas, delta_psi_deg)
      model_y = LR.slope()*two_thetas + LR.y_intercept()
      plt.plot(two_thetas, model_y, "k-")

      plt.title("ML: mosaicity FW=%4.2f deg, Dsize=%5.0fA on %d spots"%(M.x[1]*180./math.pi, 2./M.x[0], len(two_thetas)))
      plt.plot(two_thetas, tan_phi_deg_ML, "r.")
      plt.plot(two_thetas, -tan_phi_deg_ML, "r.")
      plt.plot(two_thetas, tan_outer_deg_ML, "g.")
      plt.plot(two_thetas, -tan_outer_deg_ML, "g.")
      plt.xlim([0,AD1TF7B_MAX2T])
      plt.ylim([-AD1TF7B_MAXDP,AD1TF7B_MAXDP])
      plt.show()
      plt.close()

    from xfel.mono_simulation.util import green_curve_area
    self.green_curve_area = green_curve_area(two_thetas, tan_outer_deg_ML)
    print "The green curve area is ", self.green_curve_area

    crystal._ML_half_mosaicity_deg = M.x[1]*180./(2.*math.pi)
    crystal._ML_domain_size_ang = 2./M.x[0]
    self._ML_full_mosaicity_rad = M.x[1]
    self._ML_domain_size_ang = 2./M.x[0]

    #params.refinement.mosaic.model_expansion_factor
    """The expansion factor should be initially set to 1, then expanded so that the # reflections matched becomes
    as close as possible to # of observed reflections input, in the last integration call.  Determine this by
    inspecting the output log file interactively.  Do not exceed the bare minimum threshold needed.
    The intention is to find an optimal value, global for a given dataset."""
    model_expansion_factor = 1.4
    crystal._ML_half_mosaicity_deg *= model_expansion_factor
    crystal._ML_domain_size_ang /= model_expansion_factor

    return crystal

  def ewald_proximal_volume(self):
    """computes the volume of reciprocal space (actually, half the volume, in this implementation) in which
    reciprocal lattice centroids will fall under the green curve.  In other words, this is proportional to the
    number of predicted reflections."""

    R_L = 1./self.experiments[0].beam.get_wavelength() # radius of Ewald sphere

    # TT is the outermost two-theta angle to perform the volume integration (hi-resolution cutoff)
    TT = 2. * math.asin(self.experiments[0].beam.get_wavelength() /
                         (2. * self.params.indexing.stills.ewald_proximity_resolution_cutoff))

    part_vol = math.pi * (2./3.) * (1. - math.cos(TT))
    Ewald_sphere_volume = part_vol * math.pow(R_L, 3.) # base volume of Ewald sphere segment
    R_prime = R_L + 1./self._ML_domain_size_ang
    domain_size_volume = part_vol * math.pow(R_prime, 3.) # expanded volume accomodating spot size

    # compicated integral for mosaic spread volume, must be calculated numerically
    summation = 0.
    N_terms = 100
    for x in xrange(N_terms):
      phi = (x/N_terms) * TT
      # inner integral over radius r
      integral = math.pow(R_prime + (self._ML_full_mosaicity_rad * R_L * math.sin(phi)/2.), 3.) - \
                 math.pow(R_prime, 3.)
      summation += (integral * math.sin(phi)) * (TT/N_terms)
    mosaicity_volume = (2./3.) * math.pi * summation

    return (domain_size_volume - Ewald_sphere_volume) + mosaicity_volume
