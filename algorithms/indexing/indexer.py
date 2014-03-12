#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.indexer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter, Nick Sauter, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math
from scitbx import matrix

from rstbx.array_family import flex
import copy


# we need these things

def discover_better_experimental_model(spot_positions, detector, beam,
                                       goniometer, scan, params):
  '''Given an attempt at indexing derive a more likely model for the
  experimental geometry.'''

  # Spot_positions: Centroid positions for spotfinder spots, in pixels
  # Return value: Corrected for parallax, converted to mm

  spots_mm = Indexer._map_spots_pixel_to_mm_rad(
    spots=spot_positions, detector=detector, scan=scan)

  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell = None,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  DPS.S0_vector = matrix.col(beam.get_s0())
  DPS.inv_wave = 1./beam.get_wavelength()
  if goniometer is None:
    DPS.axis = matrix.col((1,0,0))
  else:
    DPS.axis = matrix.col(goniometer.get_rotation_axis())
  DPS.set_detector(detector)

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  data = flex.vec3_double()
  for spot in spots_mm:
    data.append((spot['xyzobs.mm.value'][0],
                 spot['xyzobs.mm.value'][1],
                 spot['xyzobs.mm.value'][2]*180./math.pi))

  #from matplotlib import pyplot as plt
  #plt.plot([spot.centroid_position[0] for spot in spots_mm] , [spot.centroid_position[1] for spot in spots_mm], 'ro')
  #plt.show()

  DPS.index(raw_spot_input = data, panel_addresses = flex.int([s['panel'] for s in spots_mm]))

  # for development, we want an exhaustive plot of beam probability map:
  params.indexing.plot_search_scope = False

  # perform calculation
  if params.indexing.improve_local_scope=="origin_offset":
    new_detector = DPS.optimize_origin_offset_local_scope()
    return new_detector, beam
  elif params.indexing.improve_local_scope=="S0_vector":
    new_S0_vector = DPS.optimize_S0_local_scope()
    import copy
    new_beam = copy.copy(beam)
    new_beam.set_s0(new_S0_vector)
    return detector, new_beam

  #NKS TO DO: implement rigorous scope instead of local scope.

def candidate_basis_vectors_fft1d(spot_positions, detector, beam,
                                  goniometer, scan, params,
                                  max_cell=None):

  # Spot_positions: Centroid positions for spotfinder spots, in pixels
  # Return value: Corrected for parallax, converted to mm

  spots_mm = Indexer._map_spots_pixel_to_mm_rad(
    spots=spot_positions, detector=detector, scan=scan)

  if len(detector) > 1:
    panel_ids = [spot['panel'] for spot in spot_positions]
  else:
    panel_ids = None


  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell=max_cell,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  DPS.S0_vector = matrix.col(beam.get_s0())
  DPS.inv_wave = 1./beam.get_wavelength()
  if goniometer is None:
    DPS.axis = matrix.col((1,0,0))
  else:
    DPS.axis = matrix.col(goniometer.get_rotation_axis())
  DPS.set_detector(detector)

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  data = flex.vec3_double()
  for spot in spots_mm:
    data.append((spot['xyzobs.mm.value'][0],
                 spot['xyzobs.mm.value'][1],
                 spot['xyzobs.mm.value'][2]*180./math.pi))

  DPS.index(raw_spot_input = data, panel_addresses = panel_ids)
  solutions = DPS.getSolutions()
  return [matrix.col(s.bvec()) for s in solutions],DPS.getXyzData()

def candidate_basis_vectors_fft3d(rs_positions_xyz, params):
  # FIXME implement this
  candidate_basis_vectors_lattice_1 = []
  candidate_basis_vectors = [candidate_basis_vectors_lattice_1]
  return candidate_basis_vectors

def determine_basis_set(candidate_basis_vectors_one_lattice,
                        spot_positions, detector, beam,
                        goniometer, scan,
                        rs_positions_xyz,
                        params,
                        panelID=None
                        ):
  # given a list of 3 or more candidate basis vectors, decide on a good
  # basis set for indexing [triclinic A matrix]

  if panelID is None:
    panelID = flex.int(len(spot_positions),0)

  spots_mm = Indexer._map_spots_pixel_to_mm_rad(
    spots=spot_positions, detector=detector, scan=scan)

  #candidate vectors are direct space vectors with cell lengths in Angstroms
  # Must be presorted by KVAL; KVAL is the height of the FFT peak; specifically
  # |F(kmax)| in Figure 2 of Sauter et al (2004) 37: 399
  assert len(candidate_basis_vectors_one_lattice) > 3
  from rstbx.indexing_api.basis_choice import SelectBasisMetaprocedure as SBM
  from rstbx.indexing_api.lattice import basis_choice_adapter
  BCA = basis_choice_adapter()

  from scitbx import matrix
  BCA.S0_vector = matrix.col(beam.get_s0())
  BCA.inv_wave = 1./beam.get_wavelength()
  if goniometer is None:
    BCA.axis = matrix.col((1,0,0))
  else:
    BCA.axis = matrix.col(goniometer.get_rotation_axis())
  BCA.set_detector(detector)
  BCA.panelID = panelID

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  data = flex.vec3_double()
  for spot in spots_mm:
    data.append((spot['xyzobs.mm.value'][0],
                 spot['xyzobs.mm.value'][1],
                 spot['xyzobs.mm.value'][2]*180./math.pi))
  BCA.raw_spot_input = data

  from rstbx.dps_core import Direction
  directions = [Direction(cvec) for cvec in candidate_basis_vectors_one_lattice]
  flexlist = flex.Direction()
  for ij,d in enumerate(directions):
    d.kval = len(directions) - ij
    flexlist.append(d)
  BCA.set_presorted_solutions(flexlist)
  BCA.setXyzData(rs_positions_xyz)
  M = SBM(input_index_engine = BCA,input_dictionary = {}, horizon_phil = params)
  A_matrix = BCA.getOrientation().reduced_cell()

  triclinic_a_matrix = [A_matrix]
  return triclinic_a_matrix

# this is some kind of bodge to define a container class type

def bravais_lattice_factory(bravais_lattice_type, unit_cell,
                            a_axis, b_axis, c_axis, penalty):
  from collections import namedtuple
  BravaisLattice = namedtuple(
    'bravais_lattice_type unit_cell a_axis b_axis c_axis penalty')
  return BravaisLattice(bravais_lattice_type, unit_cell,
                        a_axis, b_axis, c_axis, penalty)

def possible_lattices_given_basis_set(rs_positions_xyz,
                                      basis_set):
  # given rs_positions_xyz decide list of possible Bravais lattice, unit
  # cell constant combinations from an input triclinic basis

  # FIXME generate list of Bravais lattice types from the
  # bravais_lattice_factory

  return

def refine_lattice_from_raw_spot_positions(indexer, spot_positions, lattice):
  # FIXME interface for doing the refinement
  positions_mm_rad = indexer._map_spots_pixel_to_mm_strategy(spot_positions)
  # do some refinement
  return


# up to here...

class toy_validate_spots_detector(object):
  def __init__(self):
    pass

  def __call__(self, spots, detector):
    '''Make sure that the input is sensible.'''

    # FIXME make this work as follows:
    # - work through spot list; for each spot verify that panel is defined
    #   and that the position is within the bounds of the panel

    return

class Indexer(object):
  ''' The indexer base class. '''

  def __init__(self, strategies, parameters):
    ''' Initialise the indexer base class.

    Params:
    strategies: TBD; strategies for e.g. basis discovery
    parameters: TBD; the set of phil parameters

    '''
    # FIXME make these work properly by naming the strategies etc.
    self.strategies = strategies
    self.parameters = parameters
    # will work like this:
    # self.do_this = do_this
    # however N.B. that this will probably mean you cannot pickle an
    # instance of this

    self._refined = False
    self._analysed = False
    self._indexed = False
    self._setuped = False

    return

  def __call__(self, spots, detector, beam, goniometer = None, scan = None):
    ''' Call to index.

    Params:
    spots: The spot list inc. panel number, space for lattice ID
    detectors: dxtbx detector
    beam: beam information
    goniometer: the goniometer; optional (for e.g. still images)
    scan: the scan information; optional (for e.g. still images)

    Returns:
    TBD

    '''

    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

    # structured loops within loops to employ input strategies to -
    #
    # - validate input (setup) *** this is not a strategy optional ***
    # - discover beam centre (setup)
    # - map spots pixels to mm (setup) *** this will depend on geometry ***
    # - map spots to RS (setup)
    # - determine candidate basis vectors (index)
    # - determine basis sets ([P1_matrix], spots_and_lattice_id) (analyse)
    # - score possible lattice for each solution (analyse)
    # - refine lattice for each solution (refine) *1
    # - reject outliers for each solution (refine)
    #
    # *1 this may also be performed with or without the lattice constraints
    #    so will need as input a BravaisLattice type...

    while not self._refined:
      while not self._analysed:
        while not self._indexed:
          while not self._setuped:
            self._setup()
          self._index()
        self._analyse()
      self._refine()

    return

  def _setup(self):
    # FIXME one day we should implement this
    # self._discover_beam_centre_strategy()
    # FIXME this should probably be delegated to the detector object
    self.spots_mm = self._map_spots_pixel_to_mm_rad(
        self.spots, self.detector, self.scan)
    self._map_centroids_to_reciprocal_space(
        self.spots_mm, self.detector, self.beam, self.goniometer)
    self._setuped = True

    return

  def _index(self):
    self._index_strategy(self)
    # if lots of unindexed reflections:
    #     consider indexing unindexed spots
    #     store results in self._lattices
    self._indexed = True
    return

  def _analyse(self):
    for lattice in self._lattices:
      self._analyse_strategy(lattice)
    self._analysed = True
    return

  def _refine(self):
    for lattice in self._lattices:
      for bravais_lattice in lattice.bravais_lattices:
        self._refine_strategy(bravais_lattice, spots)

    if not self._refined:
      # perhaps need to wind back to the mapping to reciprocal space and
      # try reindexing
      return

    # FIXME may well have a think about which Bravais lattice is correct
    # here

    for lattice in self._lattices:
      # this will perform outlier searching on each lattice and on every
      # spot - we may start going back to indexing at this stage to
      # mop up unindexed reflections
      self._outlier_strategy(lattice, spots)
    self._refined = True
    return

  def set_target_cell_lattice(self, cell, lattice):
    self._indexer_cell = cell
    self._indexer_lattice = lattice
    return

  def set_max_primitive_cell(self, max_primitive_cell):
    self._indexer_max_primitive_cell = max_primitive_cell
    return

  @staticmethod
  def _map_spots_pixel_to_mm_rad(spots, detector, scan):
    from dials.algorithms.centroid import centroid_px_to_mm_panel

    if not spots.has_key('xyzobs.mm.value'):    spots['xyzobs.mm.value']    = flex.vec3_double(len(spots))
    if not spots.has_key('xyzobs.mm.variance'): spots['xyzobs.mm.variance'] = flex.vec3_double(len(spots))
    if not spots.has_key('frame_number'):       spots['frame_number']       = flex.double(len(spots))
    if not spots.has_key('rotation_angle'):     spots['rotation_angle']     = flex.double(len(spots))

    for i_spot, spot in enumerate(spots):
      # just a quick check for now that the reflections haven't come from
      # somewhere else
      #assert spot['xyzobs.mm.value'] == (0,0,0)

      # set reflection properties that might be needed by the dials refinement
      # engine, and convert values from pixels and image number to mm/rads
      spot['frame_number'] = spot['xyzobs.px.value'][2]
      # XXX nasty hack - why are the centroid variances ever zero in the
      # first place?
      centroid_variance = list(spot['xyzobs.px.variance'])
      for i in range(3):
        if centroid_variance[i] == 0:
          centroid_variance[i] = 0.25
      spot['xyzobs.px.variance'] = centroid_variance
      centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
        detector[spot['panel']], scan,
        spot['xyzobs.px.value'],
        spot['xyzobs.px.variance'],
        (1,1,1))
      spot['xyzobs.mm.value'] = centroid_position
      spot['xyzobs.mm.variance'] = centroid_variance
      spot['rotation_angle'] = centroid_position[2]
      spots[i_spot] = spot
    return spots

  @staticmethod
  def _map_centroids_to_reciprocal_space(spots_mm, detector, beam, goniometer):
    assert(len(detector) == 1)
    x, y, _ = spots_mm.centroid_position().parts()
    s1 = detector[0].get_lab_coord(flex.vec2_double(x,y))
    s1 = s1/s1.norms() * (1/beam.get_wavelength())
    spots_mm.set_beam_vector(s1) # needed by refinement
    S = s1 - beam.get_s0()
    reciprocal_space_points = S.rotate_around_origin(
      goniometer.get_rotation_axis(),
      -spots_mm.rotation_angle())
    return reciprocal_space_points


# etc.

class IndexerFactory(object):
  ''' Factory class to create indexers '''

  @staticmethod
  def get_from_somewhere_else(params):
    '''Get a different indexer implementation, which for example may
    overload __call__ internally, or something'''

    # FIXME in here check with the registry for one of these based on the
    # input PHIL parameters

    return False

  @staticmethod
  def from_parameters(params):
    ''' Given a set of parameters, construct the indexer

    Params:
    params The input parameters

    Returns:
    The indexer instance

    '''

    one_from_somewhere_else = IndexerFactory.get_from_somewhere_else(params)
    if one_from_somewhere_else:
      return one_from_somewhere_else

    # else configure a standard one with strategies

    # this is deliberately not implemented
    strategies = IndexerFactory.get_strategies_from_somewhere(params)

    # Return the indexer with the given strategies
    return Indexer(strategies, params)

  @staticmethod
  def get_strategies_from_somewhere(params):
    '''Get the strategies from somewhere, for example a registry.'''

    indexing = params.indexing

    return { } # or whatever; TBD
