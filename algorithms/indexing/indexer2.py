#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.indexer2.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import cPickle as pickle
import math

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError, e:
  pass

import iotbx.phil
from scitbx import matrix

from cctbx.array_family import flex
from cctbx import crystal, sgtbx, xray

from cctbx.crystal.crystal_model import crystal_model as Crystal

import libtbx.load_env
dials_path = libtbx.env.dist_path('dials')

master_phil_scope = iotbx.phil.parse("""
reference {
  detector = None
    .type = path
    .help = "Use detector model from the given reference sweep."
  beam = None
    .type = path
    .help = "Use beam model from the given reference sweep."
}

min_cell = 20
  .type = float(value_min=0)
  .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
max_cell = 160
  .type = float(value_min=0)
  .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
reciprocal_space_grid {
  n_points = 256
    .type = int(value_min=0)
  d_min = 4
    .type = float(value_min=0)
    .help = "The high resolution limit in Angstrom for spots to include in "
            "the initial indexing."
}
b_iso = 200
  .type = float(value_min=0)
rmsd_cutoff = 15
  .type = float(value_min=0)
scan_range = None
  .help = "The range of images to use in indexing. Number of arguments"
          "must be a factor of two. Specifying \"0 0\" will use all images"
          "by default. The given range follows C conventions"
          "(e.g. j0 <= j < j1)."
  .type = ints(size=2)
  .multiple = True
known_symmetry {
  space_group = None
    .type = space_group
  unit_cell = None
    .type = unit_cell
  relative_length_tolerance = 0.1
    .type = float
    .help = "Relative tolerance for unit cell lengths in unit cell comparision."
  absolute_angle_tolerance = 10
    .type = float
    .help = "Angular tolerance (in degrees) in unit cell comparison."
}
optimise_initial_basis_vectors = False
  .type = bool
debug = False
  .type = bool
debug_plots = False
  .type = bool
  .help = "Requires matplotlib"
show_timing = False
  .type = bool
include file %s/data/refinement.phil
refinement_protocol {
  weight_outlier_n_sigma = 5
    .type = float(value_min=0)
  n_macro_cycles = 3
    .type = int(value_min=0)
  d_min_step = 0.5
    .type = float(value_min=0.0)
    .help = "Reduction per step in d_min for reflections to include in refinement."
  d_min_start = 3.0
    .type = float(value_min=0.0)
  d_min_final = None
    .type = float(value_min=0.0)
    .help = "Do not ever include reflections below this value in refinement."
  verbosity = 1
    .type = int(value_min=0)
  outlier_rejection {
    maximum_spot_error = None
      .type = float(value_min=0)
      .help = "Reject reflections whose predicted and observed centroids differ "
              "by more than the given multiple of the pixel size."
              "No outlier rejection is performed in the first macro cycle, and "
              "in the second macro cycle twice the given multiple is used."
    hkl_tolerance = 0.3
      .type = float(value_min=0, value_max=0.5)
  }
}
export_xds_files = False
  .type = bool
  .help = "Export results as XDS.INP, XPARM.XDS for integration with XDS."
multiple_lattice_search = False
  .type = bool
max_lattices = None
  .type = int
method = *fft3d fft1d real_space_grid_search
  .type = choice
cluster_analysis {
  method = *dbscan hcluster
    .type = choice
  hcluster {
    linkage {
      method = *ward
        .type = choice
      metric = *euclidean
        .type = choice
    }
    cutoff = 15
      .type = float(value_min=0)
    cutoff_criterion = *distance inconsistent
      .type = choice
  }
  dbscan {
    eps = 0.05
      .type = float(value_min=0.0)
    min_samples = 30
      .type = int(value_min=1)
  }
  min_cluster_size = 20
    .type = int(value_min=0)
  intersection_union_ratio_cutoff = 0.4
    .type = float(value_min=0.0, value_max=1.0)
}
""" %dials_path, process_includes=True)

master_params = master_phil_scope.fetch().extract()


class vector_group(object):
  def __init__(self):
    self.vectors = []
    self.lengths = []
    self._mean = None

  def append(self, vector, length):
    self.vectors.append(vector)
    self.lengths.append(length)
    self._mean = self.compute_mean()

  def mean(self):
    if self._mean is None:
      self._mean = self.compute_mean()
    return self._mean

  def compute_mean(self):
    sum_x = 0
    sum_y = 0
    sum_z = 0
    for v in self.vectors:
      sum_x += v.elems[0]
      sum_y += v.elems[1]
      sum_z += v.elems[2]
    return matrix.col((sum_x, sum_y, sum_z))/len(self.vectors)


def is_approximate_integer_multiple(vec_a, vec_b,
                                    relative_tolerance=0.2,
                                    angular_tolerance=5.0):
  length_a = vec_a.length()
  length_b = vec_b.length()
  #assert length_b >= length_a
  if length_a > length_b:
    vec_a, vec_b = vec_b, vec_a
    length_a, length_b = length_b, length_a
  angle = vec_a.angle(vec_b, deg=True)
  if angle < angular_tolerance or abs(180-angle) < angular_tolerance:
    n = length_b/length_a
    if abs(round(n) - n) < relative_tolerance:
      return True
  return False


deg_to_radians = math.pi/180


class indexer_base(object):

  def __init__(self, reflections, sweep, params=None):
    # XXX should use ReflectionTable internally instead of ReflectionList
    from dials.model.data import ReflectionList
    reflections = ReflectionList.from_table(reflections)

    self.reflections_raw = reflections
    self.reflections = reflections.deep_copy()
    # the lattice a given reflection belongs to: a value of -1 indicates
    # that a reflection doesn't belong to any lattice so far
    self.reflections_i_lattice = flex.int(self.reflections.size(), -1)
    self.sweep = sweep

    from dxtbx.serialize import load
    if params.reference.detector is not None:
      imageset = load.imageset(params.reference.detector)
      print "Replacing detector:"
      print self.sweep.get_detector()
      print "with:"
      print imageset.get_detector()
      self.sweep.set_detector(imageset.get_detector())
    if params.reference.beam is not None:
      imageset = load.imageset(params.reference.beam)
      print "Replacing beam:"
      print self.sweep.get_beam()
      print "with:"
      print imageset.get_beam()
      self.sweep.set_beam(imageset.get_beam())

    self.goniometer = sweep.get_goniometer()
    self.detector = sweep.get_detector()
    self.scan = sweep.get_scan()
    self.beam = sweep.get_beam()
    if params is None: params = master_params
    self.params = params

    from libtbx.utils import time_log
    self._index_reflections_timer = time_log("index_reflections")
    self._refine_timer = time_log("refinement")
    self._refine_core_timer = time_log("refinement_core")
    self._map_centroids_timer = time_log("map_centroids")
    self._map_to_grid_timer = time_log("map_to_grid")
    self._fft_timer = time_log("fft")
    self._find_peaks_timer = time_log("find_peaks")
    self._cluster_analysis_timer = time_log("cluster_analysis")
    self._ray_intersection_timer = time_log("ray_intersection")

    self.target_symmetry_primitive = None
    self.target_symmetry_centred = None
    if (self.params.known_symmetry.space_group is not None or
        self.params.known_symmetry.unit_cell is not None):
      is_centred = False
      if self.params.known_symmetry.space_group is not None:
        space_group_info = self.params.known_symmetry.space_group
        is_centred = space_group_info.group().conventional_centring_type_symbol() != 'P'
        cb_op_to_primitive = space_group_info.change_of_basis_op_to_primitive_setting()
        sgi_primitive = space_group_info.change_basis(cb_op_to_primitive)
        space_group_primitive = sgi_primitive.group()
      else:
        space_group_primitive = sgtbx.space_group("P 1")
      self.target_symmetry_primitive = crystal.symmetry(
        unit_cell=self.params.known_symmetry.unit_cell,
        space_group=space_group_primitive,
        assert_is_compatible_unit_cell=False)
      if is_centred and self.params.known_symmetry.unit_cell is not None:
        if not self.target_symmetry_primitive.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell):
          self.target_symmetry_centred = crystal.symmetry(
            unit_cell=self.params.known_symmetry.unit_cell,
            space_group=self.params.known_symmetry.space_group.group())
          self.target_symmetry_primitive = self.target_symmetry_centred.change_basis(
            self.target_symmetry_centred.change_of_basis_op_to_primitive_setting())
        else:
          self.target_symmetry_centred = self.target_symmetry_primitive.change_basis(
            self.target_symmetry_primitive.change_of_basis_op_to_reference_setting())
      if self.params.known_symmetry.unit_cell is not None:
        assert (self.target_symmetry_primitive.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell) or self.target_symmetry_centred.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell))

  def index(self):
    self.reflections = self.map_spots_pixel_to_mm_rad(
      self.reflections, self.detector, self.scan)
    self.filter_reflections_by_scan_range()
    self.reciprocal_space_points = self.map_centroids_to_reciprocal_space(
      self.reflections, self.detector, self.beam, self.goniometer)

    if self.params.debug:
      self.debug_write_reciprocal_lattice_points_as_pdb()

    self.reflections.set_crystal(flex.int(self.reflections.size(), -1))

    crystal_models = []

    had_refinement_error = False

    while True:
      if had_refinement_error:
        break
      if self.params.max_lattices is not None and len(crystal_models) >= self.params.max_lattices:
        break
      min_reflections_for_indexing = 40
      if (self.reflections.crystal() == -1).count(True) < min_reflections_for_indexing:
        break

      n_lattices_previous_cycle = len(crystal_models)

      crystal_models.extend(self.find_lattices())
      if len(crystal_models) == n_lattices_previous_cycle:
        # no more lattices found
        break

      self.refined_crystal_models = []

      # for now refine a separate sweep object for each lattice - once we have
      # true multi-lattice refinement we can just refine a single sweep object
      # XXX copy.deepcopy(sweep) does not currently work
      import copy
      sweeps = [copy.deepcopy(self.sweep) for i in range(len(crystal_models))]

      for i_cycle in range(self.params.refinement_protocol.n_macro_cycles):
        if i_cycle > 0:
          self.d_min -= self.params.refinement_protocol.d_min_step
          self.d_min = max(self.d_min, self.params.refinement_protocol.d_min_final)
          if self.d_min < 0:
            break
          print "Increasing resolution to %.1f Angstrom" %self.d_min

        # reset reflection lattice flags
        self.reflections_i_lattice = flex.int(self.reflections.size(), -1)
        self.reflections.set_crystal(self.reflections_i_lattice)

        if i_cycle == 0 and self.target_symmetry_primitive is not None:
          # if a target cell is given make sure that we match any permutation
          # of the cell dimensions
          crystal_models = [self.apply_symmetry(
            crystal_model, self.target_symmetry_primitive, cell_only=True)
                           for crystal_model in crystal_models]

        hkl_tolerance = self.params.refinement_protocol.outlier_rejection.hkl_tolerance
        self.index_reflections(crystal_models, tolerance=hkl_tolerance)

        if (i_cycle == 0 and self.target_symmetry_primitive is not None
            and self.target_symmetry_primitive.space_group() is not None):
          # now apply the space group symmetry only after the first indexing
          # need to make sure that the symmetrized orientation is similar to the P1 model
          crystal_models = [self.apply_symmetry(
            crystal_model, self.target_symmetry_primitive)
                           for crystal_model in crystal_models]

        print
        print "#" * 80
        print "Starting refinement (macro-cycle %i)" %(i_cycle+1)
        print "#" * 80
        print

        self.refined_reflections = []

        for i_lattice, crystal_model in enumerate(crystal_models):
          self.sweep = sweeps[i_lattice] # XXX


          self.i_lattice = i_lattice

          print
          print "Starting refinement of crystal model %i" %(i_lattice+1)
          print "Starting crystal model:"
          print crystal_model
          print

          self.reflections_i_lattice = self.reflections.crystal()
          self.indexed_reflections = (self.reflections_i_lattice == i_lattice)

          if self.params.debug:
            sel = flex.bool(self.reflections.size(), False)
            lengths = 1/self.reciprocal_space_points.norms()
            isel = (lengths >= self.d_min).iselection()
            sel.set_selected(isel, True)
            sel.set_selected(self.reflections_i_lattice > -1, False)
            unindexed = self.reflections_raw.select(sel)
            with open("unindexed.pickle", 'wb') as f:
              pickle.dump(unindexed.to_table(), f)

          maximum_spot_error \
            = self.params.refinement_protocol.outlier_rejection.maximum_spot_error
          if i_cycle == 0:
            maximum_spot_error = None
          elif i_cycle == 1:
            if maximum_spot_error is not None:
              maximum_spot_error *= 2

          try:
            from dials.model.experiment.experiment_list \
                 import Experiment, ExperimentList
            experiments = ExperimentList([Experiment(
              imageset=self.sweep,
              beam=self.sweep.get_beam(),
              detector=self.sweep.get_detector(),
              scan=self.sweep.get_scan(),
              goniometer=self.sweep.get_goniometer(),
              crystal=crystal_model)])
            refined_experiments, refined_reflections = self.refine(
              experiments, maximum_spot_error=maximum_spot_error)
            crystal_model = refined_experiments.crystals()[0]
          except RuntimeError, e:
            s = str(e)
            if "below the configured limit" in s:
              had_refinement_error = True
              print "Refinement failed:"
              print s
              del crystal_models[i_lattice]
              break
            raise

          self.refined_reflections.append(refined_reflections)
          crystal_models[i_lattice] = crystal_model

          self.sweep.set_detector(refined_experiments[0].detector)
          self.sweep.set_beam(refined_experiments[0].beam)
          self.sweep.set_goniometer(refined_experiments[0].goniometer)
          self.sweep.set_scan(refined_experiments[0].scan)
          # these may have been updated in refinement
          # XXX once david has implemented multi-lattice refinement there
          # should be one and only one sweep object to refine
          #self.sweep = sweeps[0]
          self.detector = self.sweep.get_detector()
          self.beam = self.sweep.get_beam()
          self.goniometer = self.sweep.get_goniometer()
          self.scan = self.sweep.get_scan()

          if not (self.params.refinement.parameterisation.beam.fix == 'all'
                  and self.params.refinement.parameterisation.detector.fix == 'all'):
            # Experimental geometry may have changed - re-map centroids to
            # reciprocal space
            self.reciprocal_space_points = self.map_centroids_to_reciprocal_space(
              self.reflections, self.detector, self.beam, self.goniometer)

          if self.d_min == self.params.refinement_protocol.d_min_final:
            print "Target d_min_final reached: finished with refinement"
            break

      if not (self.params.multiple_lattice_search and
              self.params.method == "real_space_grid_search"):
        break

    for i_lattice, crystal_model in enumerate(crystal_models):
      self.refined_crystal_models.append(crystal_model)
      suffix = ""
      if len(crystal_models) > 1:
        suffix = "_%i" %(i_lattice+1)
      self.export_as_json(crystal_model, sweeps[i_lattice], suffix=suffix)
      if self.params.export_xds_files:
        self.export_xds_files(crystal_model, sweeps[i_lattice], suffix=suffix)
      self.export_reflections(
        #self.reflections_raw.select(self.reflections_i_lattice == i_lattice),
        self.refined_reflections[i_lattice],
        file_name='indexed%s.pickle' %suffix)

    if 1 and self.params.debug and self.goniometer is not None:
      for i_lattice, cm in enumerate(self.refined_crystal_models):
        suffix = ""
        if len(crystal_models) > 1:
          suffix = "_%i" %(i_lattice+1)
        self.predict_reflections(cm)
        self.export_predicted_reflections(file_name='predictions%s.pickle' %suffix)

    print "Final refined crystal models:"
    for i, crystal_model in enumerate(self.refined_crystal_models):
      print "model %i (%i reflections):" %(
        i+1, (self.reflections_i_lattice == i).count(True))
      print crystal_model

    if self.params.show_timing:
      print self._index_reflections_timer.legend
      print self._map_centroids_timer.report()
      print self._map_to_grid_timer.report()
      print self._fft_timer.report()
      print self._find_peaks_timer.report()
      print self._cluster_analysis_timer.report()
      print self._index_reflections_timer.report()
      print self._refine_timer.report()
      print self._refine_core_timer.report()
      print self._ray_intersection_timer.report()

  def filter_reflections_by_scan_range(self):
    reflections_in_scan_range = flex.size_t()
    for i_ref, refl in enumerate(self.reflections):
      frame_number = refl.frame_number

      if len(self.params.scan_range):
        reflections_in_range = False
        for scan_range in self.params.scan_range:
          if scan_range is None: continue
          range_start, range_end = scan_range
          if frame_number >= range_start and frame_number < range_end:
            reflections_in_range = True
            break
        if not reflections_in_range:
          continue
      reflections_in_scan_range.append(i_ref)
    self.reflections = self.reflections.select(reflections_in_scan_range)
    self.reflections_raw = self.reflections_raw.select(reflections_in_scan_range)

  @staticmethod
  def map_spots_pixel_to_mm_rad(spots, detector, scan):
    """Reflections that come from dials.spotfinder only have the centroid
    position and variance set, """

    from dials.algorithms.centroid import centroid_px_to_mm_panel
    # ideally don't copy, but have separate spot attributes for mm and pixel
    spots_mm = spots.deep_copy()
    for i_spot, spot in enumerate(spots_mm):
      # just a quick check for now that the reflections haven't come from
      # somewhere else
      assert spot.image_coord_mm == (0,0)

      # set reflection properties that might be needed by the dials refinement
      # engine, and convert values from pixels and image number to mm/rads
      spot.frame_number = spot.centroid_position[2]
      centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
        detector[spot.panel_number], scan,
        spot.centroid_position,
        spot.centroid_variance,
        (1,1,1))
      spot.centroid_position = centroid_position
      spot.centroid_variance = centroid_variance
      spot.rotation_angle = centroid_position[2]
    return spots_mm

  @staticmethod
  def map_centroids_to_reciprocal_space(spots_mm, detector, beam, goniometer):
    panel_numbers = flex.size_t(spot.panel_number for spot in spots_mm)
    reciprocal_space_points = flex.vec3_double()
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      isel = sel.iselection()
      spots_panel = spots_mm.select(panel_numbers == i_panel)
      x, y, _ = spots_panel.centroid_position().parts()
      s1 = detector[i_panel].get_lab_coord(flex.vec2_double(x,y))
      s1 = s1/s1.norms() * (1/beam.get_wavelength())
      for i in range(len(s1)):
        spots_mm[isel[i]].beam_vector = s1[i] # needed by refinement
      #spots_panel.set_beam_vector(s1) # needed by refinement
      S = s1 - beam.get_s0()
      # XXX what about if goniometer fixed rotation is not identity?
      if goniometer is not None:
        reciprocal_space_points.extend(S.rotate_around_origin(
          goniometer.get_rotation_axis(),
          -spots_panel.rotation_angle()))
      else:
        reciprocal_space_points.extend(S)
    return reciprocal_space_points

  def find_candidate_orientation_matrices(self, candidate_basis_vectors,
                                          return_first=False,
                                          apply_symmetry=True):
    candidate_crystal_models = []
    vectors = candidate_basis_vectors

    min_angle = 20 # degrees, arbitrary cutoff
    for i in range(len(vectors)):
      a = vectors[i]
      for j in range(i, len(vectors)):
        b = vectors[j]
        angle = a.angle(b, deg=True)
        if angle < min_angle or (180-angle) < min_angle:
          continue
        a_cross_b = a.cross(b)
        gamma = a.angle(b, deg=True)
        if gamma < 90:
          # all angles obtuse if possible please
          b = -b
          gamma = 180 - gamma
          a_cross_b = -a_cross_b
        for k in range(j, len(vectors)):
          c = vectors[k]
          if abs(90-a_cross_b.angle(c, deg=True)) < min_angle:
            continue
          alpha = b.angle(c, deg=True)
          if alpha < 90:
            c = -c
          beta = c.angle(a, deg=True)
          if a_cross_b.dot(c) < 0:
            # we want right-handed basis set, therefore invert all vectors
            a = -a
            b = -b
            c = -c
            #assert a.cross(b).dot(c) > 0
          model = Crystal(a, b, c, space_group_symbol="P 1")
          uc = model.get_unit_cell()
          if self.target_symmetry_primitive is not None:
            symmetrized_model = self.apply_symmetry(
              model, self.target_symmetry_primitive)
            if symmetrized_model is None:
              if self.target_symmetry_centred is not None:
                symmetrized_model = self.apply_symmetry(
                model, self.target_symmetry_centred,
                return_primitive_setting=True)
            if symmetrized_model is None:
              continue
            if apply_symmetry:
              model = symmetrized_model
              uc = model.get_unit_cell()

          params = uc.parameters()
          if uc.volume() > (params[0]*params[1]*params[2]/100):
            # unit cell volume cutoff from labelit 2004 paper
            candidate_crystal_models.append(model)
            if return_first:
              return candidate_crystal_models
    return candidate_crystal_models

  def apply_symmetry(self, crystal_model, target_symmetry,
                     return_primitive_setting=False,
                     cell_only=False):
    unit_cell = crystal_model.get_unit_cell()
    target_unit_cell = target_symmetry.unit_cell()
    if cell_only:
      target_space_group = sgtbx.space_group('P 1')
    else:
      target_space_group = target_symmetry.space_group()
    A = crystal_model.get_A()
    A_inv = A.inverse()
    unit_cell_is_similar = False
    real_space_a = A_inv[:3]
    real_space_b = A_inv[3:6]
    real_space_c = A_inv[6:9]
    basis_vectors = [real_space_a, real_space_b, real_space_c]
    min_bmsd = 1e8
    best_perm = None
    # for non-cyclic permutations one axis needs inverting to keep system right-handed
    for perm, sign in zip(
        ((0,1,2), (1,2,0), (2,0,1), (1,0,2), (0,2,1), (2,1,0)),
        (1, 1, 1, -1, -1, -1)):
      crystal_model = Crystal(
        basis_vectors[perm[0]],
        basis_vectors[perm[1]],
        [i * sign for i in basis_vectors[perm[2]]],
        space_group=target_space_group)
      unit_cell = crystal_model.get_unit_cell()
      uc = target_unit_cell
      if uc is None:
        uc = unit_cell
      # XXX what about permuting the target_unit_cell (if not None)?
      symm_target_sg = crystal.symmetry(
        unit_cell=uc,
        space_group=target_space_group,
        assert_is_compatible_unit_cell=False)
      # this assumes that the initial basis vectors are good enough that
      # we can tell which should be the unique axis - probably not a robust
      # solution
      if unit_cell.is_similar_to(
        symm_target_sg.unit_cell(),
        relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
        absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance):
        bmsd = unit_cell.bases_mean_square_difference(
          symm_target_sg.unit_cell())
        eps = 1e-8
        if (bmsd+eps) < min_bmsd:
          min_bmsd = bmsd
          best_perm = list(perm)
    if best_perm is None:
      return None
    crystal_model = Crystal(
      basis_vectors[best_perm[0]],
      basis_vectors[best_perm[1]],
      basis_vectors[best_perm[2]],
      space_group=target_space_group)
    model = crystal_model
    #cb_op_target_ref = symm_target_sg.space_group_info().type().cb_op()
    #symm_target_sg_ref = symm_target_sg.change_basis(cb_op_target_ref)
    from rstbx.symmetry.constraints import parameter_reduction
    s = parameter_reduction.symmetrize_reduce_enlarge(target_space_group)
    s.set_orientation(crystal_model.get_A())
    s.symmetrize()
    #direct_matrix = s.orientation.change_basis(cb_op_target_ref).direct_matrix()
    if return_primitive_setting:
      sgi = sgtbx.space_group_info(group=target_space_group)
      cb_op_to_primitive = sgi.change_of_basis_op_to_primitive_setting()
      direct_matrix = s.orientation.change_basis(
        cb_op_to_primitive).direct_matrix()
    else:
      direct_matrix = s.orientation.direct_matrix()
    a = matrix.col(direct_matrix[:3])
    b = matrix.col(direct_matrix[3:6])
    c = matrix.col(direct_matrix[6:9])
    ## verify it is still right-handed basis set
    #assert a.cross(b).dot(c) > 0
    model = Crystal(
      a, b, c, space_group=target_space_group)
    return model

  def index_reflections_given_orientation_matix(
      self, crystal_model, tolerance=0.2, verbose=0):

    self._index_reflections_timer.start()

    if verbose > 1:
      print "Candidate crystal model:"
      print crystal_model

    n_rejects = 0

    miller_indices = flex.miller_index()
    indexed_reflections = flex.size_t()

    A = crystal_model.get_A()
    A_inv = A.inverse()

    d_spacings = 1/self.reciprocal_space_points.norms()
    inside_resolution_limit = d_spacings > self.d_min
    sel = inside_resolution_limit & (self.reflections_i_lattice == -1)
    isel = sel.iselection()
    rlps = self.reciprocal_space_points.select(isel)
    hkl_float = tuple(A_inv) * rlps
    hkl_int = hkl_float.iround()

    for i_hkl in range(hkl_int.size()):
      max_difference = max([abs(hkl_float[i_hkl][i] - hkl_int[i_hkl][i]) for i in range(3)])
      if max_difference > tolerance:
        n_rejects += 1
        continue
      miller_index = hkl_int[i_hkl]
      miller_indices.append(miller_index)
      i_ref = isel[i_hkl]
      self.reflections[i_ref].miller_index = miller_index
      indexed_reflections.append(i_ref)

    self._index_reflections_timer.stop()

    return indexed_reflections, miller_indices

  def index_reflections(self, crystal_models, tolerance=0.3):
    self._index_reflections_timer.start()
    from dials.algorithms.indexing import index_reflections
    index_reflections(self.reflections, self.reciprocal_space_points,
                      crystal_models, self.d_min, tolerance=tolerance,
                      verbosity=self.params.refinement_protocol.verbosity)
    self._index_reflections_timer.stop()

  def refine(self, experiments, maximum_spot_error=None):
    self._refine_timer.start()
    from dials.algorithms.indexing.refinement import refine
    reflections_for_refinement = self.reflections.select(
      self.indexed_reflections)
    refiner, refined = refine(
      self.params, reflections_for_refinement, experiments,
      maximum_spot_error=maximum_spot_error,
      verbosity=self.params.refinement_protocol.verbosity,
      debug_plots=self.params.debug_plots)
    used_reflections = refiner.get_reflections()
    verbosity = self.params.refinement_protocol.verbosity

    self._refine_timer.stop()
    return refiner.get_experiments(), used_reflections

  def predict_reflections(self, crystal_model):
    from dials.algorithms.integration import ReflectionPredictor
    predictor = ReflectionPredictor()

    sigma_divergence = self.beam.get_sigma_divergence()
    mosaicity = crystal_model.get_mosaicity()

    if sigma_divergence == 0.0:
      self.beam.set_sigma_divergence(0.02) # degrees
    if mosaicity == 0.0:
      crystal_model.set_mosaicity(0.139) # degrees

    reflections = predictor(self.sweep, crystal_model)
    self.predicted_reflections = reflections
    return self.predicted_reflections

  def export_predicted_reflections(self, file_name='predictions.pickle'):
    from dials.model.serialize import dump
    dump.reflections(self.predicted_reflections.to_table(), file_name)

  def debug_show_candidate_basis_vectors(self):

    vectors = self.candidate_basis_vectors

    for i, v in enumerate(vectors):
      print i, v.length()# , vector_heights[i]

    # print a table of the angles between each pair of vectors

    angles = flex.double(len(vectors)**2)
    angles.reshape(flex.grid(len(vectors), len(vectors)))

    for i in range(len(vectors)):
      v_i = vectors[i]
      for j in range(i+1, len(vectors)):
        v_j = vectors[j]
        angles[i,j] = v_i.angle(v_j, deg=True)

    print (" "*7),
    for i in range(len(vectors)):
      print "%7.3f" % vectors[i].length(),
    print
    for i in range(len(vectors)):
      print "%7.3f" % vectors[i].length(),
      for j in range(len(vectors)):
        if j <= i:
          print (" "*7),
        else:
          print "%5.1f  " %angles[i,j],
      print

  def debug_write_reciprocal_lattice_points_as_pdb(
      self, file_name='reciprocal_lattice.pdb'):
    from cctbx import crystal, xray
    cs = crystal.symmetry(unit_cell=(1000,1000,1000,90,90,90), space_group="P1")
    xs = xray.structure(crystal_symmetry=cs)
    for site in self.reciprocal_space_points:
      xs.add_scatterer(xray.scatterer("C", site=site))

    xs.sites_mod_short()
    with open(file_name, 'wb') as f:
      print >> f, xs.as_pdb_file()

  def debug_write_ccp4_map(self, map_data, file_name):
    from iotbx import ccp4_map
    gridding_first = (0,0,0)
    gridding_last = map_data.all()
    labels = ["cctbx.miller.fft_map"]
    ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=self.fft_cell,
      space_group=sgtbx.space_group("P1"),
      gridding_first=gridding_first,
      gridding_last=gridding_last,
      map_data=map_data,
      labels=flex.std_string(labels))

  def export_as_json(self, crystal_model, sweep, suffix=None, compact=False):
    from dials.model.serialize import dump
    from dials.model.experiment.experiment_list import ExperimentList, Experiment
    experiments = ExperimentList()
    experiments.append(Experiment(
      imageset=sweep, beam=sweep.get_beam(), detector=sweep.get_detector(),
      goniometer=sweep.get_goniometer(), scan=sweep.get_scan(),
      crystal=crystal_model))
    assert experiments.is_consistent()
    dump.experiment_list(experiments, 'experiments%s.json' %suffix)

  def export_reflections(self, reflections, file_name="reflections.pickle"):
    with open(file_name, 'wb') as f:
      pickle.dump(reflections.to_table(), f)

  def export_xds_files(self, crystal_model, sweep, suffix=None):
    from dxtbx.serialize import xds
    if suffix is None:
      suffix = ''
    crystal_model = crystal_model.change_basis(
      crystal_model.get_space_group().info().change_of_basis_op_to_reference_setting())
    A = crystal_model.get_A()
    A_inv = A.inverse()
    real_space_a = A_inv.elems[:3]
    real_space_b = A_inv.elems[3:6]
    real_space_c = A_inv.elems[6:9]
    to_xds = xds.to_xds(sweep)
    with open('XDS%s.INP' %suffix, 'wb') as f:
      to_xds.XDS_INP(out=f, job_card="XYCORR INIT DEFPIX INTEGRATE CORRECT")
    with open('XPARM%s.XDS' %suffix, 'wb') as f:
      to_xds.xparm_xds(
        real_space_a, real_space_b, real_space_c,
        crystal_model.get_space_group().type().number(),
        out=f)

  def find_lattice(self):
    raise NotImplementedError()

  def find_candidate_basis_vectors_nks(self, vectors):
    '''Find the candidate basis vectors from the Patterson peaks using code
    from NKS which will search for the basis which best describes the list of
    input spot positions. Based on using indexer.determine_basis_set.'''

    from dials.algorithms.indexing import indexer

    # hmm... conventionally the basis selection works on around 30 possible
    # basis vectors, the code above appears to generate 200+ for one example

    from rstbx.phil.phil_preferences import libtbx_defs,iotbx_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
      input_string=iotbx_defs+libtbx_defs).extract()

    # do we really need all of these parameters? surely some of them seem
    # a little ... redundant

    triclinic_crystal = indexer.determine_basis_set(
      candidate_basis_vectors_one_lattice=vectors,
      spot_positions=self.reflections_raw,
      detector=self.detector,
      beam=self.beam,
      goniometer=self.goniometer,
      scan=self.scan,
      rs_positions_xyz=self.reciprocal_space_points,
      params=hardcoded_phil)

    direct_matrix = triclinic_crystal[0].direct_matrix()
    return [Crystal(direct_matrix[0:3],
                    direct_matrix[3:6],
                    direct_matrix[6:9],
                    space_group_symbol="P 1")]






def optimise_basis_vectors(reciprocal_space_points, vectors):
  optimised = flex.vec3_double()
  for vector in vectors:
    minimised = basis_vector_minimser(reciprocal_space_points, vector)
    optimised.append(tuple(minimised.x))
  return optimised


from scitbx import lbfgs

# Optimise the initial basis vectors as per equation 11.4.3.4 of
# Otwinowski et al, International Tables Vol. F, chapter 11.4 pp. 282-295
class basis_vector_target(object):

  def __init__(self, reciprocal_space_points):
    self.reciprocal_space_points = reciprocal_space_points
    self._xyz_parts = self.reciprocal_space_points.parts()

  def compute_functional_and_gradients(self, vector):
    assert len(vector) == 3
    two_pi_S_dot_v = 2 * math.pi * self.reciprocal_space_points.dot(vector)
    f = - flex.sum(flex.cos(two_pi_S_dot_v))
    sin_part = flex.sin(two_pi_S_dot_v)
    g = flex.double([flex.sum(2 * math.pi * self._xyz_parts[i] * sin_part)
                     for i in range(3)])
    return f, g


class basis_vector_minimser(object):
  def __init__(self, reciprocal_space_points, vector,
               lbfgs_termination_params=None,
               lbfgs_core_params=lbfgs.core_parameters(m=20)):
    self.reciprocal_space_points = reciprocal_space_points
    if not isinstance(vector, flex.double):
      self.x = flex.double(vector)
    else:
      self.x = vector.deep_copy()
    self.n = len(self.x)
    assert self.n == 3
    self.target = basis_vector_target(self.reciprocal_space_points)
    self.minimizer = lbfgs.run(target_evaluator=self,
                               termination_params=lbfgs_termination_params,
                               core_params=lbfgs_core_params)
    #print "number of iterations:", self.minimizer.iter()

  def compute_functional_and_gradients(self):
    f, g = self.target.compute_functional_and_gradients(tuple(self.x))
    #g_fd = _gradient_fd(self.target, tuple(self.x))
    #from libtbx.test_utils import approx_equal
    #assert approx_equal(g, g_fd, eps=1e-3)
    return f, g

  def callback_after_step(self, minimizer):
    #print tuple(self.x)
    return


def _gradient_fd(target, vector, eps=1e-6):
  grads = []
  for i in range(len(vector)):
    v = list(vector)
    v[i] -= eps
    tm, _ = target.compute_functional_and_gradients(v)
    v[i] += 2 * eps
    tp, _ = target.compute_functional_and_gradients(v)
    grads.append((tp-tm)/(2*eps))
  return grads


def reject_weight_outliers_selection(reflections, sigma_cutoff=5):
  from scitbx.math import basic_statistics
  variances = flex.vec3_double([r.centroid_variance for r in reflections])
  selection = None
  for v in variances.parts():
    w = 1/v
    ln_w = flex.log(w)
    stats = basic_statistics(ln_w)
    sel = ln_w < (sigma_cutoff * stats.bias_corrected_standard_deviation
                  + stats.mean)
    if selection is None:
      selection = sel
    else:
      selection &= sel
  return selection


def hist_outline(hist):

  step_size = hist.slot_width()
  half_step_size = 0.5 * step_size
  n_slots = len(hist.slots())

  bins = flex.double(n_slots * 2 + 2, 0)
  data = flex.double(n_slots * 2 + 2, 0)
  for i in range(n_slots):
    bins[2 * i + 1] = hist.slot_centers()[i] - half_step_size
    bins[2 * i + 2] = hist.slot_centers()[i] + half_step_size
    data[2 * i + 1] = hist.slots()[i]
    data[2 * i + 2] = hist.slots()[i]

  bins[0] = bins[1] - step_size
  bins[-1] = bins[-2] + step_size
  data[0] = 0
  data[-1] = 0

  return (bins, data)


def plot_centroid_weights_histograms(reflections, n_slots=50):
  from matplotlib import pyplot
  variances = flex.vec3_double([r.centroid_variance for r in reflections])
  vx, vy, vz = variances.parts()
  wx = 1/vx
  wy = 1/vy
  wz = 1/vz
  #hx = flex.histogram(vx, n_slots=n_slots)
  #hy = flex.histogram(vy, n_slots=n_slots)
  #hz = flex.histogram(vz, n_slots=n_slots)
  wx = flex.log(wx)
  wy = flex.log(wy)
  wz = flex.log(wz)
  hx = flex.histogram(wx, n_slots=n_slots)
  hy = flex.histogram(wy, n_slots=n_slots)
  hz = flex.histogram(wz, n_slots=n_slots)
  fig = pyplot.figure()

  #outliers = reflections.select(wx > 50)
  #for refl in outliers:
    #print refl

  for i, h in enumerate([hx, hy, hz]):
    ax = fig.add_subplot(311+i)

    slots = h.slots().as_double()
    bins, data = hist_outline(h)
    log_scale = True
    if log_scale:
      data.set_selected(data == 0, 0.1) # otherwise lines don't get drawn when we have some empty bins
      ax.set_yscale("log")
    ax.plot(bins, data, '-k', linewidth=2)
    #pyplot.suptitle(title)
    data_min = min([slot.low_cutoff for slot in h.slot_infos() if slot.n > 0])
    data_max = max([slot.low_cutoff for slot in h.slot_infos() if slot.n > 0])
    ax.set_xlim(data_min, data_max+h.slot_width())
  pyplot.show()

