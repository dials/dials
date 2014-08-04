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
from libtbx import easy_pickle
import math

from libtbx.utils import Sorry
import iotbx.phil
from scitbx import matrix

from dials.array_family import flex
from cctbx import crystal, sgtbx, uctbx, xray

from dxtbx.model.crystal import crystal_model as Crystal
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList

master_phil_scope = iotbx.phil.parse("""
reference {
  detector = None
    .type = path
    .help = "Use detector model from the given reference sweep."
  beam = None
    .type = path
    .help = "Use beam model from the given reference sweep."
}
discover_better_experimental_model = False
  .type = bool
min_cell = 20
  .type = float(value_min=0)
  .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
max_cell = Auto
  .type = float(value_min=0)
  .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
fft3d {
  peak_search = *flood_fill clean
    .type = choice
  reciprocal_space_grid {
    n_points = 256
      .type = int(value_min=0)
    d_max = None
      .type = float(value_min=0)
    d_min = 4
      .type = float(value_min=0)
      .help = "The high resolution limit in Angstrom for spots to include in "
              "the initial indexing."
  }
}
sigma_phi_deg = None
  .type = float(value_min=0)
  .help = "Override the phi sigmas for refinement. Mainly intended for single-shot"
          "rotation images where the phi sigma is almost certainly incorrect."
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
basis_vector_combinations {
  max_try = 10
    .type = int(value_min=1)
  metric = *model_likelihood n_indexed
    .type = choice
}
index_assignment {
  method = *simple local
    .type = choice
  simple {
    hkl_tolerance = 0.3
      .type = float(value_min=0, value_max=0.5)
  }
  local {
    epsilon = 0.05
      .type = float
    delta = 8
      .type = int
    l_min = 0.8
      .type = float
    nearest_neighbours = 20
      .type = int(value_min=1)
  }
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
include scope dials.data.refinement.phil_scope
refinement_protocol {
  weight_outlier_n_sigma = 5
    .type = float(value_min=0)
  n_macro_cycles = 3
    .type = int(value_min=0)
  d_min_step = 1.0
    .type = float(value_min=0.0)
    .help = "Reduction per step in d_min for reflections to include in refinement."
  d_min_start = 4.0
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
    maximum_phi_error = None
      .type = float(value_min=0)
      .help = "Reject reflections whose predicted and observed phi centroids "
              "differ by more than the given value (degress)."
              "No outlier rejection is performed in the first macro cycle, and "
              "in the second macro cycle twice the given multiple is used."
  }
}
method = *fft3d fft1d real_space_grid_search
  .type = choice
multiple_lattice_search {
  cluster_analysis_search = False
    .type = bool
    .help = "Perform cluster analysis search for multiple lattices."
  recycle_unindexed_reflections = False
    .type = bool
    .help = "Attempt another cycle of indexing on the unindexed reflections "
            "if enough reflections are unindexed."
  recycle_unindexed_reflections_cutoff = 0.1
    .type = float(value_min=0, value_max=1)
    .help = "Attempt another cycle of indexing on the unindexed reflections "
            "if more than the fraction of input reflections are unindexed."
  minimum_angular_separation = 5
    .type = float(value_min=0)
    .help = "The minimum angular separation (in degrees) between two lattices."
  max_lattices = None
    .type = int
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
}
output {
  experiments = experiments.json
    .type = path
  reflections = indexed.pickle
    .type = path
}

""", process_includes=True)

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

  def __init__(self, reflections, imagesets, params=None):
    # XXX should use ReflectionTable internally instead of ReflectionList
    from dials.model.data import ReflectionList

    self.reflections = reflections
    self.imagesets = imagesets
    for imageset in imagesets[1:]:
      assert imageset.get_detector() == self.imagesets[0].get_detector()
      assert imageset.get_beam() == self.imagesets[0].get_beam()
      assert imageset.get_goniometer() == self.imagesets[0].get_goniometer()
    sweep = self.imagesets[0]
    self.sweep = sweep

    from dxtbx.serialize import load
    if params is None: params = master_params
    if params.reference.detector is not None:
      try:
        experiments = load.experiment_list(
          params.reference.detector, check_format=False)
        assert len(experiments.detectors()) == 1
        reference_detector = experiments.detectors()[0]
      except Exception, e:
        imageset = load.imageset(params.reference.detector)
        reference_detector = imageset.get_detector()
      print "Replacing detector:"
      print self.sweep.get_detector()
      print "with:"
      print reference_detector
      for imageset in self.imagesets:
        imageset.set_detector(reference_detector)
    if params.reference.beam is not None:
      try:
        experiments = load.experiment_list(
          params.reference.detector, check_format=False)
        assert len(experiments.beams()) == 1
        reference_beam = experiments.beams()[0]
      except Exception, e:
        imageset = load.imageset(params.reference.beam)
        reference_beam = imageset.get_beam()
      print "Replacing beam:"
      print self.sweep.get_beam()
      print "with:"
      print reference_beam
      for imageset in self.imagesets:
        imageset.set_beam(reference_beam)

    self.goniometer = sweep.get_goniometer()
    self.detector = sweep.get_detector()
    #self.scan = sweep.get_scan()
    self.beam = sweep.get_beam()
    self.params = params
    self.refined_experiments = None

    self.target_symmetry_primitive = None
    self.target_symmetry_minimum_cell = None
    self.target_symmetry_centred = None
    self.cb_op_centred_to_primitive = None
    self.cb_op_minimum_cell_to_primitive = None
    self.cb_op_primitive_to_given = None
    if (self.params.known_symmetry.space_group is not None or
        self.params.known_symmetry.unit_cell is not None):
      is_centred = False
      if self.params.known_symmetry.space_group is not None:
        space_group_info = self.params.known_symmetry.space_group
        is_centred = space_group_info.group().conventional_centring_type_symbol() != 'P'
        cb_op_given_to_primitive \
          = space_group_info.change_of_basis_op_to_primitive_setting()
        sgi_primitive = space_group_info.change_basis(
          cb_op_given_to_primitive)
        self.cb_op_primitive_to_given = cb_op_given_to_primitive.inverse()
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
          self.cb_op_centred_to_primitive \
            = self.target_symmetry_centred.change_of_basis_op_to_primitive_setting()
          self.target_symmetry_primitive = self.target_symmetry_centred.change_basis(
            self.cb_op_centred_to_primitive)
          self.cb_op_primitive_to_given = self.cb_op_centred_to_primitive.inverse()
        else:
          self.target_symmetry_centred = self.target_symmetry_primitive.change_basis(
            self.target_symmetry_primitive.change_of_basis_op_to_reference_setting())
          self.cb_op_centred_to_primitive \
            = self.target_symmetry_centred.change_of_basis_op_to_primitive_setting()
          self.cb_op_primitive_to_given = None # given symmetry was primitive
      if self.params.known_symmetry.unit_cell is not None:
        assert (self.target_symmetry_primitive.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell) or self.target_symmetry_centred.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell))
        self.target_symmetry_minimum_cell = self.target_symmetry_primitive.minimum_cell()
        self.cb_op_minimum_cell_to_primitive \
          = self.target_symmetry_primitive.change_of_basis_op_to_minimum_cell().inverse()
        if self.cb_op_primitive_to_given is not None:
          assert self.target_symmetry_primitive.change_basis(
            self.cb_op_primitive_to_given).unit_cell().is_similar_to(
              self.params.known_symmetry.unit_cell)

  def index(self):
    import libtbx
    self.reflections_input = self.reflections
    self.reflections = flex.reflection_table()
    for i, imageset in enumerate(self.imagesets):
      sel = (self.reflections_input['id'] == i)
      self.reflections.extend(self.map_spots_pixel_to_mm_rad(
        self.reflections_input.select(sel),
        imageset.get_detector(), imageset.get_scan()))
    self.filter_reflections_by_scan_range()
    self.map_centroids_to_reciprocal_space(
      self.reflections, self.detector, self.beam, self.goniometer)

    self.find_max_cell()

    if self.params.sigma_phi_deg is not None:
      var_x, var_y, _ = self.reflections['xyzobs.mm.variance'].parts()
      var_phi_rad = flex.double(
        var_x.size(), (math.pi/180 * self.params.sigma_phi_deg)**2)
      self.reflections['xyzobs.mm.variance'] = flex.vec3_double(
        var_x, var_y, var_phi_rad)

    if self.params.debug:
      self.debug_write_reciprocal_lattice_points_as_pdb()

    self.reflections['id'] = flex.int(len(self.reflections), -1)

    experiments = ExperimentList()

    had_refinement_error = False
    have_similar_crystal_models = False

    if self.params.discover_better_experimental_model:
      opt_detector, opt_beam = self.discover_better_experimental_model(
        self.reflections, self.detector, self.beam, self.goniometer,
        self.imagesets[0].get_scan())
      self.sweep.set_detector(opt_detector)
      self.sweep.set_beam(opt_beam)
      self.detector = opt_detector
      self.beam = opt_beam

    while True:
      self.d_min = self.params.refinement_protocol.d_min_start
      if had_refinement_error or have_similar_crystal_models:
        break
      max_lattices = self.params.multiple_lattice_search.max_lattices
      if max_lattices is not None and len(experiments) >= max_lattices:
        break
      cutoff_fraction = \
        self.params.multiple_lattice_search.recycle_unindexed_reflections_cutoff
      d_spacings = 1/self.reflections['rlp'].norms()
      min_reflections_for_indexing = \
        cutoff_fraction * len(self.reflections.select(d_spacings > self.d_min))
      crystal_ids = self.reflections.select(d_spacings > self.d_min)['id']
      if (crystal_ids == -1).count(True) < min_reflections_for_indexing:
        print "Finish searching for more lattices: %i unindexed reflections remaining." %(
          min_reflections_for_indexing)
        break

      n_lattices_previous_cycle = len(experiments)

      experiments.extend(self.find_lattices())
      if len(experiments) == 0:
        raise Sorry("No suitable lattice could be found.")
      elif len(experiments) == n_lattices_previous_cycle:
        # no more lattices found
        break

      for i_cycle in range(self.params.refinement_protocol.n_macro_cycles):
        if i_cycle > 0:
          self.d_min -= self.params.refinement_protocol.d_min_step
          self.d_min = max(self.d_min, self.params.refinement_protocol.d_min_final)
          if self.d_min < 0:
            break
          print "Increasing resolution to %.1f Angstrom" %self.d_min

        # reset reflection lattice flags
        # the lattice a given reflection belongs to: a value of -1 indicates
        # that a reflection doesn't belong to any lattice so far
        self.reflections['id'] = flex.int(len(self.reflections), -1)

        if (i_cycle == 0 and self.target_symmetry_primitive is not None
            and self.target_symmetry_primitive.unit_cell() is not None):
          # if a target cell is given make sure that we match any permutation
          # of the cell dimensions
          for i_expt, expt in enumerate(experiments):
            if i_expt >= n_lattices_previous_cycle:
              expt.crystal, _ = self.apply_symmetry(
                expt.crystal, self.target_symmetry_primitive,
                return_primitive_setting=True,
                cell_only=True)

        self.index_reflections(experiments)

        if (i_cycle == 0 and self.target_symmetry_primitive is not None
            and self.target_symmetry_primitive.space_group() is not None):
          # now apply the space group symmetry only after the first indexing
          # need to make sure that the symmetrized orientation is similar to the P1 model
          for i_expt, expt in enumerate(experiments):
            if i_expt >= n_lattices_previous_cycle:
              expt.crystal, cb_op_to_primitive = self.apply_symmetry(
                expt.crystal, self.target_symmetry_primitive,
                return_primitive_setting=True,
                space_group_only=True)
              if not cb_op_to_primitive.is_identity_op():
                miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)
                miller_indices = cb_op_to_primitive.apply(miller_indices)
                self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)
              if self.cb_op_primitive_to_given is not None:
                expt.crystal = expt.crystal.change_basis(self.cb_op_primitive_to_given)
                miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)
                miller_indices = self.cb_op_primitive_to_given.apply(miller_indices)
                self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)

        if len(experiments) > 1:
          from dials.algorithms.indexing.compare_orientation_matrices \
               import difference_rotation_matrix_and_euler_angles
          cryst_b = experiments.crystals()[-1]
          have_similar_crystal_models = False
          for i_a, cryst_a in enumerate(experiments.crystals()[:-1]):
            R_ab, euler_angles, cb_op_ab = \
              difference_rotation_matrix_and_euler_angles(cryst_a, cryst_b)
            min_angle = self.params.multiple_lattice_search.minimum_angular_separation
            if max([abs(ea) for ea in euler_angles]) < min_angle: # degrees
              print "Crystal models too similar, rejecting crystal %i:" %(
                len(experiments))
              print "Rotation matrix to transform crystal %i to crystal %i" %(
                i_a+1, len(experiments))
              print R_ab
              print "Euler angles (xyz): %.2f, %.2f, %.2f" %euler_angles
              #show_rotation_matrix_differences([cryst_a, cryst_b])
              have_similar_crystal_models = True
              del experiments[-1]
              break
          if have_similar_crystal_models:
            break

        print
        print "#" * 80
        print "Starting refinement (macro-cycle %i)" %(i_cycle+1)
        print "#" * 80
        print
        self.indexed_reflections = (self.reflections['id'] > -1)

        if self.params.debug:
          sel = flex.bool(len(self.reflections), False)
          lengths = 1/self.reflections['rlp'].norms()
          isel = (lengths >= self.d_min).iselection()
          sel.set_selected(isel, True)
          sel.set_selected(self.reflections['id'] > -1, False)
          unindexed = self.reflections.select(sel)
          self.export_reflections(
            unindexed, "unindexed.pickle")

        maximum_spot_error \
          = self.params.refinement_protocol.outlier_rejection.maximum_spot_error
        maximum_phi_error \
          = self.params.refinement_protocol.outlier_rejection.maximum_phi_error
        if 1 and i_cycle == 0:
          maximum_spot_error = None
          maximum_phi_error = None
        elif i_cycle == 1:
          if maximum_spot_error is not None:
            maximum_spot_error *= 2
          if maximum_phi_error is not None:
            maximum_phi_error *= 2

        try:
          refined_experiments, refined_reflections = self.refine(
            experiments, maximum_spot_error=maximum_spot_error,
            maximum_phi_error=maximum_phi_error)
        except RuntimeError, e:
          s = str(e)
          if ("below the configured limit" in s or
              "Insufficient matches for crystal" in s):
            if len(experiments) == 1:
              raise Sorry(e)
            had_refinement_error = True
            print "Refinement failed:"
            print s
            del experiments[-1]
            break
          raise

        self.refined_reflections = refined_reflections

        self.sweep.set_detector(refined_experiments[0].detector)
        self.sweep.set_beam(refined_experiments[0].beam)
        self.sweep.set_goniometer(refined_experiments[0].goniometer)
        self.sweep.set_scan(refined_experiments[0].scan)
        # these may have been updated in refinement
        self.detector = self.sweep.get_detector()
        self.beam = self.sweep.get_beam()
        self.goniometer = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        if not (self.params.refinement.parameterisation.beam.fix == 'all'
                and self.params.refinement.parameterisation.detector.fix == 'all'):
          # Experimental geometry may have changed - re-map centroids to
          # reciprocal space
          self.map_centroids_to_reciprocal_space(
            self.reflections, self.detector, self.beam, self.goniometer)

        if self.d_min == self.params.refinement_protocol.d_min_final:
          print "Target d_min_final reached: finished with refinement"
          break

        # update for next cycle
        experiments = refined_experiments
        self.refined_experiments = refined_experiments

      if not self.params.multiple_lattice_search.recycle_unindexed_reflections:
        break

    # XXX currently need to store the imageset otherwise we can't read the experiment list back in
    for expt in refined_experiments:
      expt.imageset = self.sweep

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import show_rotation_matrix_differences
      show_rotation_matrix_differences(self.refined_experiments.crystals())

    print "Final refined crystal models:"
    for i, crystal_model in enumerate(self.refined_experiments.crystals()):
      print "model %i (%i reflections):" %(
        i+1, (self.reflections['id'] == i).count(True))
      print crystal_model

  def find_max_cell(self):
    import libtbx
    if self.params.max_cell is libtbx.Auto:
      if self.params.known_symmetry.unit_cell is not None:
        reference = self.target_symmetry_primitive.as_reference_setting()
        uc_params  = reference.unit_cell().parameters()
        self.params.max_cell = 1.5 * max(uc_params[:3])
      else:
        # The nearest neighbour analysis gets fooled when the same part of
        # reciprocal space has been measured twice as this introduced small
        # random differences in position between reflections measured twice.
        # Therefore repeat the nearest neighbour analysis several times in small
        # wedges where there shouldn't be any overlap in reciprocal space
        from rstbx.indexing_api.nearest_neighbor import neighbor_analysis
        phi_deg = self.reflections['xyzobs.mm.value'].parts()[2] * (180/math.pi)
        if (flex.max(phi_deg) - flex.min(phi_deg)) < 1e-3:
          NN = neighbor_analysis(self.reflections['rlp'])
          self.params.max_cell = NN.max_cell
        else:
          phi_min = flex.min(phi_deg)
          phi_max = flex.max(phi_deg)
          step_size = 5 #degrees
          d_phi = phi_max - phi_min
          n_steps = int(math.ceil(d_phi / step_size))
          max_cell = flex.double()
          for n in range(n_steps):
            sel = (phi_deg > (phi_min+n*step_size)) & (phi_deg < (phi_min+(n+1)*step_size))
            rlp = self.reflections['rlp'].select(sel)
            if len(rlp) == 0:
              continue
            NN = neighbor_analysis(rlp)
            max_cell.append(NN.max_cell)
            #print NN.max_cell
          self.params.max_cell = flex.mean(max_cell) # or max or median?
        print "Found max_cell: %.1f Angstrom" %(self.params.max_cell)

  def filter_reflections_by_scan_range(self):
    reflections_in_scan_range = flex.size_t()
    for i_ref, refl in enumerate(self.reflections):
      frame_number = refl['xyzobs.px.value'][2]

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

  @staticmethod
  def map_spots_pixel_to_mm_rad(spots, detector, scan):
    """Reflections that come from dials.spotfinder only have the centroid
    position and variance set, """
    from dials.algorithms.centroid import centroid_px_to_mm_panel
    ## ideally don't copy, but have separate spot attributes for mm and pixel
    import copy
    spots = copy.deepcopy(spots)
    spots['xyzobs.mm.value'] = flex.vec3_double(len(spots))
    spots['xyzobs.mm.variance'] = flex.vec3_double(len(spots))
    panel_numbers = flex.size_t(spots['panel'])
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      isel = sel.iselection()
      spots_panel = spots.select(panel_numbers == i_panel)
      centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
        detector[i_panel], scan,
        spots_panel['xyzobs.px.value'],
        spots_panel['xyzobs.px.variance'],
        flex.vec3_double(len(spots_panel), (1,1,1)))
      spots['xyzobs.mm.value'].set_selected(sel, centroid_position)
      spots['xyzobs.mm.variance'].set_selected(sel, centroid_variance)
    return spots

  @staticmethod
  def map_centroids_to_reciprocal_space(spots_mm, detector, beam, goniometer):
    if 's1' not in spots_mm: spots_mm['s1'] = flex.vec3_double(len(spots_mm))
    spots_mm['rlp'] = flex.vec3_double(len(spots_mm))
    panel_numbers = flex.size_t(spots_mm['panel'])
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      spots_panel = spots_mm.select(panel_numbers == i_panel)
      x, y, rot_angle = spots_panel['xyzobs.mm.value'].parts()
      s1 = detector[i_panel].get_lab_coord(flex.vec2_double(x,y))
      s1 = s1/s1.norms() * (1/beam.get_wavelength())
      S = s1 - beam.get_s0()
      # XXX what about if goniometer fixed rotation is not identity?
      if goniometer is not None:
        spots_mm['rlp'].set_selected(sel, S.rotate_around_origin(
          goniometer.get_rotation_axis(),
          -rot_angle))
      else:
        spots_mm['rlp'].set_selected(sel, S)

  def discover_better_experimental_model(self, reflections, detector, beam,
                                         goniometer, scan):
    from rstbx.phil.phil_preferences import indexing_api_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
      input_string=indexing_api_defs).extract()
    params = hardcoded_phil
    import copy

    from dials.algorithms.indexing import indexer
    opt_detector, opt_beam = indexer.discover_better_experimental_model(
      copy.deepcopy(reflections), detector, beam,
      goniometer=goniometer, scan=scan, params=hardcoded_phil)
    print "DISCOVERED BETTER MODEL:"
    print opt_detector
    print opt_beam
    return opt_detector, opt_beam

  def find_candidate_orientation_matrices(self, candidate_basis_vectors,
                                          max_combinations=1,
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
          model_orig = model
          cb_op_to_niggli = uc.change_of_basis_op_to_niggli_cell()
          model = model.change_basis(cb_op_to_niggli)
          uc = model.get_unit_cell()
          if (self.target_symmetry_centred is not None or
              self.target_symmetry_primitive is not None or
              self.target_symmetry_minimum_cell is not None):
            symmetrized_model = None
            if self.target_symmetry_minimum_cell is not None:
              symmetrized_model, cb_op_to_primitive = self.apply_symmetry(
              model, self.target_symmetry_minimum_cell,
              return_primitive_setting=True)
            if symmetrized_model is None and self.target_symmetry_primitive is not None:
              symmetrized_model, cb_op_to_primitive = self.apply_symmetry(
                model, self.target_symmetry_primitive,
              return_primitive_setting=True)
            if symmetrized_model is None and self.target_symmetry_centred is not None:
              symmetrized_model, cb_op_to_primitive = self.apply_symmetry(
              model, self.target_symmetry_centred,
              return_primitive_setting=True)
            if symmetrized_model is None:
              continue
            if apply_symmetry:
              model = symmetrized_model
              uc = model.get_unit_cell()
            elif cb_op_to_primitive is not None:
              model = model.change_basis(cb_op_to_primitive)
              uc = model.get_unit_cell()

          params = uc.parameters()
          if uc.volume() > (params[0]*params[1]*params[2]/100):
            # unit cell volume cutoff from labelit 2004 paper
            candidate_crystal_models.append(model)
            if len(candidate_crystal_models) == max_combinations:
              return candidate_crystal_models
    return candidate_crystal_models

  def choose_best_orientation_matrix(self, candidate_orientation_matrices):

    from libtbx import group_args
    class Solution(group_args):
      pass

    # Tracker for solutions based on code in rstbx/dps_core/basis_choice.py
    class SolutionTracker(object):
      def __init__(self):
        self.all_solutions = []
        self.volume_filtered = []

      def append(self, item):
        self.all_solutions.append(item)
        self.update_analysis()

      def __len__(self):
        return len(self.volume_filtered)

      def update_analysis(self):
        self.best_likelihood = max(
          s.model_likelihood for s in self.all_solutions)
        self.close_solutions = [
        s for s in self.all_solutions
        if s.model_likelihood >= (0.9 * self.best_likelihood)]
        self.min_volume = min(s.crystal.get_unit_cell().volume()
                              for s in self.close_solutions)
        self.volume_filtered = [
          s for s in self.close_solutions
          if s.crystal.get_unit_cell().volume() < (1.25 * self.min_volume)]
        self.best_volume_filtered_liklihood = max(
          s.model_likelihood for s in self.volume_filtered)

      def best_solution(self):
        solutions = [s for s in self.volume_filtered
                    if s.model_likelihood == self.best_volume_filtered_liklihood]
        return solutions[0]

    # Tracker for solutions based on code in rstbx/dps_core/basis_choice.py
    class SolutionTrackerSimple(object):
      def __init__(self):
        self.all_solutions = []

      def append(self, item):
        self.all_solutions.append(item)
        self.update_analysis()

      def __len__(self):
        return len(self.all_solutions)

      def update_analysis(self):
        self.best_n_indexed = max(
          s.n_indexed for s in self.all_solutions)

      def best_solution(self):
        solutions = [s for s in self.all_solutions
                    if s.n_indexed == self.best_n_indexed]
        return solutions[0]

    if self.params.basis_vector_combinations.metric == "model_likelihood":
      solutions = SolutionTracker()
    else:
      solutions = SolutionTrackerSimple()

    n_indexed = flex.int()
    min_likelihood = 0.3
    best_likelihood = min_likelihood

    import copy
    params = copy.deepcopy(self.params)
    params.refinement.parameterisation.crystal.scan_varying = False
    params.refinement.parameterisation.detector.fix = "all"
    params.refinement.parameterisation.beam.fix = "all"
    params.refinement.refinery.max_iterations = 1
    params.refinement.reflections.minimum_number_of_reflections = 1

    from dials.algorithms.indexing.compare_orientation_matrices \
         import difference_rotation_matrix_and_euler_angles
    for cm in candidate_orientation_matrices:
      if (self.refined_experiments is not None and
          len(self.refined_experiments) > 0):
        orientation_too_similar = False
        cryst_b = cm
        if self.target_symmetry_primitive is not None:
          cryst_b, cb_op_to_primitive = self.apply_symmetry(
            cryst_b, self.target_symmetry_primitive,
            return_primitive_setting=True,
            cell_only=True)
          cryst_b, cb_op_to_primitive = self.apply_symmetry(
            cryst_b, self.target_symmetry_primitive,
            return_primitive_setting=True,
            space_group_only=True)
          if self.cb_op_primitive_to_given is not None:
            cryst_b = cryst_b.change_basis(self.cb_op_primitive_to_given)
        for i_a, cryst_a in enumerate(self.refined_experiments.crystals()):
          R_ab, euler_angles, cb_op_ab = \
            difference_rotation_matrix_and_euler_angles(cryst_a, cryst_b)
          min_angle = self.params.multiple_lattice_search.minimum_angular_separation
          #print euler_angles
          if max([abs(ea) for ea in euler_angles]) < min_angle: # degrees
            orientation_too_similar = True
            break
        if orientation_too_similar:
          if params.debug:
            print "skipping crystal: too similar to other crystals"
          continue

      sel = ((self.reflections['id'] == -1) &
             (1/self.reflections['rlp'].norms() > self.d_min))
      refl = self.reflections.select(sel)
      experiment = Experiment(beam=self.beam,
                              detector=self.detector,
                              goniometer=self.goniometer,
                              scan=self.imagesets[0].get_scan(),
                              crystal=cm)
      if self.params.index_assignment.method == 'local':
        params_local = self.params.index_assignment.local
        from dials.algorithms.indexing import index_reflections_local
        index_reflections_local(
          refl, ExperimentList([experiment]),
          self.d_min, epsilon=params_local.epsilon,
          delta=params_local.delta, l_min=params_local.l_min,
          nearest_neighbours=params_local.nearest_neighbours)
      else:
        params_simple = self.params.index_assignment.simple
        from dials.algorithms.indexing import index_reflections
        index_reflections(refl, ExperimentList([experiment]), self.d_min,
                          tolerance=params_simple.hkl_tolerance)
      n_indexed.append((refl['id'] > -1).count(True))
      from dials.algorithms.refinement import RefinerFactory

      if (self.target_symmetry_primitive is not None
          and self.target_symmetry_primitive.space_group() is not None):
        experiment.crystal, cb_op_to_primitive = self.apply_symmetry(
          experiment.crystal, self.target_symmetry_primitive,
          return_primitive_setting=True,
          space_group_only=True)
        if not cb_op_to_primitive.is_identity_op():
          miller_indices = refl['miller_index'].select(refl['id'] == 0)
          miller_indices = cb_op_to_primitive.apply(miller_indices)
          refl['miller_index'].set_selected(refl['id'] == 0, miller_indices)
        if self.cb_op_primitive_to_given is not None:
          experiment.crystal = experiment.crystal.change_basis(self.cb_op_primitive_to_given)
          miller_indices = refl['miller_index'].select(refl['id'] == 0)
          miller_indices = self.cb_op_primitive_to_given.apply(miller_indices)
          refl['miller_index'].set_selected(refl['id'] == 0, miller_indices)

      refiner = RefinerFactory.from_parameters_data_experiments(
        params, refl.select(refl['id'] > -1), ExperimentList([experiment]),
        verbosity=0)
      refiner.run()
      rmsds = refiner.rmsds()
      xy_rmsds = math.sqrt(rmsds[0]**2 + rmsds[1]**2)
      model_likelihood = 1.0 - xy_rmsds
      if model_likelihood < 0:
        continue
      solutions.append(
        Solution(model_likelihood=model_likelihood,
                 crystal=cm,
                 rmsds=rmsds,
                 n_indexed=n_indexed[-1]))
      if self.params.debug:
        print "model_likelihood: %.2f" %model_likelihood
        print "n_indexed: %i" %n_indexed[-1]
      if model_likelihood > best_likelihood:
        best_likelihood = model_likelihood

    if len(solutions):
      best_solution = solutions.best_solution()
      if self.params.debug:
        print "best model_likelihood: %.2f" %best_solution.model_likelihood
        print "best n_indexed: %i" %best_solution.n_indexed
      return best_solution.crystal, best_solution.n_indexed
    else:
      return None, None

  def apply_symmetry(self, crystal_model, target_symmetry,
                     return_primitive_setting=False,
                     cell_only=False,
                     space_group_only=False):
    assert [cell_only, space_group_only].count(True) < 2
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
    if space_group_only:
      best_perm = (0,1,2)
      best_sign = 1
    else:
      best_perm = None
      # for non-cyclic permutations one axis needs inverting to keep system right-handed
      for perm, sign in zip(
          ((0,1,2), (1,2,0), (2,0,1), (1,0,2), (0,2,1), (2,1,0)),
          (1, 1, 1, -1, -1, -1)):
        a = matrix.col(basis_vectors[perm[0]])
        b = matrix.col(basis_vectors[perm[1]])
        c = matrix.col([i * sign for i in basis_vectors[perm[2]]])
        unit_cell = uctbx.unit_cell((a.length(),
                                     b.length(),
                                     c.length(),
                                     b.angle(c, deg=True),
                                     c.angle(a, deg=True),
                                     a.angle(b, deg=True)))
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
          bmsd = unit_cell.change_basis(
            symm_target_sg.change_of_basis_op_to_reference_setting())\
            .bases_mean_square_difference(
              symm_target_sg.as_reference_setting().unit_cell())
          eps = 1e-8
          if (bmsd+eps) < min_bmsd:
            min_bmsd = bmsd
            best_perm = list(perm)
            best_sign = sign
      if best_perm is None:
        return None, None
    direct_matrix = matrix.sqr(
      basis_vectors[best_perm[0]] + basis_vectors[best_perm[1]] +
      tuple([i * best_sign for i in basis_vectors[best_perm[2]]]))

    from cctbx.crystal_orientation import crystal_orientation
    from rstbx import dps_core # import dependency
    from rstbx.dps_core.lepage import iotbx_converter
    from dials.algorithms.indexing.symmetry import dials_crystal_from_orientation

    orientation = crystal_orientation(direct_matrix.inverse(), True)
    items = iotbx_converter(orientation.unit_cell(), max_delta=5.0)
    target_sg_reference_setting \
      = target_space_group.info().reference_setting().group()
    if target_space_group.crystal_system().lower() != 'triclinic':
      for item in items:
        if (target_space_group.crystal_system() !=
            item['best_group'].crystal_system()):
          continue
        if (target_sg_reference_setting.conventional_centring_type_symbol() !=
            item['best_group'].conventional_centring_type_symbol()):
          continue
        cb_op = item['cb_op_inp_best'].c().as_double_array()[0:9]
        cb_op_inv = item['cb_op_inp_best'].inverse().c().as_double_array()[0:9]
        orient = crystal_orientation(direct_matrix.inverse(), True)
        orient_best = orient.change_basis(matrix.sqr(cb_op).transpose())
        constrain_orient = orient_best.constrain(item['system'])
        orientation = constrain_orient.change_basis(
          matrix.sqr(cb_op_inv).transpose())
        target_space_group = target_sg_reference_setting.change_basis(
          item['cb_op_inp_best'].inverse())
        break

    direct_matrix = orientation.direct_matrix()
    a = matrix.col(direct_matrix[:3])
    b = matrix.col(direct_matrix[3:6])
    c = matrix.col(direct_matrix[6:9])
    ## verify it is still right-handed basis set
    #assert a.cross(b).dot(c) > 0
    cb_op_to_primitive = None
    model = Crystal(
      a, b, c, space_group=target_space_group)
    if return_primitive_setting:
      # for some reason I do not understand, it is necessary to go via the
      # reference setting to ensure we get the correct primitive setting
      cb_op_to_reference \
        = model.get_space_group().info().change_of_basis_op_to_reference_setting()
      cb_op_to_primitive = cb_op_to_reference
      cb_op_to_primitive.update(
        model.get_space_group().info().reference_setting().change_of_basis_op_to_primitive_setting())
      model = model.change_basis(cb_op_to_primitive)
    return model, cb_op_to_primitive

  def index_reflections(self, experiments):
    if self.params.index_assignment.method == 'local':
      params_local = self.params.index_assignment.local
      from dials.algorithms.indexing import index_reflections_local
      index_reflections_local(
        self.reflections,
        experiments, self.d_min, epsilon=params_local.epsilon,
        delta=params_local.delta, l_min=params_local.l_min,
        nearest_neighbours=params_local.nearest_neighbours,
        verbosity=self.params.refinement_protocol.verbosity)
    else:
      params_simple = self.params.index_assignment.simple
      from dials.algorithms.indexing import index_reflections
      index_reflections(self.reflections,
                        experiments, self.d_min,
                        tolerance=params_simple.hkl_tolerance,
                        verbosity=self.params.refinement_protocol.verbosity)

  def refine(self, experiments, maximum_spot_error=None,
             maximum_phi_error=None):
    from dials.algorithms.indexing.refinement import refine
    reflections_for_refinement = self.reflections.select(
      self.indexed_reflections)
    refiner, refined, outliers = refine(
      self.params, reflections_for_refinement, experiments,
      maximum_spot_error=maximum_spot_error,
      maximum_phi_error=maximum_phi_error,
      verbosity=self.params.refinement_protocol.verbosity,
      debug_plots=self.params.debug_plots)
    if outliers is not None:
      self.reflections['id'].set_selected(
        self.indexed_reflections.iselection().select(outliers), -1)
    used_reflections = refiner.get_reflections()
    verbosity = self.params.refinement_protocol.verbosity
    matches = refiner.get_matches()
    xyzcal_mm = flex.vec3_double(len(used_reflections))
    xyzcal_mm.set_selected(matches['iobs'], matches['xyzcal.mm'])
    used_reflections['xyzcal.mm'] = xyzcal_mm
    used_reflections.set_flags(
      matches['iobs'], used_reflections.flags.used_in_refinement)
    return refiner.get_experiments(), used_reflections

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
    for site in self.reflections['rlp']:
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

  def export_as_json(self, experiments, file_name="experiments.json",
                     compact=False):
    from dxtbx.serialize import dump
    assert experiments.is_consistent()
    dump.experiment_list(experiments, file_name)

  def export_reflections(self, reflections, file_name="reflections.pickle"):
    easy_pickle.dump(file_name, reflections)

  def find_lattices(self):
    raise NotImplementedError()




def optimise_basis_vectors(reciprocal_lattice_points, vectors):
  optimised = flex.vec3_double()
  for vector in vectors:
    minimised = basis_vector_minimser(reciprocal_lattice_points, vector)
    optimised.append(tuple(minimised.x))
  return optimised


from scitbx import lbfgs

# Optimise the initial basis vectors as per equation 11.4.3.4 of
# Otwinowski et al, International Tables Vol. F, chapter 11.4 pp. 282-295
class basis_vector_target(object):

  def __init__(self, reciprocal_lattice_points):
    self.reciprocal_lattice_points = reciprocal_lattice_points
    self._xyz_parts = self.reciprocal_lattice_points.parts()

  def compute_functional_and_gradients(self, vector):
    assert len(vector) == 3
    two_pi_S_dot_v = 2 * math.pi * self.reciprocal_lattice_points.dot(vector)
    f = - flex.sum(flex.cos(two_pi_S_dot_v))
    sin_part = flex.sin(two_pi_S_dot_v)
    g = flex.double([flex.sum(2 * math.pi * self._xyz_parts[i] * sin_part)
                     for i in range(3)])
    return f, g


class basis_vector_minimser(object):
  def __init__(self, reciprocal_lattice_points, vector,
               lbfgs_termination_params=None,
               lbfgs_core_params=lbfgs.core_parameters(m=20)):
    self.reciprocal_lattice_points = reciprocal_lattice_points
    if not isinstance(vector, flex.double):
      self.x = flex.double(vector)
    else:
      self.x = vector.deep_copy()
    self.n = len(self.x)
    assert self.n == 3
    self.target = basis_vector_target(self.reciprocal_lattice_points)
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

