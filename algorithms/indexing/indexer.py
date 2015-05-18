#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.indexer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math
import logging
from logging import info, debug
from dials.util import log

debug_handle = log.debug_handle()
info_handle = log.info_handle()

from libtbx.utils import Sorry
import iotbx.phil
from scitbx import matrix

from dials.array_family import flex
from cctbx import crystal, sgtbx, xray

from dxtbx.model.crystal import crystal_model as Crystal
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList

index_only_phil_str = """\
indexing {
  nproc = 1
    .type = int(value_min=1)
    .help = "The number of processes to use."
  discover_better_experimental_model = False
    .type = bool
    .expert_level = 1
  mm_search_scope = 4.0
    .help = "Global radius of origin offset search."
    .type = float(value_min=0)
    .expert_level = 1
  wide_search_binning = 2
    .help = "Modify the coarseness of the wide grid search for the beam centre."
    .type = float(value_min=0)
    .expert_level = 1
  min_cell = 20
    .type = float(value_min=0)
    .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
    .expert_level = 1
  max_cell = Auto
    .type = float(value_min=0)
    .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
  max_cell_multiplier = 1.3
    .type = float(value_min=0)
    .help = "Multiply the estimated maximum basis vector length by this value."
  nearest_neighbor_percentile = 0.05
    .type = float(value_min=0)
    .help = "Percentile of NN histogram to use for max cell determination."
    .expert_level = 1
  fft3d {
    peak_search = *flood_fill clean
      .type = choice
      .expert_level = 2
    reciprocal_space_grid {
      n_points = 256
        .type = int(value_min=0)
        .expert_level = 1
      d_min = Auto
        .type = float(value_min=0)
        .help = "The high resolution limit in Angstrom for spots to include in "
                "the initial indexing."
    }
  }
  sigma_phi_deg = None
    .type = float(value_min=0)
    .help = "Override the phi sigmas for refinement. Mainly intended for single-shot"
            "rotation images where the phi sigma is almost certainly incorrect."
    .expert_level = 2
  b_iso = 200
    .type = float(value_min=0)
    .expert_level = 2
  rmsd_cutoff = 15
    .type = float(value_min=0)
    .expert_level = 1
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
      .help = "Target space group for indexing."
    unit_cell = None
      .type = unit_cell
      .help = "Target unit cell for indexing."
    relative_length_tolerance = 0.1
      .type = float
      .help = "Relative tolerance for unit cell lengths in unit cell comparision."
      .expert_level = 1
    absolute_angle_tolerance = 5
      .type = float
      .help = "Angular tolerance (in degrees) in unit cell comparison."
      .expert_level = 1
    max_delta = 5
      .type = float(value_min=0)
      .help = "Maximum allowed Le Page delta used in searching for basis vector"
              "combinations that are consistent with the given symmetry."
      .expert_level = 1
  }
  basis_vector_combinations {
    max_try = 50
      .type = int(value_min=1)
      .help = "Number of putative basis vector combinations to try."
      .expert_level = 1
    filter
      .expert_level = 1
    {
      check_doubled_cell = True
        .type = bool
      likelihood_cutoff = 0.8
        .type = float(value_min=0, value_max=1)
      volume_cutoff = 1.25
        .type = float(value_min=1)
      n_indexed_cutoff = 0.9
        .type = float(value_min=0, value_max=1)
    }
  }
  index_assignment {
    method = *simple local
      .type = choice
      .help = "Choose between simple 'global' index assignment and xds-style "
              "'local' index assignment."
      .expert_level = 1
    simple {
      hkl_tolerance = 0.3
        .type = float(value_min=0, value_max=0.5)
        .help = "Maximum allowable deviation from integer-ness for assigning "
                "a miller index to a reciprocal lattice vector."
    }
    local
      .expert_level = 1
    {
      epsilon = 0.05
        .type = float
        .help = "This corresponds to the xds parameter INDEX_ERROR="
      delta = 8
        .type = int
        .help = "This corresponds to the xds parameter INDEX_MAGNITUDE="
      l_min = 0.8
        .type = float
        .help = "This corresponds to the xds parameter INDEX_QUALITY="
      nearest_neighbours = 20
        .type = int(value_min=1)
    }
  }
  optimise_initial_basis_vectors = False
    .type = bool
    .expert_level = 2
  debug = False
    .type = bool
    .expert_level = 1
  debug_plots = False
    .type = bool
    .help = "Requires matplotlib"
    .expert_level = 1
  refinement_protocol {
    n_macro_cycles = 5
      .type = int(value_min=1)
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
    disable_unit_cell_volume_sanity_check = False
      .type = bool
      .help = "Disable sanity check on unrealistic increases in unit cell volume"
              "during refinement."
      .expert_level = 1
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
                "differ by more than the given value (degrees)."
                "No outlier rejection is performed in the first macro cycle, and "
                "in the second macro cycle twice the given multiple is used."
    }
  }
  method = *fft3d fft1d real_space_grid_search
    .type = choice
  multiple_lattice_search
    .expert_level = 1
  {
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
  real_space_grid_search
    .expert_level = 1
  {
    characteristic_grid = 0.02
      .type = float(value_min=0)
  }
}
"""

index_only_phil_scope = iotbx.phil.parse(index_only_phil_str, process_includes=True)

master_phil_scope = iotbx.phil.parse("""
%s
include scope dials.algorithms.refinement.refiner.phil_scope
""" %index_only_phil_str, process_includes=True)

master_params = master_phil_scope.fetch().extract()


def filter_reflections_by_scan_range(reflections, scan_range):
  reflections_in_scan_range = flex.bool(len(reflections), False)
  frame_number = reflections['xyzobs.px.value'].parts()[2]

  for scan_range in scan_range:
    if scan_range is None: continue
    range_start, range_end = scan_range
    reflections_in_scan_range.set_selected(
      (frame_number >= range_start) & (frame_number < range_end), True)
  return reflections.select(reflections_in_scan_range)


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
    self.reflections = reflections
    self.imagesets = imagesets
    for imageset in imagesets[1:]:
      assert imageset.get_detector() == self.imagesets[0].get_detector()
      assert imageset.get_beam() == self.imagesets[0].get_beam()
      assert imageset.get_goniometer() == self.imagesets[0].get_goniometer()
    sweep = self.imagesets[0]
    self.sweep = sweep

    if params is None: params = master_params

    self.goniometer = sweep.get_goniometer()
    self.detector = sweep.get_detector()
    #self.scan = sweep.get_scan()
    self.beam = sweep.get_beam()
    self.params = params.indexing
    self.all_params = params
    self.refined_experiments = None

    if 'flags' in self.reflections:
      strong_sel = self.reflections.get_flags(self.reflections.flags.strong)
      if strong_sel.count(True) > 0:
        self.reflections = self.reflections.select(strong_sel)
    if 'flags' not in self.reflections or strong_sel.count(True) == 0:
      # backwards compatibility for testing
      self.reflections.set_flags(
        flex.size_t_range(len(self.reflections)), self.reflections.flags.strong)

    self._setup_symmetry()

    # now actually do the indexing
    self.index()


  def _setup_symmetry(self):
    self.target_symmetry_primitive = None
    self.target_symmetry_reference_setting = None
    self.cb_op_inp_ref = None

    target_unit_cell = self.params.known_symmetry.unit_cell
    target_space_group = self.params.known_symmetry.space_group
    if target_space_group is not None:
      target_space_group = target_space_group.group()

    if target_unit_cell is not None or target_space_group is not None:

      if target_unit_cell is not None and target_space_group is not None:
        from cctbx.sgtbx.bravais_types import bravais_lattice
        target_bravais_t = bravais_lattice(
          group=target_space_group.info().reference_setting().group())
        best_subgroup = None
        best_angular_difference = 1e8
        from cctbx.sgtbx import lattice_symmetry
        space_groups = [target_space_group]
        if target_space_group.conventional_centring_type_symbol() != 'P':
          space_groups.append(sgtbx.space_group())
        for target_space_group in space_groups:
          cs = crystal.symmetry(
            unit_cell=target_unit_cell, space_group=target_space_group,
            assert_is_compatible_unit_cell=False)
          target_best_cell = cs.best_cell().unit_cell()
          subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=0.1)
          for subgroup in subgroups.result_groups:
            bravais_t = bravais_lattice(
              group=subgroup['ref_subsym'].space_group())
            if bravais_t == target_bravais_t:
              # allow for the cell to be given as best cell, reference setting
              # primitive settings, or minimum cell
              best_subsym = subgroup['best_subsym']
              ref_subsym = best_subsym.as_reference_setting()
              if not (
                best_subsym.unit_cell().is_similar_to(target_unit_cell) or
                ref_subsym.unit_cell().is_similar_to(target_unit_cell) or
                ref_subsym.primitive_setting().unit_cell().is_similar_to(target_unit_cell) or
                best_subsym.primitive_setting().unit_cell().is_similar_to(target_unit_cell) or
                best_subsym.minimum_cell().unit_cell().is_similar_to(
                  target_unit_cell.minimum_cell()) or
                best_subsym.unit_cell().is_similar_to(target_best_cell)):
                continue
              if subgroup['max_angular_difference'] < best_angular_difference:
                best_subgroup = subgroup
                best_angular_difference = subgroup['max_angular_difference']

        if best_subgroup is None:
          raise Sorry("Unit cell incompatible with space group")

        cb_op_inp_best = best_subgroup['cb_op_inp_best']
        best_subsym = best_subgroup['best_subsym']
        cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
        self.cb_op_inp_ref = cb_op_best_ref * cb_op_inp_best
        self.target_symmetry_reference_setting = crystal.symmetry(
          unit_cell=target_unit_cell.change_basis(self.cb_op_inp_ref),
          space_group=self.params.known_symmetry.space_group.as_reference_setting().group())

      elif target_unit_cell is not None:
        self.target_symmetry_reference_setting = crystal.symmetry(
          unit_cell=target_unit_cell,
          space_group=sgtbx.space_group())
        self.cb_op_inp_ref = sgtbx.change_of_basis_op()

      elif target_space_group is not None:
        self.cb_op_inp_ref = target_space_group.info().change_of_basis_op_to_reference_setting()
        self.target_symmetry_reference_setting = crystal.symmetry(
          space_group=target_space_group.change_basis(self.cb_op_inp_ref))

      self.cb_op_reference_to_primitive \
        = self.target_symmetry_reference_setting.change_of_basis_op_to_primitive_setting()
      if target_unit_cell is not None:
        self.target_symmetry_primitive \
          = self.target_symmetry_reference_setting.change_basis(
            self.cb_op_reference_to_primitive)
      else:
        self.target_symmetry_primitive = crystal.symmetry(
          space_group=self.target_symmetry_reference_setting.space_group()\
          .change_basis(self.cb_op_reference_to_primitive))
      self.cb_op_ref_inp = self.cb_op_inp_ref.inverse()
      self.cb_op_primitive_inp \
        = self.cb_op_ref_inp * self.cb_op_reference_to_primitive.inverse()

      if self.target_symmetry_reference_setting is not None:
        debug("Target symmetry (reference setting):")
        self.target_symmetry_reference_setting.show_summary(f=debug_handle)
      if self.target_symmetry_primitive is not None:
        debug("Target symmetry (primitive cell):")
        self.target_symmetry_primitive.show_summary(f=debug_handle)
      debug("cb_op reference->primitive: " + str(self.cb_op_reference_to_primitive))
      debug("cb_op primitive->input: " + str(self.cb_op_primitive_inp))

  def index(self):
    self.reflections_input = self.reflections
    self.reflections = flex.reflection_table()
    for i, imageset in enumerate(self.imagesets):
      if 'imageset_id' in self.reflections_input:
        sel = (self.reflections_input['imageset_id'] == i)
      else:
        sel = (self.reflections_input['id'] == i)
      self.reflections.extend(self.map_spots_pixel_to_mm_rad(
        self.reflections_input.select(sel),
        imageset.get_detector(), imageset.get_scan()))
    self.filter_reflections_by_scan_range()
    if len(self.reflections) == 0:
      raise Sorry("No reflections left to index!")

    if self.params.discover_better_experimental_model:

      from dials.command_line.discover_better_experimental_model \
           import discover_better_experimental_model

      from rstbx.phil.phil_preferences import indexing_api_defs
      import iotbx.phil
      hardcoded_phil = iotbx.phil.parse(
        input_string=indexing_api_defs).extract()
      hardcoded_phil.indexing.mm_search_scope = self.params.mm_search_scope

      new_detector, new_beam = discover_better_experimental_model(
        [self.imagesets[0]], [self.reflections], hardcoded_phil, nproc=self.params.nproc,
        wide_search_binning=self.params.wide_search_binning)

      self.sweep.set_detector(new_detector)
      self.sweep.set_beam(new_beam)
      self.detector = new_detector
      self.beam = new_beam

    self.map_centroids_to_reciprocal_space(
      self.reflections, self.detector, self.beam, self.goniometer)

    try:
      self.find_max_cell()
    except AssertionError, e:
      if "too few spots" in str(e).lower():
        raise Sorry(e)

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

    while True:
      self.d_min = self.params.refinement_protocol.d_min_start
      if had_refinement_error or have_similar_crystal_models:
        break
      max_lattices = self.params.multiple_lattice_search.max_lattices
      if max_lattices is not None and len(experiments) >= max_lattices:
        break
      if len(experiments) > 0:
        cutoff_fraction = \
          self.params.multiple_lattice_search.recycle_unindexed_reflections_cutoff
        d_spacings = 1/self.reflections['rlp'].norms()
        d_min_indexed = flex.min(d_spacings.select(self.indexed_reflections))
        min_reflections_for_indexing = \
          cutoff_fraction * len(self.reflections.select(d_spacings > d_min_indexed))
        crystal_ids = self.reflections.select(d_spacings > d_min_indexed)['id']
        if (crystal_ids == -1).count(True) < min_reflections_for_indexing:
          info("Finish searching for more lattices: %i unindexed reflections remaining." %(
            min_reflections_for_indexing))
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
          d_min = self.d_min - self.params.refinement_protocol.d_min_step
          d_min = max(d_min, 0)
          d_min = max(d_min, self.params.refinement_protocol.d_min_final)
          if d_min >= 0:
            self.d_min = d_min
            info("Increasing resolution to %.1f Angstrom" %d_min)

        # reset reflection lattice flags
        # the lattice a given reflection belongs to: a value of -1 indicates
        # that a reflection doesn't belong to any lattice so far
        self.reflections['id'] = flex.int(len(self.reflections), -1)

        #if (i_cycle == 0 and self.target_symmetry_primitive is not None
            #and self.target_symmetry_primitive.unit_cell() is not None):
          ## if a target cell is given make sure that we match any permutation
          ## of the cell dimensions
          #for i_cryst, cryst in enumerate(experiments.crystals()):
            #if i_cryst >= n_lattices_previous_cycle:
              #new_cryst, _ = self.apply_symmetry(
                #cryst, self.target_symmetry_primitive,
                #return_primitive_setting=True,
                #cell_only=True)
              #experiments.crystals()[i_cryst].update(new_cryst)

        self.index_reflections(
          experiments, self.reflections,
          verbosity=self.params.refinement_protocol.verbosity)

        if (i_cycle == 0 and self.target_symmetry_primitive is not None
            and self.target_symmetry_primitive.space_group() is not None):
          # now apply the space group symmetry only after the first indexing
          # need to make sure that the symmetrized orientation is similar to the P1 model
          for i_expt, expt in enumerate(experiments):
            if i_expt >= n_lattices_previous_cycle:
              new_cryst, cb_op_to_primitive = self.apply_symmetry(
                expt.crystal, self.target_symmetry_primitive,
                space_group_only=True)
              if not cb_op_to_primitive.is_identity_op():
                miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)
                miller_indices = cb_op_to_primitive.apply(miller_indices)
                self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)
              if self.cb_op_primitive_inp is not None:
                new_cryst = new_cryst.change_basis(self.cb_op_primitive_inp)
                info(new_cryst.get_space_group().info())
                miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)
                miller_indices = self.cb_op_primitive_inp.apply(miller_indices)
                self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)
              expt.crystal.update(new_cryst)

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
              info("Crystal models too similar, rejecting crystal %i:" %(
                len(experiments)))
              info("Rotation matrix to transform crystal %i to crystal %i" %(
                i_a+1, len(experiments)))
              info(R_ab)
              info("Euler angles (xyz): %.2f, %.2f, %.2f" %euler_angles)
              #show_rotation_matrix_differences([cryst_a, cryst_b])
              have_similar_crystal_models = True
              del experiments[-1]
              break
          if have_similar_crystal_models:
            break

        info("")
        info("#" * 80)
        info("Starting refinement (macro-cycle %i)" %(i_cycle+1))
        info("#" * 80)
        info("")
        self.indexed_reflections = (self.reflections['id'] > -1)

        sel = flex.bool(len(self.reflections), False)
        lengths = 1/self.reflections['rlp'].norms()
        isel = (lengths >= self.d_min).iselection()
        sel.set_selected(isel, True)
        sel.set_selected(self.reflections['id'] > -1, False)
        self.unindexed_reflections = self.reflections.select(sel)

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

        reflections_for_refinement = self.reflections.select(
          self.indexed_reflections)
        try:
          refined_experiments, refined_reflections = self.refine(
            experiments, reflections_for_refinement,
            maximum_spot_error=maximum_spot_error,
            maximum_phi_error=maximum_phi_error)
        except RuntimeError, e:
          s = str(e)
          if ("below the configured limit" in s or
              "Insufficient matches for crystal" in s):
            if len(experiments) == 1:
              raise Sorry(e)
            had_refinement_error = True
            info("Refinement failed:")
            info(s)
            del experiments[-1]
            break
          raise

        # sanity check for unrealistic unit cell volume increase during refinement
        # usually this indicates too many parameters are being refined given the
        # number of observations provided.
        if not self.params.refinement_protocol.disable_unit_cell_volume_sanity_check:
          for orig_expt, refined_expt in zip(experiments, refined_experiments):
            uc1 = orig_expt.crystal.get_unit_cell()
            uc2 = refined_expt.crystal.get_unit_cell()
            volume_change = abs(uc1.volume()-uc2.volume())/uc1.volume()
            cutoff = 0.5
            if volume_change > cutoff:
              msg = "\n".join((
                "Unrealistic unit cell volume increase during refinement of %.1f%%.",
                "Please try refining fewer parameters, either by enforcing symmetry",
                "constraints (space_group=) and/or disabling experimental geometry",
                "refinement (detector.fix=all and beam.fix=all). To disable this",
                "sanity check set disable_unit_cell_volume_sanity_check=True.")) %(
                100*volume_change)
              raise Sorry(msg)

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

        if not (self.all_params.refinement.parameterisation.beam.fix == 'all'
                and self.all_params.refinement.parameterisation.detector.fix == 'all'):
          # Experimental geometry may have changed - re-map centroids to
          # reciprocal space
          self.map_centroids_to_reciprocal_space(
            self.reflections, self.detector, self.beam, self.goniometer)

        # update for next cycle
        experiments = refined_experiments
        self.refined_experiments = refined_experiments

        if (i_cycle >=2 and
            self.d_min == self.params.refinement_protocol.d_min_final or
            self.d_min <= 0):
          info("Target d_min_final reached: finished with refinement")
          break

      if not self.params.multiple_lattice_search.recycle_unindexed_reflections:
        break

    if not 'refined_experiments' in locals():
      raise Sorry("None of the experiments could refine.")

    # XXX currently need to store the imageset otherwise we can't read the experiment list back in
    for expt in refined_experiments:
      expt.imageset = self.sweep

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import show_rotation_matrix_differences
      show_rotation_matrix_differences(
        self.refined_experiments.crystals(), out=info_handle)

    info("Final refined crystal models:")
    for i, crystal_model in enumerate(self.refined_experiments.crystals()):
      info("model %i (%i reflections):" %(
        i+1, (self.reflections['id'] == i).count(True)))
      info(crystal_model)

    # set xyzcal.px field in self.refined_reflections
    panel_numbers = flex.size_t(self.refined_reflections['panel'])
    xyzcal_mm = self.refined_reflections['xyzcal.mm']
    x_mm, y_mm, z_rad = xyzcal_mm.parts()
    xy_cal_mm = flex.vec2_double(x_mm, y_mm)
    xy_cal_px = flex.vec2_double(len(xy_cal_mm))
    for i_panel in range(len(self.detector)):
      panel = self.detector[i_panel]
      sel = (panel_numbers == i_panel)
      isel = sel.iselection()
      ref_panel = self.refined_reflections.select(panel_numbers == i_panel)
      xy_cal_px.set_selected(
        sel, panel.millimeter_to_pixel(xy_cal_mm.select(sel)))
    x_px, y_px = xy_cal_px.parts()
    scan = self.sweep.get_scan()
    if scan is not None:
      z_px = scan.get_array_index_from_angle(z_rad, deg=False)
    else:
      # must be a still image, z centroid not meaningful
      z_px = z_rad
    xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
    self.refined_reflections['xyzcal.px'] = xyzcal_px

  def find_max_cell(self):
    import libtbx
    if self.params.max_cell is libtbx.Auto:
      if self.params.known_symmetry.unit_cell is not None:
        uc_params = self.target_symmetry_primitive.unit_cell().parameters()
        self.params.max_cell = self.params.max_cell_multiplier * max(uc_params[:3])
        info("Using max_cell: %.1f Angstrom" %(self.params.max_cell))
      else:
        # The nearest neighbour analysis gets fooled when the same part of
        # reciprocal space has been measured twice as this introduced small
        # random differences in position between reflections measured twice.
        # Therefore repeat the nearest neighbour analysis several times in small
        # wedges where there shouldn't be any overlap in reciprocal space
        #from rstbx.indexing_api.nearest_neighbor import neighbor_analysis
        from dials.algorithms.indexing.nearest_neighbor import neighbor_analysis
        phi_deg = self.reflections['xyzobs.mm.value'].parts()[2] * (180/math.pi)
        if (flex.max(phi_deg) - flex.min(phi_deg)) < 1e-3:
          NN = neighbor_analysis(self.reflections['rlp'],
                                 tolerance=self.params.max_cell_multiplier,
                                 percentile=self.params.nearest_neighbor_percentile)
          self.params.max_cell = NN.max_cell
        else:
          phi_min = flex.min(phi_deg)
          phi_max = flex.max(phi_deg)
          step_size = 5 #degrees
          d_phi = phi_max - phi_min
          n_steps = int(math.ceil(d_phi / step_size))
          max_cell = flex.double()
          for n in range(n_steps):
            sel = (phi_deg >= (phi_min+n*step_size)) & (phi_deg < (phi_min+(n+1)*step_size))
            rlp = self.reflections['rlp'].select(sel)
            if len(rlp) == 0:
              continue
            try:
              NN = neighbor_analysis(
                rlp, tolerance=self.params.max_cell_multiplier,
                percentile=self.params.nearest_neighbor_percentile)
              max_cell.append(NN.max_cell)
            except AssertionError:
              continue
            debug("%s %s %s"  %(
              phi_min+n*step_size, phi_min+(n+1)*step_size, NN.max_cell))
          debug(list(max_cell))
          debug("median: %s" %flex.median(max_cell))
          debug("mean: %s" %flex.mean(max_cell))
          self.params.max_cell = flex.median(max_cell) # mean or max or median?
        info("Found max_cell: %.1f Angstrom" %(self.params.max_cell))

  def filter_reflections_by_scan_range(self):
    if len(self.params.scan_range):
      self.reflections = filter_reflections_by_scan_range(
        self.reflections, self.params.scan_range)

  @staticmethod
  def map_spots_pixel_to_mm_rad(spots, detector, scan):
    """Map spot centroids from pixel/image number to mm/radian.

    Used to convert spot centroids coming from e.g. dials.find_spots which are in
    pixel/image number units to mm/radian units as required for indexing and refinement.

    :param spots: a reflection table containing the columns 'xyzobs.px.value',
                  'xyzobs.px.variance' and 'panel'
    :type spots: dials.array_family.flex.reflection_table
    :param detector: a dxtbx detector object
    :type detector: dxtbx.model.detector.Detector
    :param scan: a dxtbx scan object. May be None, e.g. for a still image
    :type scan: dxtbx.model.scan.Scan
    :returns: A copy of the input reflection table containing the additional keys
              'xyzobs.mm.value' and 'xyzobs.mm.variance'
    :rtype: dials.array_family.flex.reflection_table
    """

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
    """Map mm/radian spot centroids to reciprocal space.

    Used to convert spot centroids provided in mm/radian units to reciprocal space
    as required for indexing. Adds the column 'rlp' to the reflection table, which
    contains a :py:class:`.flex.vec3_double` array of the reciprocal lattice vectors.

    :param spots_mm: a reflection table containing the column 'xyzobs.mm.value'
    :type spots_mm: dials.array_family.flex.reflection_table
    :param detector: a dxtbx detector object
    :type detector: dxtbx.model.detector.Detector
    :param beam: a dxtbx beam object
    :type beam: dxtbx.model.beam.Beam
    :param goniometer: a dxtbx goniometer object. May be None, e.g. for a still image
    :type goniometer: dxtbx.model.goniometer.Goniometer
    """

    if 'imageset_id' not in spots_mm:
      spots_mm['imageset_id'] = spots_mm['id']
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

  def find_candidate_orientation_matrices(self, candidate_basis_vectors,
                                          max_combinations=1):
    candidate_crystal_models = []
    vectors = candidate_basis_vectors

    # select unique combinations of input vectors to test
    # the order of combinations is such that combinations comprising vectors
    # nearer the beginning of the input list will appear before combinations
    # comprising vectors towards the end of the list
    n = len(vectors)
    combinations = flex.vec3_int(flex.nested_loop((n,n,n)))
    combinations = combinations.select(
      flex.sort_permutation(combinations.as_vec3_double().norms()))

    # select only those combinations where j > i and k > j
    i, j, k = combinations.as_vec3_double().parts()
    sel = flex.bool(len(combinations), True)
    sel &= j > i
    sel &= k > j
    combinations = combinations.select(sel)

    min_angle = 20 # degrees, arbitrary cutoff
    for i, j, k in combinations:
      a = vectors[i]
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
      c = vectors[k]
      if abs(90-a_cross_b.angle(c, deg=True)) < min_angle:
        continue
      alpha = b.angle(c, deg=True)
      if alpha < 90:
        c = -c
      #beta = c.angle(a, deg=True)
      if a_cross_b.dot(c) < 0:
        # we want right-handed basis set, therefore invert all vectors
        a = -a
        b = -b
        c = -c
        #assert a.cross(b).dot(c) > 0
      model = Crystal(a, b, c, space_group_symbol="P 1")
      uc = model.get_unit_cell()
      model_orig = model
      best_model = None
      min_bmsd = 1e8
      cb_op_to_niggli = uc.change_of_basis_op_to_niggli_cell()
      model = model.change_basis(cb_op_to_niggli)
      #print model.get_unit_cell()
      uc = model.get_unit_cell()
      if self.target_symmetry_primitive is not None:
        max_delta = self.params.known_symmetry.max_delta
        from dials.algorithms.indexing.symmetry import find_matching_symmetry
        best_subgroup = find_matching_symmetry(
          uc, self.target_symmetry_primitive.space_group(),
          max_delta=max_delta)
        cb_op_extra = None
        if best_subgroup is None:
          if self.target_symmetry_reference_setting is not None:
            # if we have been told we have a centred unit cell check that
            # indexing hasn't found the centred unit cell instead of the
            # primitive cell
            best_subgroup = find_matching_symmetry(
              uc, self.target_symmetry_reference_setting.space_group().build_derived_point_group(),
              max_delta=max_delta)
            cb_op_extra = self.cb_op_reference_to_primitive
            if best_subgroup is None:
              continue
          else:
            continue
        cb_op_inp_best = best_subgroup['cb_op_inp_best']
        best_subsym = best_subgroup['best_subsym']
        cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
        ref_subsym = best_subsym.change_basis(cb_op_best_ref)
        cb_op_ref_primitive = ref_subsym.change_of_basis_op_to_primitive_setting()
        cb_op_to_primitive = cb_op_ref_primitive * cb_op_best_ref * cb_op_inp_best
        if cb_op_extra is not None:
          cb_op_to_primitive = cb_op_extra * cb_op_to_primitive
        best_model = model.change_basis(cb_op_to_primitive)
        if (self.target_symmetry_primitive.unit_cell () is not None
            and not best_model.get_unit_cell().is_similar_to(
          self.target_symmetry_primitive.unit_cell(),
          relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
            absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance)):
          continue
      else:
        best_model = model

      if best_model is None:
        continue

      uc = best_model.get_unit_cell()
      params = uc.parameters()
      if uc.volume() > (params[0]*params[1]*params[2]/100):
        # unit cell volume cutoff from labelit 2004 paper
        candidate_crystal_models.append(best_model)
        if len(candidate_crystal_models) == max_combinations:
          return candidate_crystal_models
    return candidate_crystal_models

  def choose_best_orientation_matrix(self, candidate_orientation_matrices):

    filter_params = self.params.basis_vector_combinations.filter
    solutions = SolutionTracker(
      check_doubled_cell=filter_params.check_doubled_cell,
      likelihood_cutoff=filter_params.likelihood_cutoff,
      volume_cutoff=filter_params.volume_cutoff,
      n_indexed_cutoff=filter_params.n_indexed_cutoff)

    def run_one_refinement(args):
      params, reflections, experiments = args
      from dials.algorithms.refinement import RefinerFactory
      try:
        logger = logging.getLogger()
        disabled = logger.disabled
        logger.disabled = True
        refiner = RefinerFactory.from_parameters_data_experiments(
          params, reflections, experiments,
          verbosity=0)
        refiner.run()
      except RuntimeError, e:
        #print e
        return
      else:
        rmsds = refiner.rmsds()
        xy_rmsds = math.sqrt(rmsds[0]**2 + rmsds[1]**2)
        model_likelihood = 1.0 - xy_rmsds
        soln = Solution(model_likelihood=model_likelihood,
                        crystal=experiments.crystals()[0],
                        rmsds=rmsds,
                        n_indexed=len(reflections))
        return soln
      finally:
        logger.disabled = disabled

    import copy
    params = copy.deepcopy(self.all_params)
    params.refinement.parameterisation.crystal.scan_varying = False
    #params.refinement.parameterisation.detector.fix = "all"
    #params.refinement.parameterisation.beam.fix = "all"
    params.refinement.refinery.max_iterations = 2
    params.refinement.reflections.minimum_number_of_reflections = 1

    args = []

    from dials.algorithms.indexing.compare_orientation_matrices \
         import difference_rotation_matrix_and_euler_angles
    for cm in candidate_orientation_matrices:
      sel = ((self.reflections['id'] == -1) &
             (1/self.reflections['rlp'].norms() > self.d_min))
      refl = self.reflections.select(sel)
      experiment = Experiment(beam=self.beam,
                              detector=self.detector,
                              goniometer=self.goniometer,
                              scan=self.imagesets[0].get_scan(),
                              crystal=cm)
      self.index_reflections(ExperimentList([experiment]), refl)

      if (self.target_symmetry_primitive is not None
          and self.target_symmetry_primitive.space_group() is not None):
        new_crystal, cb_op_to_primitive = self.apply_symmetry(
          experiment.crystal, self.target_symmetry_primitive,
          space_group_only=True)
        if new_crystal is None:
          continue
        experiment.crystal = new_crystal
        if not cb_op_to_primitive.is_identity_op():
          miller_indices = refl['miller_index'].select(refl['id'] == 0)
          miller_indices = cb_op_to_primitive.apply(miller_indices)
          refl['miller_index'].set_selected(refl['id'] == 0, miller_indices)
        if 0 and self.cb_op_primitive_to_given is not None:
          experiment.crystal = experiment.crystal.change_basis(self.cb_op_primitive_to_given)
          miller_indices = refl['miller_index'].select(refl['id'] == 0)
          miller_indices = self.cb_op_primitive_to_given.apply(miller_indices)
          refl['miller_index'].set_selected(refl['id'] == 0, miller_indices)

      if (self.refined_experiments is not None and
          len(self.refined_experiments) > 0):
        orientation_too_similar = False
        cryst_b = experiment.crystal
        for i_a, cryst_a in enumerate(self.refined_experiments.crystals()):
          R_ab, euler_angles, cb_op_ab = \
            difference_rotation_matrix_and_euler_angles(cryst_a, cryst_b)
          min_angle = self.params.multiple_lattice_search.minimum_angular_separation
          #print euler_angles
          if max([abs(ea) for ea in euler_angles]) < min_angle: # degrees
            orientation_too_similar = True
            break
        if orientation_too_similar:
          debug("skipping crystal: too similar to other crystals")
          continue

      args.append(
        (params, refl.select(refl['id'] > -1), ExperimentList([experiment])))

    from libtbx import easy_mp
    results = easy_mp.parallel_map(
      run_one_refinement,
      args,
      processes=self.params.nproc,
      preserve_exception_message=True,
    )

    for soln in results:
      if soln is None:
        continue
      solutions.append(soln)
      debug("unit_cell: " + str(soln.crystal.get_unit_cell()))
      debug("model_likelihood: %.2f" %soln.model_likelihood)
      debug("n_indexed: %i" %soln.n_indexed)

    if len(solutions):
      best_solution = solutions.best_solution()
      debug("best model_likelihood: %.2f" %best_solution.model_likelihood)
      debug("best n_indexed: %i" %best_solution.n_indexed)
      return best_solution.crystal, best_solution.n_indexed
    else:
      return None, None

  def apply_symmetry(self, crystal_model, target_symmetry,
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

    from cctbx.crystal_orientation import crystal_orientation
    from cctbx.sgtbx.bravais_types import bravais_lattice
    from rstbx import dps_core # import dependency
    from rstbx.dps_core.lepage import iotbx_converter

    max_delta = self.params.known_symmetry.max_delta
    items = iotbx_converter(crystal_model.get_unit_cell(), max_delta=max_delta)
    target_sg_ref = target_space_group.info().reference_setting().group()
    best_angular_difference = 1e8
    best_subgroup = None
    for item in items:
      if (bravais_lattice(group=target_sg_ref) !=
          bravais_lattice(group=item['ref_subsym'].space_group())):
        continue
      if item['max_angular_difference'] < best_angular_difference:
        best_angular_difference = item['max_angular_difference']
        best_subgroup = item

    if best_subgroup is None:
      return None, None

    cb_op_inp_best = best_subgroup['cb_op_inp_best']
    orient = crystal_orientation(A, True)
    orient_best = orient.change_basis(
      matrix.sqr(cb_op_inp_best.c().as_double_array()[0:9]).transpose())
    constrain_orient = orient_best.constrain(best_subgroup['system'])
    #best_orientation = constrain_orient

    best_subsym = best_subgroup['best_subsym']
    cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
    target_sg_best = target_sg_ref.change_basis(cb_op_best_ref.inverse())
    ref_subsym = best_subsym.change_basis(cb_op_best_ref)
    cb_op_ref_primitive = ref_subsym.change_of_basis_op_to_primitive_setting()
    primitive_subsym = ref_subsym.change_basis(cb_op_ref_primitive)
    cb_op_best_primitive = cb_op_ref_primitive * cb_op_best_ref
    cb_op_inp_primitive = cb_op_ref_primitive * cb_op_best_ref * cb_op_inp_best

    direct_matrix = constrain_orient.direct_matrix()

    a = matrix.col(direct_matrix[:3])
    b = matrix.col(direct_matrix[3:6])
    c = matrix.col(direct_matrix[6:9])
    model = Crystal(
      a, b, c, space_group=target_sg_best)

    model = model.change_basis(cb_op_best_primitive)
    return model, cb_op_inp_primitive

  def index_reflections(self, experiments, reflections, verbosity=0):
    if self.params.index_assignment.method == 'local':
      params_local = self.params.index_assignment.local
      from dials.algorithms.indexing import index_reflections_local
      index_reflections_local(
        reflections,
        experiments, self.d_min, epsilon=params_local.epsilon,
        delta=params_local.delta, l_min=params_local.l_min,
        nearest_neighbours=params_local.nearest_neighbours,
        verbosity=verbosity)
    else:
      params_simple = self.params.index_assignment.simple
      from dials.algorithms.indexing import index_reflections
      index_reflections(reflections,
                        experiments, self.d_min,
                        tolerance=params_simple.hkl_tolerance,
                        verbosity=verbosity)

  def refine(self, experiments, reflections, maximum_spot_error=None,
             maximum_phi_error=None):
    from dials.algorithms.indexing.refinement import refine
    refiner, refined, outliers = refine(
      self.all_params, reflections, experiments,
      maximum_spot_error=maximum_spot_error,
      maximum_phi_error=maximum_phi_error,
      verbosity=self.params.refinement_protocol.verbosity,
      debug_plots=self.params.debug_plots)
    if outliers is not None:
      reflections['id'].set_selected(outliers, -1)
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
      debug("%s %s" %(i, v.length()))# , vector_heights[i]

    # print a table of the angles between each pair of vectors

    angles = flex.double(len(vectors)**2)
    angles.reshape(flex.grid(len(vectors), len(vectors)))

    for i in range(len(vectors)):
      v_i = vectors[i]
      for j in range(i+1, len(vectors)):
        v_j = vectors[j]
        angles[i,j] = v_i.angle(v_j, deg=True)

    print >> debug_handle, (" "*7),
    for i in range(len(vectors)):
      print >> debug_handle, "%7.3f" % vectors[i].length(),
    print >> debug_handle
    for i in range(len(vectors)):
      print >> debug_handle, "%7.3f" % vectors[i].length(),
      for j in range(len(vectors)):
        if j <= i:
          print >> debug_handle, (" "*7),
        else:
          print >> debug_handle, "%5.1f  " %angles[i,j],
      print >> debug_handle

  def debug_plot_candidate_basis_vectors(self):
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D # import dependency
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([0], [0], [0], marker='+', s=50)

    #http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
    #draw a vector
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
            FancyArrowPatch.draw(self, renderer)
    for v in self.candidate_basis_vectors:
      x, y, z = v.elems
      a = Arrow3D(
        [0,x],[0,y],[0,z], mutation_scale=10, lw=1, arrowstyle="-|>", color="k")
      ax.add_artist(a)
      a = Arrow3D(
        [0,-x],[0,-y],[0,-z], mutation_scale=10, lw=1, arrowstyle="-|>", color="k")
      ax.add_artist(a)

    x, y, z = zip(*self.candidate_basis_vectors)
    ax.scatter(x, y, z, marker='.', s=1)
    ax.scatter([-i for i in x], [-i for i in y], [-i for i in z], marker='.', s=1)
    pyplot.show()

  def debug_write_reciprocal_lattice_points_as_pdb(
      self, file_name='reciprocal_lattice.pdb'):
    from cctbx import crystal, xray
    cs = crystal.symmetry(unit_cell=(1000,1000,1000,90,90,90), space_group="P1")
    for i_panel in range(len(self.detector)):
      if len(self.detector) > 1:
        file_name = 'reciprocal_lattice_%i.pdb' %i_panel
      with open(file_name, 'wb') as f:
        xs = xray.structure(crystal_symmetry=cs)
        reflections = self.reflections.select(self.reflections['panel'] == i_panel)
        for site in reflections['rlp']:
          xs.add_scatterer(xray.scatterer("C", site=site))
        xs.sites_mod_short()
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
    from libtbx import easy_pickle
    easy_pickle.dump(file_name, reflections)

  def find_lattices(self):
    raise NotImplementedError()


from libtbx import group_args
class Solution(group_args):
  pass

# Tracker for solutions based on code in rstbx/dps_core/basis_choice.py
class SolutionTracker(object):
  def __init__(self, check_doubled_cell=True, likelihood_cutoff=0.8,
               volume_cutoff=1.25, n_indexed_cutoff=0.9):
    self.check_doubled_cell = check_doubled_cell
    self.likelihood_cutoff = likelihood_cutoff
    self.volume_cutoff = volume_cutoff
    self.n_indexed_cutoff = n_indexed_cutoff
    self.all_solutions = []
    self.filtered_solutions = []

  def append(self, item):
    self.all_solutions.append(item)
    self.update_analysis()

  def __len__(self):
    return len(self.filtered_solutions)

  def filter_doubled_cell(self, solutions):
    from dials.algorithms.indexing.compare_orientation_matrices import difference_rotation_matrix_and_euler_angles
    accepted_solutions = []
    for i1, s1 in enumerate(solutions):
      doubled_cell = False
      for (m1,m2,m3) in ((2,1,1), (1,2,1), (1,1,2), (2,2,1), (2,1,2), (1,2,2), (2,2,2)):
        if doubled_cell:
          break
        a, b, c = s1.crystal.get_real_space_vectors()
        new_cryst = Crystal(real_space_a=1/m1*a,
                            real_space_b=1/m2*b,
                            real_space_c=1/m3*c,
                            space_group=s1.crystal.get_space_group())
        new_unit_cell = new_cryst.get_unit_cell()
        for s2 in solutions:
          if s2 is s1:
            continue
          if new_unit_cell.is_similar_to(s2.crystal.get_unit_cell(),
                                         relative_length_tolerance=0.05):
            R, ea, cb = difference_rotation_matrix_and_euler_angles(new_cryst, s2.crystal)
            if ((max(ea) < 1) and
                (s1.model_likelihood < (1.5 * s2.model_likelihood))):
              doubled_cell = True
              break

      if not doubled_cell:
        accepted_solutions.append(s1)

    return accepted_solutions

  def filter_by_likelihood(self, solutions):
    best_likelihood = max(s.model_likelihood for s in solutions)
    offset = 0
    while (best_likelihood + offset) <= 0:
      offset += 1
    return [
      s for s in solutions
      if (s.model_likelihood+offset) >= (self.likelihood_cutoff * (best_likelihood+offset))]

  def filter_by_volume(self, solutions):
    # filter by volume - prefer solutions with a smaller unit cell
    min_volume = min(s.crystal.get_unit_cell().volume()
                          for s in solutions)
    return [
      s for s in solutions
      if s.crystal.get_unit_cell().volume() < (self.volume_cutoff * min_volume)]

  def filter_by_n_indexed(self, solutions, n_indexed_cutoff=None):
    if n_indexed_cutoff is None:
      n_indexed_cutoff = self.n_indexed_cutoff
    # filter by number of indexed reflections - prefer solutions that
    # account for more of the diffracted spots
    max_n_indexed = max(s.n_indexed for s in solutions)
    return [s for s in solutions
            if s.n_indexed >= n_indexed_cutoff * max_n_indexed]

  def update_analysis(self):
    # pre-filter out solutions that only account for a very small
    # percentage of the indexed spots relative to the best one
    self.filtered_solutions = self.filter_by_n_indexed(
      self.all_solutions, n_indexed_cutoff=0.05) # 5 percent

    self.filtered_solutions = self.filter_doubled_cell(self.filtered_solutions)

    self.filtered_solutions = self.filter_by_likelihood(self.filtered_solutions)

    self.filtered_solutions = self.filter_by_volume(self.filtered_solutions)

    self.filtered_solutions = self.filter_by_n_indexed(self.filtered_solutions)

    return

  def best_solution(self):
    self.best_filtered_liklihood = max(
      s.model_likelihood for s in self.filtered_solutions)

    solutions = [s for s in self.filtered_solutions
                 if s.model_likelihood == self.best_filtered_liklihood]
    return solutions[0]


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

