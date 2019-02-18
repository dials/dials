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

from __future__ import absolute_import, division
from __future__ import print_function
import math
import logging
logger = logging.getLogger(__name__)

from dials.util import log

debug_handle = log.debug_handle(logger)

import libtbx
from dials.util import Sorry
import iotbx.phil
from scitbx import matrix

from dials.array_family import flex
from cctbx import crystal, sgtbx, xray

from dxtbx.model import Crystal
from dxtbx.model.experiment_list import Experiment, ExperimentList

max_cell_phil_str = '''\
max_cell_estimation
  .expert_level = 1
{
  filter_ice = True
    .type = bool
    .help = "Filter out reflections at typical ice ring resolutions"
            "before max_cell estimation."
  filter_overlaps = True
    .type = bool
    .help = "Filter out reflections with overlapping bounding boxes before"
            "max_cell estimation."
  overlaps_border = 0
    .type = int(value_min=0)
    .help = "Optionally add a border around the bounding boxes before finding"
            "overlaps."
  multiplier = 1.3
    .type = float(value_min=0)
    .help = "Multiply the estimated maximum basis vector length by this value."
    .expert_level = 2
  step_size = 45
    .type = float(value_min=0)
    .help = "Step size, in degrees, of the blocks used to peform the max_cell "
            "estimation."
    .expert_level = 2
  nearest_neighbor_percentile = None
    .type = float(value_min=0, value_max=1)
    .help = "Percentile of NN histogram to use for max cell determination."
    .expert_level = 2
  histogram_binning = linear *log
    .type = choice
    .help = "Choose between linear or logarithmic bins for nearest neighbour"
            "histogram analysis."
  nn_per_bin = 5
    .type = int(value_min=1)
    .help = "Target number of nearest neighbours per histogram bin."
  max_height_fraction = 0.25
    .type = float(value_min=0, value_max=1)
    .expert_level=2
}
'''

index_only_phil_str = '''\
indexing {
  nproc = 1
    .type = int(value_min=1)
    .help = "The number of processes to use."
  mm_search_scope = 4.0
    .help = "Global radius of origin offset search."
    .type = float(value_min=0)
    .expert_level = 1
  wide_search_binning = 2
    .help = "Modify the coarseness of the wide grid search for the beam centre."
    .type = float(value_min=0)
    .expert_level = 1
  min_cell_volume = 25
    .type = float(value_min=0)
    .help = "Minimum unit cell volume (in Angstrom^3)."
    .expert_level = 1
  min_cell = 3
    .type = float(value_min=0)
    .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
    .expert_level = 1
  max_cell = Auto
    .type = float(value_min=0)
    .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
    .expert_level = 1
  %s
  fft3d {
    peak_search = *flood_fill clean
      .type = choice
      .expert_level = 2
    peak_volume_cutoff = 0.15
      .type = float
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
  b_iso = Auto
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
  basis_vector_combinations
    .expert_level = 1
  {
    max_combinations = None
      .type = int(value_min=1)
      .help = "Maximum number of basis vector combinations to test for agreement"
              "with input symmetry."
    max_refine = Auto
      .type = int(value_min=1)
      .help = "Maximum number of putative crystal models to test. Default"
              "for rotation sweeps: 50, for still images: 5"
      .expert_level = 1
    sys_absent_threshold = 0.9
      .type = float(value_min=0.0, value_max=1.0)
    solution_scorer = filter *weighted
      .type = choice
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
    weighted
      .expert_level = 1
    {
      power = 1
        .type = int(value_min=1)
      volume_weight = 1
        .type = float(value_min=0)
      n_indexed_weight = 1
        .type = float(value_min=0)
      rmsd_weight = 1
        .type = float(value_min=0)
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
  check_misindexing {
    grid_search_scope = 0
      .type = int
      .help = "Search scope for testing misindexing on h, k, l."
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
  combine_scans = False
    .type = bool
    .expert_level = 1
  refinement_protocol {
    mode = *refine_shells repredict_only
      .type = choice
      .expert_level = 1
      .help = "refine_shells: refine in increasing resolution cutoffs after indexing."
              "repredict_only: do not refine after indexing, just update spot"
              "predictions."
    n_macro_cycles = 5
      .type = int(value_min=1)
      .help = "Maximum number of macro cycles of refinement, reindexing all"
              "reflections using updated geometry at the beginning of each"
              "cycle. Does not apply to stills.indexer=stills."
    d_min_step = Auto
      .type = float(value_min=0.0)
      .help = "Reduction per step in d_min for reflections to include in refinement."
    d_min_start = None
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
  }
  method = *fft3d fft1d real_space_grid_search
    .type = choice
  multiple_lattice_search
    .expert_level = 1
  {
    cluster_analysis_search = False
      .type = bool
      .help = "Perform cluster analysis search for multiple lattices."
    recycle_unindexed_reflections_cutoff = 0.1
      .type = float(value_min=0, value_max=1)
      .help = "Attempt another cycle of indexing on the unindexed reflections "
              "if more than the fraction of input reflections are unindexed."
    minimum_angular_separation = 5
      .type = float(value_min=0)
      .help = "The minimum angular separation (in degrees) between two lattices."
    max_lattices = 1
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
  fft1d
    .expert_level = 1
  {
    characteristic_grid = None
      .help = Sampling frequency in radians. See Steller 1997. If None, \
              determine a grid sampling automatically using the input \
              reflections, using at most 0.029 radians.
      .type = float(value_min=0)
  }
  real_space_grid_search
    .expert_level = 1
  {
    characteristic_grid = 0.02
      .type = float(value_min=0)
  }
  stills {
    indexer = *Auto stills sweeps
      .type = choice
      .help = Use the stills or sweeps indexer.  Auto: choose based on the input \
              imagesets (stills or sweeps).
      .expert_level = 1
    ewald_proximity_resolution_cutoff = 2.0
      .type = float
      .help = For still images, this high-resolution cutoff is used to calculate
      .help = the acceptable volume of reciprocal space for spot prediction
    refine_all_candidates = True
      .type = bool
      .help = If False, no attempt is made to refine the model from initial basis \
              vector selection. The indexing solution with the best RMSD is chosen.
    candidate_outlier_rejection = True
      .type = bool
      .expert_level = 1
      .help = If True, while refining candiate basis solutions, also apply Sauter/ \
              Poon (2010) outlier rejection
    refine_candidates_with_known_symmetry = False
      .type = bool
      .expert_level = 2
      .help = If False, when choosing the best set of candidate basis solutions, \
              refine the candidates in the P1 setting. If True, after indexing \
              in P1, convert the candidates to the known symmetry and apply the \
              corresponding change of basis to the indexed reflections.
    rmsd_min_px = 2
      .type = float
      .help = Minimum acceptable RMSD for choosing candidate basis solutions \
              (in pixels)
    ewald_proximal_volume_max = 0.0025
      .type = float
      .help = Maximum acceptable ewald proximal volume when choosing candidate \
              basis solutions
    isoforms
      .help = Constrain the unit cell to specific values during refinement after initial indexing.
      .multiple=True
    {
      name=None
        .type=str
      cell=None
        .type=unit_cell
      lookup_symbol=None
        .type=str
        .help=The sgtbx lookup symbol of the reflections pointgroup
      rmsd_target_mm=None
        .type=float
        .help=Maximum acceptable DIALS positional rmsd, in mm
      beam_restraint=None
        .type=floats(size=2)
        .help=Known beam position in mm X,Y, rmsd_target_mm is reused here as a circle of confusion
        .help=to assure that no images are accepted where the lattice is misindexed by a unit shift.
    }
  }
}
''' %max_cell_phil_str

index_only_phil_scope = iotbx.phil.parse(index_only_phil_str, process_includes=True)

master_phil_scope = iotbx.phil.parse("""
%s
include scope dials.algorithms.refinement.refiner.phil_scope
""" %index_only_phil_str, process_includes=True)

# override default refinement parameters
master_phil_scope = master_phil_scope.fetch(source=iotbx.phil.parse(
  """\
refinement {
  reflections {
    reflections_per_degree=100
  }
}
"""))

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
    self.volumes = []
    self._mean = None

  def append(self, vector, length, volume):
    self.vectors.append(vector)
    self.lengths.append(length)
    self.volumes.append(volume)
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

    if params is None: params = master_params

    self.params = params.indexing
    self.all_params = params
    self.refined_experiments = None
    self.hkl_offset = None

    if self.all_params.refinement.reflections.outlier.algorithm in ('auto', libtbx.Auto):
      if self.imagesets[0].get_goniometer() is None:
        self.all_params.refinement.reflections.outlier.algorithm = 'sauter_poon'
      else:
        # different default to dials.refine
        # tukey is faster and more appropriate at the indexing step
        self.all_params.refinement.reflections.outlier.algorithm = 'tukey'

    for imageset in imagesets[1:]:
      if imageset.get_detector().is_similar_to(self.imagesets[0].get_detector()):
        imageset.set_detector(self.imagesets[0].get_detector())
      if imageset.get_goniometer().is_similar_to(self.imagesets[0].get_goniometer()):
        imageset.set_goniometer(self.imagesets[0].get_goniometer())
        # can only share a beam if we share a goniometer?
        if imageset.get_beam().is_similar_to(self.imagesets[0].get_beam()):
          imageset.set_beam(self.imagesets[0].get_beam())
        if (self.params.combine_scans and
            imageset.get_scan() == self.imagesets[0].get_scan()):
          imageset.set_scan(self.imagesets[0].get_scan())

    if 'flags' in self.reflections:
      strong_sel = self.reflections.get_flags(self.reflections.flags.strong)
      if strong_sel.count(True) > 0:
        self.reflections = self.reflections.select(strong_sel)
    if 'flags' not in self.reflections or strong_sel.count(True) == 0:
      # backwards compatibility for testing
      self.reflections.set_flags(
        flex.size_t_range(len(self.reflections)), self.reflections.flags.strong)

    self._setup_symmetry()
    self.d_min = None

    self.setup_indexing()

  @staticmethod
  def from_parameters(reflections, imagesets,
                      known_crystal_models=None, params=None):

    if params is None:
      params = master_params

    if known_crystal_models is not None:
      from dials.algorithms.indexing.known_orientation \
           import indexer_known_orientation
      if params.indexing.known_symmetry.space_group is None:
        params.indexing.known_symmetry.space_group \
          = known_crystal_models[0].get_space_group().info()
      idxr = indexer_known_orientation(
        reflections, imagesets, params, known_crystal_models)
    else:
      has_stills = False
      has_sweeps = False
      for imageset in imagesets:
        if imageset.get_goniometer() is None or imageset.get_scan() is None or \
            imageset.get_scan().get_oscillation()[1] == 0:
          if has_sweeps:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_stills = True
        else:
          if has_stills:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_sweeps = True
      assert not (has_stills and has_sweeps)
      use_stills_indexer = has_stills

      if not (params.indexing.stills.indexer is libtbx.Auto or params.indexing.stills.indexer.lower() == 'auto'):
        if params.indexing.stills.indexer == 'stills':
          use_stills_indexer = True
        elif params.indexing.stills.indexer == 'sweeps':
          use_stills_indexer = False
        else:
          assert False

      if params.indexing.basis_vector_combinations.max_refine is libtbx.Auto:
        if use_stills_indexer:
          params.indexing.basis_vector_combinations.max_refine = 5
        else:
          params.indexing.basis_vector_combinations.max_refine = 50

      if use_stills_indexer:
        # Ensure the indexer and downstream applications treat this as set of stills
        from dxtbx.imageset import ImageSet#, MemImageSet
        reset_sets = []
        for i in xrange(len(imagesets)):
          imagesweep = imagesets.pop(0)
          imageset = ImageSet(imagesweep.data(), imagesweep.indices())
          # if isinstance(imageset, MemImageSet):
          #   imageset = MemImageSet(imagesweep._images, imagesweep.indices())
          # else:
          #   imageset = ImageSet(imagesweep.reader(), imagesweep.indices())
          #   imageset._models = imagesweep._models
          imageset.set_scan(None)
          imageset.set_goniometer(None)
          reset_sets.append(imageset)
        imagesets.extend(reset_sets)

      if params.indexing.method == "fft3d":
        if use_stills_indexer:
          from dials.algorithms.indexing.stills_indexer \
            import stills_indexer_fft3d as indexer_fft3d
        else:
          from dials.algorithms.indexing.fft3d import indexer_fft3d
        idxr = indexer_fft3d(reflections, imagesets, params=params)
      elif params.indexing.method == "fft1d":
        if use_stills_indexer:
          from dials.algorithms.indexing.stills_indexer \
            import stills_indexer_fft1d as indexer_fft1d
        else:
          from dials.algorithms.indexing.fft1d import indexer_fft1d
        idxr = indexer_fft1d(reflections, imagesets, params=params)
      elif params.indexing.method == "real_space_grid_search":
        if use_stills_indexer:
          from dials.algorithms.indexing.stills_indexer import \
            stills_indexer_real_space_grid_search as indexer_real_space_grid_search
        else:
          from dials.algorithms.indexing.real_space_grid_search \
            import indexer_real_space_grid_search
        idxr = indexer_real_space_grid_search(reflections, imagesets, params=params)

    return idxr

  def _setup_symmetry(self):
    self.target_symmetry_primitive = None
    self.target_symmetry_reference_setting = None
    self.cb_op_inp_ref = None

    target_unit_cell = self.params.known_symmetry.unit_cell
    target_space_group = self.params.known_symmetry.space_group
    if target_space_group is not None:
      target_space_group = target_space_group.group()
      target_space_group = target_space_group.build_derived_patterson_group()

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
        for target in space_groups:
          cs = crystal.symmetry(
            unit_cell=target_unit_cell, space_group=target,
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
          space_group=target_space_group.info().as_reference_setting().group())

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
        logger.debug("Target symmetry (reference setting):")
        self.target_symmetry_reference_setting.show_summary(f=debug_handle)
      if self.target_symmetry_primitive is not None:
        logger.debug("Target symmetry (primitive cell):")
        self.target_symmetry_primitive.show_summary(f=debug_handle)
      logger.debug("cb_op reference->primitive: " + str(self.cb_op_reference_to_primitive))
      logger.debug("cb_op primitive->input: " + str(self.cb_op_primitive_inp))

  def setup_indexing(self):
    reflections_input = self.reflections
    self.reflections = flex.reflection_table()
    for i, imageset in enumerate(self.imagesets):
      if 'imageset_id' not in reflections_input:
        reflections_input['imageset_id'] = reflections_input['id']
      sel = (reflections_input['imageset_id'] == i)
      self.reflections.extend(self.map_spots_pixel_to_mm_rad(
        reflections_input.select(sel),
        imageset.get_detector(), imageset.get_scan()))
    self.filter_reflections_by_scan_range()
    if len(self.reflections) == 0:
      raise Sorry("No reflections left to index!")

    spots_mm = self.reflections
    self.reflections = flex.reflection_table()

    for i, imageset in enumerate(self.imagesets):
      spots_sel = spots_mm.select(spots_mm['imageset_id'] == i)
      self.map_centroids_to_reciprocal_space(
        spots_sel, imageset.get_detector(), imageset.get_beam(),
        imageset.get_goniometer())
      spots_sel['entering'] = self.calculate_entering_flags(
        spots_sel, beam=imageset.get_beam(),
        goniometer=imageset.get_goniometer())
      self.reflections.extend(spots_sel)

    try:
      self.find_max_cell()
    except AssertionError as e:
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

  def index(self):

    experiments = ExperimentList()

    had_refinement_error = False
    have_similar_crystal_models = False

    while True:
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
          logger.info("Finish searching for more lattices: %i unindexed reflections remaining." %(
            min_reflections_for_indexing))
          break

      n_lattices_previous_cycle = len(experiments)

      if self.d_min is None:
        self.d_min = self.params.refinement_protocol.d_min_start

      if len(experiments) == 0:
        experiments.extend(self.find_lattices())
      else:
        try:
          new = self.find_lattices()
          experiments.extend(new)
        except Sorry:
          logger.info("Indexing remaining reflections failed")

      if self.params.refinement_protocol.d_min_step is libtbx.Auto:
        n_cycles = self.params.refinement_protocol.n_macro_cycles
        if self.d_min is None or n_cycles == 1:
          self.params.refinement_protocol.d_min_step = 0
        else:
          d_spacings = 1/self.reflections['rlp'].norms()
          d_min_all = flex.min(d_spacings)
          self.params.refinement_protocol.d_min_step \
            = (self.d_min - d_min_all)/(n_cycles-1)
          logger.info("Using d_min_step %.1f"
               %self.params.refinement_protocol.d_min_step)

      if len(experiments) == 0:
        raise Sorry("No suitable lattice could be found.")
      elif len(experiments) == n_lattices_previous_cycle:
        # no more lattices found
        break

      for i_cycle in range(self.params.refinement_protocol.n_macro_cycles):
        if (i_cycle > 0 and self.d_min is not None and
            self.params.refinement_protocol.d_min_step > 0):
          d_min = self.d_min - self.params.refinement_protocol.d_min_step
          d_min = max(d_min, 0)
          d_min = max(d_min, self.params.refinement_protocol.d_min_final)
          if d_min >= 0:
            self.d_min = d_min
            logger.info("Increasing resolution to %.2f Angstrom" %d_min)

        # reset reflection lattice flags
        # the lattice a given reflection belongs to: a value of -1 indicates
        # that a reflection doesn't belong to any lattice so far
        self.reflections['id'] = flex.int(len(self.reflections), -1)

        self.index_reflections(experiments, self.reflections)

        if i_cycle == 0 and self.params.known_symmetry.space_group is not None:
          # now apply the space group symmetry only after the first indexing
          # need to make sure that the symmetrized orientation is similar to the P1 model
          target_space_group = self.target_symmetry_primitive.space_group()
          for i_cryst, cryst in enumerate(experiments.crystals()):
            if i_cryst >= n_lattices_previous_cycle:
              new_cryst, cb_op_to_primitive = self.apply_symmetry(
                cryst, target_space_group)
              if self.cb_op_primitive_inp is not None:
                new_cryst = new_cryst.change_basis(self.cb_op_primitive_inp)
              cryst.update(new_cryst)
              cryst.set_space_group(
                self.params.known_symmetry.space_group.group())
              for i_expt, expt in enumerate(experiments):
                if expt.crystal is not cryst:
                  continue
                if not cb_op_to_primitive.is_identity_op():
                  miller_indices = self.reflections['miller_index'].select(
                    self.reflections['id'] == i_expt)
                  miller_indices = cb_op_to_primitive.apply(miller_indices)
                  self.reflections['miller_index'].set_selected(
                    self.reflections['id'] == i_expt, miller_indices)
                if self.cb_op_primitive_inp is not None:
                  miller_indices = self.reflections['miller_index'].select(
                    self.reflections['id'] == i_expt)
                  miller_indices = self.cb_op_primitive_inp.apply(miller_indices)
                  self.reflections['miller_index'].set_selected(
                    self.reflections['id'] == i_expt, miller_indices)
          logger.info("\nIndexed crystal models:")
          self.show_experiments(experiments, self.reflections, d_min=self.d_min)

        if len(experiments) > 1:
          from dials.algorithms.indexing.compare_orientation_matrices \
               import difference_rotation_matrix_axis_angle
          cryst_b = experiments.crystals()[-1]
          have_similar_crystal_models = False
          for i_a, cryst_a in enumerate(experiments.crystals()[:-1]):
            R_ab, axis, angle, cb_op_ab = \
              difference_rotation_matrix_axis_angle(cryst_a, cryst_b)
            min_angle = self.params.multiple_lattice_search.minimum_angular_separation
            if abs(angle) < min_angle: # degrees
              logger.info("Crystal models too similar, rejecting crystal %i:" %(
                len(experiments)))
              logger.info("Rotation matrix to transform crystal %i to crystal %i" %(
                i_a+1, len(experiments)))
              logger.info(R_ab)
              logger.info("Rotation of %.3f degrees" %angle + " about axis (%.3f, %.3f, %.3f)" %axis)
              #show_rotation_matrix_differences([cryst_a, cryst_b])
              have_similar_crystal_models = True
              del experiments[-1]
              break
          if have_similar_crystal_models:
            break

        logger.info("")
        logger.info("#" * 80)
        logger.info("Starting refinement (macro-cycle %i)" %(i_cycle+1))
        logger.info("#" * 80)
        logger.info("")
        self.indexed_reflections = (self.reflections['id'] > -1)

        sel = flex.bool(len(self.reflections), False)
        lengths = 1/self.reflections['rlp'].norms()
        if self.d_min is not None:
          isel = (lengths <= self.d_min).iselection()
          sel.set_selected(isel, True)
        sel.set_selected(self.reflections['id'] == -1, True)
        self.reflections.unset_flags(sel, self.reflections.flags.indexed)
        self.unindexed_reflections = self.reflections.select(sel)

        reflections_for_refinement = self.reflections.select(
          self.indexed_reflections)
        if self.params.refinement_protocol.mode == 'repredict_only':
          refined_experiments, refined_reflections = experiments, reflections_for_refinement
          from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
          ref_predictor = ExperimentsPredictorFactory.from_experiments(
              experiments, spherical_relp=self.all_params.refinement.parameterisation.spherical_relp_model)
          ref_predictor(refined_reflections)
        else:
          try:
            refined_experiments, refined_reflections = self.refine(
              experiments, reflections_for_refinement)
          except Sorry as e:
            if len(experiments) == 1:
              raise
            had_refinement_error = True
            logger.info("Refinement failed:")
            logger.info(e)
            del experiments[-1]
            break

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
        self.refined_reflections.unset_flags(
          self.refined_reflections['id'] < 0,
          self.refined_reflections.flags.indexed)

        for i, imageset in enumerate(self.imagesets):
          ref_sel = self.refined_reflections.select(
            self.refined_reflections['imageset_id'] == i)
          ref_sel = ref_sel.select(ref_sel['id'] >= 0)
          for i_expt in set(ref_sel['id']):
            expt = refined_experiments[i_expt]
            imageset.set_detector(expt.detector)
            imageset.set_beam(expt.beam)
            imageset.set_goniometer(expt.goniometer)
            imageset.set_scan(expt.scan)
            expt.imageset = imageset

        if not (self.all_params.refinement.parameterisation.beam.fix == 'all'
                and self.all_params.refinement.parameterisation.detector.fix == 'all'):
          # Experimental geometry may have changed - re-map centroids to
          # reciprocal space

          spots_mm = self.reflections
          self.reflections = flex.reflection_table()
          for i, imageset in enumerate(self.imagesets):
            spots_sel = spots_mm.select(spots_mm['imageset_id'] == i)
            self.map_centroids_to_reciprocal_space(
              spots_sel, imageset.get_detector(), imageset.get_beam(),
              imageset.get_goniometer())
            self.reflections.extend(spots_sel)

        # update for next cycle
        experiments = refined_experiments
        self.refined_experiments = refined_experiments

        logger.info("\nRefined crystal models:")
        self.show_experiments(
          self.refined_experiments, self.reflections, d_min=self.d_min)

        if (i_cycle >=2 and
            self.d_min == self.params.refinement_protocol.d_min_final):
          logger.info("Target d_min_final reached: finished with refinement")
          break

    if not 'refined_experiments' in locals():
      raise Sorry("None of the experiments could refine.")

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import rotation_matrix_differences
      logger.info(
        rotation_matrix_differences(self.refined_experiments.crystals()))

    self.refined_reflections['xyzcal.px'] = flex.vec3_double(
      len(self.refined_reflections))
    for i, imageset in enumerate(self.imagesets):
      imgset_sel = self.refined_reflections['imageset_id'] == i
      # set xyzcal.px field in self.refined_reflections
      refined_reflections = self.refined_reflections.select(imgset_sel)
      panel_numbers = flex.size_t(refined_reflections['panel'])
      xyzcal_mm = refined_reflections['xyzcal.mm']
      x_mm, y_mm, z_rad = xyzcal_mm.parts()
      xy_cal_mm = flex.vec2_double(x_mm, y_mm)
      xy_cal_px = flex.vec2_double(len(xy_cal_mm))
      for i_panel in range(len(imageset.get_detector())):
        panel = imageset.get_detector()[i_panel]
        sel = (panel_numbers == i_panel)
        isel = sel.iselection()
        ref_panel = refined_reflections.select(panel_numbers == i_panel)
        xy_cal_px.set_selected(
          sel, panel.millimeter_to_pixel(xy_cal_mm.select(sel)))
      x_px, y_px = xy_cal_px.parts()
      scan = imageset.get_scan()
      if scan is not None:
        z_px = scan.get_array_index_from_angle(z_rad, deg=False)
      else:
        # must be a still image, z centroid not meaningful
        z_px = z_rad
      xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
      self.refined_reflections['xyzcal.px'].set_selected(imgset_sel, xyzcal_px)

  def show_experiments(self, experiments, reflections, d_min=None):
    if d_min is not None:
      reciprocal_lattice_points = reflections['rlp']
      d_spacings = 1/reciprocal_lattice_points.norms()
      reflections = reflections.select(d_spacings > d_min)
    for i_expt, expt in enumerate(experiments):
      logger.info("model %i (%i reflections):" %(
        i_expt+1, (reflections['id'] == i_expt).count(True)))
      logger.info(expt.crystal)

    indexed_flags = reflections.get_flags(reflections.flags.indexed)
    imageset_id = reflections['imageset_id']
    rows = [['Imageset', '#indexed', '#unindexed']]
    for i in range(flex.max(imageset_id)+1):
      imageset_indexed_flags = indexed_flags.select(imageset_id == i)
      rows.append([str(i), str(imageset_indexed_flags.count(True)),
                   str(imageset_indexed_flags.count(False))])
    from libtbx import table_utils
    logger.info(
      table_utils.format(rows, has_header=True, prefix='| ', postfix=' |'))

  def find_max_cell(self):
    params = self.params.max_cell_estimation
    if self.params.max_cell is libtbx.Auto:
      if self.params.known_symmetry.unit_cell is not None:
        uc_params = self.target_symmetry_primitive.unit_cell().parameters()
        self.params.max_cell = params.multiplier * max(uc_params[:3])
        logger.info("Using max_cell: %.1f Angstrom" %(self.params.max_cell))
      else:
        self.params.max_cell = find_max_cell(
          self.reflections, max_cell_multiplier=params.multiplier,
          step_size=params.step_size,
          nearest_neighbor_percentile=params.nearest_neighbor_percentile,
          histogram_binning=params.histogram_binning,
          nn_per_bin=params.nn_per_bin,
          max_height_fraction=params.max_height_fraction,
          filter_ice=params.filter_ice,
          filter_overlaps=params.filter_overlaps,
          overlaps_border=params.overlaps_border).max_cell
        logger.info("Found max_cell: %.1f Angstrom" %(self.params.max_cell))

  def filter_reflections_by_scan_range(self):
    if len(self.params.scan_range):
      self.reflections = filter_reflections_by_scan_range(
        self.reflections, self.params.scan_range)

  @staticmethod
  def calculate_entering_flags(reflections, beam, goniometer):
    if goniometer is None:
      return flex.bool(len(reflections), False)
    axis = matrix.col(goniometer.get_rotation_axis())
    s0 = matrix.col(beam.get_s0())
    # calculate a unit vector normal to the spindle-beam plane for this
    # experiment, such that the vector placed at the centre of the Ewald sphere
    # points to the hemisphere in which reflections cross from inside to outside
    # of the sphere (reflections are exiting). NB this vector is in +ve Y
    # direction when using imgCIF coordinate frame.
    vec = s0.cross(axis)
    entering = reflections['s1'].dot(vec) < 0.
    return entering

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
      # e.g. data imported from XDS; no variance known then; since is used
      # only for weights assign as 1 => uniform weights
      if not 'xyzobs.px.variance' in spots_panel:
        spots_panel['xyzobs.px.variance'] = flex.vec3_double(
            len(spots_panel), (1,1,1))
      centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
        detector[i_panel], scan,
        spots_panel['xyzobs.px.value'],
        spots_panel['xyzobs.px.variance'],
        flex.vec3_double(len(spots_panel), (1,1,1)))
      spots['xyzobs.mm.value'].set_selected(sel, centroid_position)
      spots['xyzobs.mm.variance'].set_selected(sel, centroid_variance)
    return spots

  @staticmethod
  def map_centroids_to_reciprocal_space(spots_mm, detector, beam, goniometer,
                                        calculated=False):
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

    if 's1' not in spots_mm: spots_mm['s1'] = flex.vec3_double(len(spots_mm))
    spots_mm['rlp'] = flex.vec3_double(len(spots_mm))
    panel_numbers = flex.size_t(spots_mm['panel'])
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      spots_panel = spots_mm.select(panel_numbers == i_panel)
      if calculated:
        x, y, rot_angle = spots_panel['xyzcal.mm'].parts()
      else:
        x, y, rot_angle = spots_panel['xyzobs.mm.value'].parts()
      s1 = detector[i_panel].get_lab_coord(flex.vec2_double(x,y))
      s1 = s1/s1.norms() * (1/beam.get_wavelength())
      spots_mm['s1'].set_selected(sel, s1)
      S = s1 - beam.get_s0()
      if goniometer is not None:
        setting_rotation = matrix.sqr(goniometer.get_setting_rotation())
        rotation_axis = goniometer.get_rotation_axis_datum()
        fixed_rotation = matrix.sqr(goniometer.get_fixed_rotation())
        spots_mm['rlp'].set_selected(
          sel, tuple(setting_rotation.inverse()) * S)
        spots_mm['rlp'].set_selected(
          sel, spots_mm['rlp'].select(sel).rotate_around_origin(
            rotation_axis, -rot_angle))
        spots_mm['rlp'].set_selected(
          sel, tuple(fixed_rotation.inverse()) * spots_mm['rlp'].select(sel))
      else:
        spots_mm['rlp'].set_selected(sel, S)

  def find_candidate_orientation_matrices(self, candidate_basis_vectors):
    candidate_crystal_models = []
    vectors = candidate_basis_vectors
    if (self.target_symmetry_primitive is not None and
        self.target_symmetry_primitive.unit_cell() is not None):
      target_min_cell = self.target_symmetry_primitive.minimum_cell().unit_cell()
    else:
      target_min_cell = None

    # select unique combinations of input vectors to test
    # the order of combinations is such that combinations comprising vectors
    # nearer the beginning of the input list will appear before combinations
    # comprising vectors towards the end of the list
    n = len(vectors)
    # hardcoded limit on number of vectors, fixes issue #72
    # https://github.com/dials/dials/issues/72
    n = min(n, 100)
    vectors = vectors[:n]
    combinations = flex.vec3_int(flex.nested_loop((n,n,n)))
    combinations = combinations.select(
      flex.sort_permutation(combinations.as_vec3_double().norms()))

    from rstbx.indexing_api import tools
    transformations = [
      sgtbx.change_of_basis_op(
        str(sgtbx.rt_mx(sgtbx.rot_mx(T['trans'].transpose().as_int()))))
      for T in tools.R if T['mod'] < 5]
    transformations = [] # XXX temporarily disable cell doubling checks
    transformations.insert(0, sgtbx.change_of_basis_op())

    # select only those combinations where j > i and k > j
    i, j, k = combinations.as_vec3_double().parts()
    sel = flex.bool(len(combinations), True)
    sel &= j > i
    sel &= k > j
    combinations = combinations.select(sel)

    max_combinations = self.params.basis_vector_combinations.max_combinations
    if max_combinations is not None and max_combinations < len(combinations):
      combinations = combinations[:max_combinations]

    half_pi = 0.5 * math.pi
    min_angle = 20/180 * math.pi # 20 degrees, arbitrary cutoff
    for i, j, k in combinations:
      a = vectors[i]
      b = vectors[j]
      angle = a.angle(b)
      if angle < min_angle or (math.pi-angle) < min_angle:
        continue
      a_cross_b = a.cross(b)
      gamma = a.angle(b)
      if gamma < half_pi:
        # all angles obtuse if possible please
        b = -b
        gamma = math.pi - gamma
        a_cross_b = -a_cross_b
      c = vectors[k]
      if abs(half_pi-a_cross_b.angle(c)) < min_angle:
        continue
      alpha = b.angle(c, deg=True)
      if alpha < half_pi:
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
      if self.target_symmetry_primitive is not None:
        max_delta = self.params.known_symmetry.max_delta
        from dials.algorithms.indexing.symmetry import find_matching_symmetry
        for T in transformations:
          uc = model.get_unit_cell()
          det = T.c().r().determinant()
          if self.target_symmetry_primitive.unit_cell() is None:
            if det > 1: break
          else:
            primitive_volume = self.target_symmetry_primitive.unit_cell().volume()
          if det > 1 and abs(uc.volume()/primitive_volume - det) < 1e-1:
            uc = uc.change_basis(T)
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
          cb_op_to_primitive = cb_op_ref_primitive * cb_op_best_ref * cb_op_inp_best * T
          if cb_op_extra is not None:
            cb_op_to_primitive = cb_op_extra * cb_op_to_primitive
          best_model = model.change_basis(cb_op_to_primitive)

          def is_similar_cell(uc1, uc2):
            return uc1.is_similar_to(
              uc2,
              relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
              absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance)

          if (self.target_symmetry_primitive.unit_cell() is not None
              and not (
                is_similar_cell(best_model.get_unit_cell(),
                                self.target_symmetry_primitive.unit_cell())
                or
                is_similar_cell(
                  best_model.get_unit_cell().minimum_cell(), target_min_cell))):
            best_model = None
            continue
          else:
            break
      else:
        best_model = model

      if best_model is None:
        continue

      uc = best_model.get_unit_cell()
      params = uc.parameters()
      if uc.volume() > (params[0]*params[1]*params[2]/100):
        # unit cell volume cutoff from labelit 2004 paper
        candidate_crystal_models.append(best_model)
        if len(candidate_crystal_models) == self.params.basis_vector_combinations.max_refine:
          return candidate_crystal_models
    return candidate_crystal_models

  def choose_best_orientation_matrix(self, candidate_orientation_matrices):

    solution_scorer = self.params.basis_vector_combinations.solution_scorer
    if solution_scorer == 'weighted':
      weighted_params = self.params.basis_vector_combinations.weighted
      solutions = SolutionTrackerWeighted(
        power=weighted_params.power,
        volume_weight=weighted_params.volume_weight,
        n_indexed_weight=weighted_params.n_indexed_weight,
        rmsd_weight=weighted_params.rmsd_weight)
    else:
      filter_params = self.params.basis_vector_combinations.filter
      solutions = SolutionTrackerFilter(
        check_doubled_cell=filter_params.check_doubled_cell,
        likelihood_cutoff=filter_params.likelihood_cutoff,
        volume_cutoff=filter_params.volume_cutoff,
        n_indexed_cutoff=filter_params.n_indexed_cutoff)

    def run_one_refinement(args):
      params, reflections, experiments = args
      indexed_reflections = reflections.select(reflections['id'] > -1)

      from dials.command_line import check_indexing_symmetry
      grid_search_scope = params.indexing.check_misindexing.grid_search_scope

      best_offset = (0,0,0)
      best_cc = 0.0
      best_nref = 0

      if grid_search_scope > 0:
        offsets, ccs, nref \
          = check_indexing_symmetry.get_indexing_offset_correlation_coefficients(
          indexed_reflections, experiments.crystals()[0],
          grid=grid_search_scope, map_to_asu=True)

        if len(offsets) > 1:
          max_nref = flex.max(nref)

          # select "best" solution - needs nref > 0.5 max nref && highest CC
          # FIXME perform proper statistical test in here do not like heuristics

          for offset, cc, n in zip(offsets, ccs, nref):
            if n < (max_nref // 2):
              continue
            if cc > best_cc:
              best_cc = cc
              best_offset = offset
              best_nref = n

          #print offsets[13], nref[13], '%.2f' %ccs[13] # (0,0,0)
          #print best_offset, best_nref, '%.2f' %best_cc

          if best_offset != (0,0,0):
            logger.debug('Applying h,k,l offset: (%i, %i, %i)' %best_offset
                 + ' [cc = %.2f]' %best_cc)
            indexed_reflections['miller_index'] = apply_hkl_offset(
              indexed_reflections['miller_index'], best_offset)

      from dials.algorithms.refinement import RefinerFactory
      reflogger = logging.getLogger('dials.algorithms.refinement')
      level = reflogger.getEffectiveLevel()
      reflogger.setLevel(logging.ERROR)
      try:
        refiner = RefinerFactory.from_parameters_data_experiments(
          params, indexed_reflections, experiments)
        refiner.run()
      except (RuntimeError, ValueError, Sorry) as e:
        return
      else:
        rmsds = refiner.rmsds()
        xy_rmsds = math.sqrt(rmsds[0]**2 + rmsds[1]**2)
        model_likelihood = 1.0 - xy_rmsds
        soln = Solution(model_likelihood=model_likelihood,
                        crystal=experiments.crystals()[0],
                        rmsds=rmsds,
                        n_indexed=len(indexed_reflections),
                        fraction_indexed=float(len(indexed_reflections))/len(reflections),
                        hkl_offset=best_offset)
        return soln
      finally:
        reflogger.setLevel(level)

    import copy
    params = copy.deepcopy(self.all_params)
    params.refinement.parameterisation.auto_reduction.action = "fix"
    params.refinement.parameterisation.scan_varying = False
    params.refinement.refinery.max_iterations = 4
    params.refinement.reflections.reflections_per_degree = min(
      params.refinement.reflections.reflections_per_degree, 20)
    if params.refinement.reflections.outlier.block_width is libtbx.Auto:
      # auto block_width determination is potentially too expensive to do at
      # this stage: instead set separate_blocks=False and increase value
      # of tukey.iqr_multiplier to be more tolerant of outliers
      params.refinement.reflections.outlier.separate_blocks = False
      params.refinement.reflections.outlier.tukey.iqr_multiplier = \
        2 * params.refinement.reflections.outlier.tukey.iqr_multiplier

    args = []

    from dials.algorithms.indexing.compare_orientation_matrices \
         import difference_rotation_matrix_axis_angle
    for cm in candidate_orientation_matrices:
      sel = (self.reflections['id'] == -1)
      if self.d_min is not None:
        sel &= (1/self.reflections['rlp'].norms() > self.d_min)
      xo, yo, zo = self.reflections['xyzobs.mm.value'].parts()
      imageset_id = self.reflections['imageset_id']
      experiments = ExperimentList()
      for i_imageset, imageset in enumerate(self.imagesets):
        scan = imageset.get_scan()
        if scan is not None:
          start, end = scan.get_oscillation_range()
          if (end - start) > 360:
            # only use reflections from the first 360 degrees of the scan
            sel.set_selected(
              (imageset_id == i_imageset) & (zo > ((start * math.pi/180) + 2 * math.pi)), False)
        experiments.append(Experiment(imageset=imageset,
                                      beam=imageset.get_beam(),
                                      detector=imageset.get_detector(),
                                      goniometer=imageset.get_goniometer(),
                                      scan=scan,
                                      crystal=cm))
      refl = self.reflections.select(sel)
      self.index_reflections(experiments, refl)
      if refl.get_flags(refl.flags.indexed).count(True) == 0:
        continue

      from rstbx.dps_core.cell_assessment import SmallUnitCellVolume
      threshold = self.params.basis_vector_combinations.sys_absent_threshold
      if (threshold and (self.target_symmetry_primitive is None or
                         self.target_symmetry_primitive.unit_cell() is None)):
        try:
          self.correct_non_primitive_basis(experiments, refl, threshold)
          if refl.get_flags(refl.flags.indexed).count(True) == 0:
            continue
        except SmallUnitCellVolume:
          logger.debug("correct_non_primitive_basis SmallUnitCellVolume error for unit cell %s:"
                       %experiments[0].crystal.get_unit_cell())
          continue
        except RuntimeError as e:
          if 'Krivy-Gruber iteration limit exceeded' in str(e):
            logger.debug("correct_non_primitive_basis Krivy-Gruber iteration limit exceeded error for unit cell %s:"
                         %experiments[0].crystal.get_unit_cell())
            continue
          raise
        if experiments[0].crystal.get_unit_cell().volume() < self.params.min_cell_volume:
          continue

      if self.params.known_symmetry.space_group is not None:
        target_space_group = self.target_symmetry_primitive.space_group()
        new_crystal, cb_op_to_primitive = self.apply_symmetry(
          experiments[0].crystal, target_space_group)
        if new_crystal is None:
          continue
        experiments[0].crystal.update(new_crystal)
        if not cb_op_to_primitive.is_identity_op():
          sel = refl['id'] > -1
          miller_indices = refl['miller_index'].select(sel)
          miller_indices = cb_op_to_primitive.apply(miller_indices)
          refl['miller_index'].set_selected(sel, miller_indices)
        if 0 and self.cb_op_primitive_to_given is not None:
          sel = refl['id'] > -1
          experiments[0].crystal.update(
            experiments[0].crystal.change_basis(self.cb_op_primitive_to_given))
          miller_indices = refl['miller_index'].select(sel)
          miller_indices = self.cb_op_primitive_to_given.apply(miller_indices)
          refl['miller_index'].set_selected(sel, miller_indices)

      if (self.refined_experiments is not None and
          len(self.refined_experiments) > 0):
        orientation_too_similar = False
        cryst_b = experiments[0].crystal
        for i_a, cryst_a in enumerate(self.refined_experiments.crystals()):
          R_ab, axis, angle, cb_op_ab = \
            difference_rotation_matrix_axis_angle(cryst_a, cryst_b)
          min_angle = self.params.multiple_lattice_search.minimum_angular_separation
          if abs(angle) < min_angle: # degrees
            orientation_too_similar = True
            break
        if orientation_too_similar:
          logger.debug("skipping crystal: too similar to other crystals")
          continue

      args.append((params, refl, experiments))

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

    if len(solutions):
      logger.info("Candidate solutions:")
      logger.info(str(solutions))
      best_solution = solutions.best_solution()
      logger.debug("best model_likelihood: %.2f" %best_solution.model_likelihood)
      logger.debug("best n_indexed: %i" %best_solution.n_indexed)
      self.hkl_offset = best_solution.hkl_offset
      return best_solution.crystal, best_solution.n_indexed
    else:
      return None, None

  def correct_non_primitive_basis(self, experiments, reflections, threshold):
    assert len(experiments.crystals()) == 1
    while True:
      sel = reflections['miller_index'] != (0,0,0)
      if sel.count(True) == 0:
        break
      T = detect_non_primitive_basis(
        reflections['miller_index'].select(sel), threshold=threshold)
      if T is None:
        break

      crystal_model = experiments.crystals()[0]
      direct_matrix = matrix.sqr(crystal_model.get_A()).inverse()
      M = T.inverse().transpose()
      new_direct_matrix = M * direct_matrix
      crystal_model.set_A(new_direct_matrix.inverse())

      from rstbx.dps_core.cell_assessment import unit_cell_too_small
      unit_cell_too_small(crystal_model.get_unit_cell())
      cb_op = crystal_model.get_unit_cell().change_of_basis_op_to_niggli_cell()
      crystal_model.update(crystal_model.change_basis(cb_op))

      reflections['id'] = flex.int(len(reflections), -1)
      reflections.unset_flags(flex.bool(len(reflections), True),
        reflections.flags.indexed)
      self.index_reflections(experiments, reflections)

  def apply_symmetry(self, crystal_model, target_space_group):
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
    assert target_sg_best.is_compatible_unit_cell(model.get_unit_cell())

    model = model.change_basis(cb_op_best_primitive)
    return model, cb_op_inp_primitive

  def index_reflections(self, experiments, reflections):
    if self.params.index_assignment.method == 'local':
      params_local = self.params.index_assignment.local
      from dials.algorithms.indexing import index_reflections_local
      index_reflections_local(
        reflections,
        experiments, self.d_min, epsilon=params_local.epsilon,
        delta=params_local.delta, l_min=params_local.l_min,
        nearest_neighbours=params_local.nearest_neighbours)
    else:
      params_simple = self.params.index_assignment.simple
      from dials.algorithms.indexing import index_reflections
      index_reflections(reflections,
                        experiments, self.d_min,
                        tolerance=params_simple.hkl_tolerance)
    if self.hkl_offset is not None and self.hkl_offset != (0,0,0):
      reflections['miller_index'] = apply_hkl_offset(
        reflections['miller_index'], self.hkl_offset)
      self.hkl_offset = None

  def refine(self, experiments, reflections):
    from dials.algorithms.indexing.refinement import refine
    refiner, refined, outliers = refine(
      self.all_params, reflections, experiments,
      verbosity=self.params.refinement_protocol.verbosity,
      debug_plots=self.params.debug_plots)
    if outliers is not None:
      reflections['id'].set_selected(outliers, -1)
    predicted = refiner.predict_for_indexed()
    verbosity = self.params.refinement_protocol.verbosity
    reflections['xyzcal.mm'] = predicted['xyzcal.mm']
    reflections['entering'] = predicted['entering']
    reflections.unset_flags(
      flex.bool(len(reflections), True), reflections.flags.centroid_outlier)
    assert reflections.get_flags(reflections.flags.centroid_outlier).count(True) == 0
    reflections.set_flags(predicted.get_flags(predicted.flags.centroid_outlier),
                          reflections.flags.centroid_outlier)
    reflections.set_flags(refiner.selection_used_for_refinement(),
                          reflections.flags.used_in_refinement)
    return refiner.get_experiments(), reflections

  def debug_show_candidate_basis_vectors(self):

    vectors = self.candidate_basis_vectors

    logger.debug("Candidate basis vectors:")
    for i, v in enumerate(vectors):
      logger.debug("%s %s" %(i, v.length()))# , vector_heights[i]

    if self.params.debug:
      # print a table of the angles between each pair of vectors

      from cStringIO import StringIO
      s = StringIO()

      angles = flex.double(len(vectors)**2)
      angles.reshape(flex.grid(len(vectors), len(vectors)))

      for i in range(len(vectors)):
        v_i = vectors[i]
        for j in range(i+1, len(vectors)):
          v_j = vectors[j]
          angles[i,j] = v_i.angle(v_j, deg=True)

      print((" "*7), end=' ', file=s)
      for i in range(len(vectors)):
        print("%7.3f" % vectors[i].length(), end=' ', file=s)
      print(file=s)
      for i in range(len(vectors)):
        print("%7.3f" % vectors[i].length(), end=' ', file=s)
        for j in range(len(vectors)):
          if j <= i:
            print((" "*7), end=' ', file=s)
          else:
            print("%5.1f  " %angles[i,j], end=' ', file=s)
        print(file=s)

      logger.debug(s.getvalue())

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
    for i_panel in range(len(self.imagesets[0].get_detector())):
      if len(self.imagesets[0].get_detector()) > 1:
        file_name = 'reciprocal_lattice_%i.pdb' %i_panel
      with open(file_name, 'wb') as f:
        xs = xray.structure(crystal_symmetry=cs)
        reflections = self.reflections.select(self.reflections['panel'] == i_panel)
        for site in reflections['rlp']:
          xs.add_scatterer(xray.scatterer("C", site=site))
        xs.sites_mod_short()
        f.write(xs.as_pdb_file())

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
    reflections.as_pickle(file_name)

  def find_lattices(self):
    raise NotImplementedError()


from libtbx import group_args
class Solution(group_args):
  pass


def filter_doubled_cell(solutions):
  from dials.algorithms.indexing.compare_orientation_matrices import difference_rotation_matrix_axis_angle
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
          R, axis, angle, cb = difference_rotation_matrix_axis_angle(new_cryst, s2.crystal)
          if ((angle < 1) and
              (s1.n_indexed < (1.1 * s2.n_indexed))):
            doubled_cell = True
            break

    if not doubled_cell:
      accepted_solutions.append(s1)

  return accepted_solutions


# Tracker for solutions based on code in rstbx/dps_core/basis_choice.py
class SolutionTrackerFilter(object):
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

    if self.check_doubled_cell:
      self.filtered_solutions = filter_doubled_cell(self.filtered_solutions)

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

  def __str__(self):
    rows = []
    rows.append(
      ['unit_cell', 'volume', 'n_indexed', 'fraction_indexed', 'likelihood'])

    for i, s in enumerate(self.all_solutions):
      s = self.all_solutions[i]
      rows.append(
        [format(s.crystal.get_unit_cell(), '{:.2f} {:.2f} {:.2f} {:.1f} {:.1f} {:.1f}'),
         "%.0f" %s.crystal.get_unit_cell().volume(),
         str(s.n_indexed), "%.0f" %(s.fraction_indexed*100), "%.2f" %s.model_likelihood])

    from libtbx import table_utils
    return table_utils.format(rows=rows, has_header=True)


class SolutionTrackerWeighted(object):
  def __init__(self, power=2, volume_weight=1, n_indexed_weight=1, rmsd_weight=1):
    self.volume_weight = volume_weight
    self.n_indexed_weight = n_indexed_weight
    self.rmsd_weight = rmsd_weight
    self.power = power
    self.all_solutions = []

  def append(self, item):
    self.all_solutions.append(item)

  def __len__(self):
    return len(self.all_solutions)

  def score_by_volume(self, reverse=False):
    # smaller volume = better
    volumes = flex.double(
      s.crystal.get_unit_cell().volume() for s in self.all_solutions)
    score = flex.log(volumes)/math.log(2)
    return self.volume_weight * (score - flex.min(score))

  def score_by_rmsd_xy(self, reverse=False):
    # smaller rmsds = better
    rmsd_x, rmsd_y, rmsd_z = flex.vec3_double(
      s.rmsds for s in self.all_solutions).parts()
    rmsd_xy = flex.sqrt(flex.pow2(rmsd_x) + flex.pow2(rmsd_y))
    score = flex.log(rmsd_xy)/math.log(2)
    return self.rmsd_weight * (score - flex.min(score))

  def score_by_rmsd_z(self, reverse=False):
    # smaller rmsds = better
    rmsd_x, rmsd_y, rmsd_z = flex.vec3_double(
      s.rmsds for s in self.all_solutions).parts()
    score = flex.log(rmsd_z)/math.log(2)
    return score - flex.min(score)

  def score_by_fraction_indexed(self, reverse=False):
    # more indexed reflections = better
    fraction_indexed = flex.double(s.fraction_indexed for s in self.all_solutions)
    score = flex.log(fraction_indexed)/math.log(2)
    return self.n_indexed_weight * (-score + flex.max(score))

  def best_solution(self):
    scores = self.solution_scores()
    perm = flex.sort_permutation(scores)
    return self.all_solutions[perm[0]]

  def solution_scores(self):
    scores = sum(flex.pow(score.as_double(), self.power)
                 for score in (
                               self.score_by_fraction_indexed(),
                               self.score_by_volume(),
                               self.score_by_rmsd_xy(),
                               #self.score_by_rmsd_z(),
                               ))
    return scores

  def __str__(self):
    rows = []
    rows.append(
      ['unit_cell', 'volume', 'volume score',
       '#indexed', '% indexed', '% indexed score',
       'rmsd_xy', 'rmsd_xy score',
       #'rmsd_z', 'rank_rmsd_z',
       'overall score'])

    score_by_fraction_indexed = self.score_by_fraction_indexed()
    score_by_volume = self.score_by_volume()
    score_by_rmsd_xy = self.score_by_rmsd_xy()
    score_by_rmsd_z = self.score_by_rmsd_z()
    overall_scores = self.solution_scores()

    perm = flex.sort_permutation(overall_scores)

    rmsd_x, rmsd_y, rmsd_z = flex.vec3_double(
      s.rmsds for s in self.all_solutions).parts()
    rmsd_xy = flex.sqrt(flex.pow2(rmsd_x) + flex.pow2(rmsd_y))
    rmsd_z *= (180/math.pi)

    for i in perm:
      s = self.all_solutions[i]
      rows.append(
        [format(s.crystal.get_unit_cell(), '{:.2f} {:.2f} {:.2f} {:.1f} {:.1f} {:.1f}'),
         "%.0f" %s.crystal.get_unit_cell().volume(), "%.2f" %score_by_volume[i],
         str(s.n_indexed), "%.0f" %(s.fraction_indexed*100),
         "%.2f" %score_by_fraction_indexed[i],
         "%.2f" %rmsd_xy[i], "%.2f" %score_by_rmsd_xy[i],
         #"%.2f" %rmsd_z[i], "%.2f" %score_by_rmsd_z[i],
         "%.2f" %overall_scores[i]])

    from libtbx import table_utils
    return table_utils.format(rows=rows, has_header=True)


def detect_non_primitive_basis(miller_indices, threshold=0.9):

  from rstbx.indexing_api import tools
  for test in tools.R:
    cum = tools.cpp_absence_test(miller_indices, test['mod'], test['vec'])
    for counter in xrange(test['mod']):
      if float(cum[counter])/miller_indices.size() > threshold and counter==0:
        # (if counter != 0 there is no obvious way to correct this)
        logger.debug("Detected exclusive presence of %dH %dK %dL = %dn, remainder %d"%(
          test['vec'][0], test['vec'][1], test['vec'][2], test['mod'],counter))
        logger.debug("%s, %s, %s" %(
          test['vec'], test['mod'], float(cum[counter])/miller_indices.size()))
        #flag = {'vec':test['vec'],'mod':test['mod'],
                #'remainder':counter, 'trans':test['trans'].elems}
        return test['trans']


def find_max_cell(reflections, max_cell_multiplier, step_size,
                  nearest_neighbor_percentile=None, histogram_binning='linear',
                  nn_per_bin=5, max_height_fraction=0.25,
                  filter_ice=True, filter_overlaps=True, overlaps_border=0):
  logger.debug('Finding suitable max_cell based on %i reflections' % len(reflections))
  # Exclude potential ice-ring spots from nearest neighbour analysis if needed
  if filter_ice:
    from dials.algorithms.spot_finding.per_image_analysis import \
     ice_rings_selection
    ice_sel = ice_rings_selection(reflections)
    reflections = reflections.select(~ice_sel)
    logger.debug('Rejecting %i reflections at ice ring resolution' %ice_sel.count(True))

  # need bounding box in reflections to find overlaps; this is not there if
  # spots are from XDS (for example)
  if filter_overlaps and 'bbox' in reflections:
    overlap_sel = flex.bool(len(reflections), False)

    overlaps = reflections.find_overlaps(border=overlaps_border)
    for item in overlaps.edges():
      i0 = overlaps.source(item)
      i1 = overlaps.target(item)
      overlap_sel[i0] = True
      overlap_sel[i1] = True
    logger.debug('Rejecting %i overlapping bounding boxes' %overlap_sel.count(True))
    reflections = reflections.select(~overlap_sel)
  logger.debug('%i reflections remain for max_cell identification' % len(reflections))

  # Histogram spot counts in resolution bins: filter out outlier bin counts
  # according to the Tukey method
  d_spacings = 1/reflections['rlp'].norms()
  d_star_sq = 1/flex.pow2(d_spacings)
  hist = flex.histogram(d_star_sq, n_slots=100)
  bin_counts = hist.slots().as_double()
  perm = flex.sort_permutation(bin_counts)
  sorted_bin_counts = bin_counts.select(perm)
  sorted_bin_counts = sorted_bin_counts.select(sorted_bin_counts > 0)
  if sorted_bin_counts.size() >=4:
    from libtbx.math_utils import nearest_integer as nint
    q1 = sorted_bin_counts[nint(len(sorted_bin_counts)/4)]
    q2 = sorted_bin_counts[nint(len(sorted_bin_counts)/2)]
    q3 = sorted_bin_counts[nint(len(sorted_bin_counts)*3/4)]
    iqr = q3-q1
    inlier_sel = (bin_counts >= (q1 - 1.5*iqr)) & (bin_counts <= (q3 + 1.5*iqr))

    sel = flex.bool(d_star_sq.size(), True)
    for i_slot in range(hist.slots().size()):
      if not inlier_sel[i_slot]:
        sel.set_selected(
          (d_star_sq > (hist.slot_centers()[i_slot] - 0.5*hist.slot_width())) &
          (d_star_sq <= (hist.slot_centers()[i_slot] + 0.5*hist.slot_width())), False)

    reflections = reflections.select(sel)
  #from libtbx import easy_pickle
  #easy_pickle.dump('filtered.pickle', reflections)

  assert len(reflections) > 0, \
         'Too few spots remaining for nearest neighbour analysis (%d)' %len(reflections)
  # The nearest neighbour analysis gets fooled when the same part of
  # reciprocal space has been measured twice as this introduced small
  # random differences in position between reflections measured twice.
  # Therefore repeat the nearest neighbour analysis several times in small
  # wedges where there shouldn't be any overlap in reciprocal space
  #from rstbx.indexing_api.nearest_neighbor import neighbor_analysis
  if 'entering' in reflections:
    entering_flags = reflections['entering']
  else:
    entering_flags = flex.bool(len(reflections), False)
  from dials.algorithms.indexing.nearest_neighbor import neighbor_analysis
  phi_deg = reflections['xyzobs.mm.value'].parts()[2] * (180/math.pi)

  NN = neighbor_analysis(reflections,
                         step_size=step_size,
                         max_height_fraction=max_height_fraction,
                         tolerance=max_cell_multiplier,
                         percentile=nearest_neighbor_percentile,
                         histogram_binning=histogram_binning,
                         nn_per_bin=nn_per_bin)

  return NN


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

def apply_hkl_offset(indices, offset):
  h, k, l = indices.as_vec3_double().parts()
  h += offset[0]
  k += offset[1]
  l += offset[2]
  return flex.miller_index(
    h.iround(), k.iround(), l.iround())
