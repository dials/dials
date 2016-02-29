#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

from __future__ import division
from dials.array_family import flex
from logging import info, debug
from dials.util import log
debug_handle = log.debug_handle()
info_handle = log.info_handle()
import libtbx
from libtbx.utils import Sorry
from dials.algorithms.indexing.indexer import indexer_base
from dials.algorithms.indexing.known_orientation import indexer_known_orientation
from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search
from dials.algorithms.indexing.fft3d import indexer_fft3d
from dials.algorithms.indexing.fft1d import indexer_fft1d
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList
from dials.algorithms.indexing.nave_parameters import nave_parameters
import math
from dials.algorithms.indexing.indexer import master_params

def calc_2D_rmsd_and_displacements(reflections):

  displacements = flex.vec2_double(reflections['xyzobs.px.value'].parts()[0], reflections['xyzobs.px.value'].parts()[1]) - \
                  flex.vec2_double(reflections['xyzcal.px'].parts()[0], reflections['xyzcal.px'].parts()[1])
  rmsd = math.sqrt(flex.mean(displacements.dot( displacements )))

  return rmsd,displacements

def plot_displacements(reflections, predictions, experiments):
    rmsd,displacements = calc_2D_rmsd_and_displacements(predictions)

    if True: #params.refinement.plot_residual_scatter:
      from matplotlib import pyplot as plt
      fig = plt.figure()
      for cv in displacements:
        plt.plot([cv[0]],[-cv[1]],"r.")
      plt.title(" %d spots, r.m.s.d. %5.2f pixels"%(len(displacements),rmsd))
      plt.axes().set_aspect("equal")
      plt.show()
      plt.close()

    if True: #params.refinement.plot_residual_map:
      from matplotlib import pyplot as plt
      PX = reflections["xyzobs.px.value"]
      fig = plt.figure()
      sz1,sz2 = experiments[0].detector[0].get_image_size()
      for item,cv in zip(predictions,displacements):
        plt.plot([item['xyzcal.px'][0]],[sz1 - item['xyzcal.px'][1]],"r.")
        plt.plot([item['xyzobs.px.value'][0]],[sz1 - item['xyzobs.px.value'][1]],"g.")
        plt.plot([item['xyzcal.px'][0], item['xyzcal.px'][0] + 10.*cv[0]],
                 [sz1 - item['xyzcal.px'][1], sz1 - item['xyzcal.px'][1] - 10.*cv[1]],'r-')

      plt.xlim([0,experiments[0].detector[0].get_image_size()[0]])
      plt.ylim([0,experiments[0].detector[0].get_image_size()[1]])
      plt.title(" %d spots, r.m.s.d. %5.2f pixels"%(len(displacements),rmsd))
      plt.axes().set_aspect("equal")
      plt.show()
      plt.close()

def e_refine(params, experiments, reflections, graph_verbose=False):
    # Stills-specific parameters we always want
    assert params.refinement.reflections.outlier.algorithm in (None, "null"), \
      "Cannot index, set refinement.reflections.outlier.algorithm=null" # we do our own outlier rejection

    from dials.algorithms.refinement.refiner import RefinerFactory
    refiner = RefinerFactory.from_parameters_data_experiments(params,
      reflections, experiments, verbosity=1)

    history = refiner.run()

    ref_sel = refiner.selection_used_for_refinement()
    assert ref_sel.count(True) == len(reflections)

    if not graph_verbose: return refiner

    RR = refiner.predict_for_reflection_table(reflections)
    plot_displacements(reflections, RR, experiments)

    return refiner


class stills_indexer(indexer_base):
  ''' Class for indexing stills '''

  @staticmethod
  def from_parameters(reflections, imagesets,
                      known_crystal_models=None, params=None):

    if params is None:
      params = master_params

    if known_crystal_models is not None:
      idxr = stills_indexer_known_orientation(
        reflections, imagesets, params, known_crystal_models)
    elif params.indexing.method == "fft3d":
      idxr = stills_indexer_fft3d(reflections, imagesets, params=params)
    elif params.indexing.method == "fft1d":
      idxr = stills_indexer_fft1d(reflections, imagesets, params=params)
    elif params.indexing.method == "real_space_grid_search":
      idxr = stills_indexer_real_space_grid_search(reflections, imagesets, params=params)

    return idxr

  def __init__(self, reflections, imagesets, params=None):
    if params.refinement.reflections.outlier.algorithm in ('auto', libtbx.Auto):
      # The stills_indexer provides it's own outlier rejection
      params.refinement.reflections.outlier.algorithm = 'null'
    indexer_base.__init__(self, reflections, imagesets, params)

  def index(self):
    # most of this is the same as dials.algorithms.indexing.indexer.indexer_base.index(), with some stills
    # specific modifications (don't re-index after choose best orientation matrix, but use the indexing from
    # choose best orientation matrix, also don't use macrocycles)
    # of refinement after indexing).
    if self.params.refinement_protocol.n_macro_cycles > 1:
      raise Sorry("For stills, please set refinement_protocol.n_macro_cycles = 1")

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

      self.imagesets[0].set_detector(new_detector)
      self.imagesets[0].set_beam(new_beam)

    spots_mm = self.reflections
    self.reflections = flex.reflection_table()

    if 'imageset_id' not in spots_mm:
      spots_mm['imageset_id'] = spots_mm['id']
    for i, imageset in enumerate(self.imagesets):
      spots_sel = spots_mm.select(spots_mm['imageset_id'] == i)
      self.map_centroids_to_reciprocal_space(
        spots_sel, imageset.get_detector(), imageset.get_beam(),
        imageset.get_goniometer())
      self.reflections.extend(spots_sel)

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

      if hasattr(self, '_best_indexed'):
        self.reflections = self._best_indexed
      else:
        self.index_reflections(
          experiments, self.reflections,
          verbosity=self.params.refinement_protocol.verbosity)

      if len(experiments) == n_lattices_previous_cycle:
        # no more lattices found
        break

      if (self.target_symmetry_primitive is not None
          and self.target_symmetry_primitive.space_group() is not None):
        # now apply the space group symmetry only after the first indexing
        # need to make sure that the symmetrized orientation is similar to the P1 model
        for i_cryst, cryst in enumerate(experiments.crystals()):
          if i_cryst >= n_lattices_previous_cycle:
            new_cryst, cb_op_to_primitive = self.apply_symmetry(
              cryst, self.target_symmetry_primitive,
              space_group_only=True)
            if self.cb_op_primitive_inp is not None:
              new_cryst = new_cryst.change_basis(self.cb_op_primitive_inp)
              info(new_cryst.get_space_group().info())
            cryst.update(new_cryst)
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

        info("")
        info("#" * 80)
        info("Starting refinement")
        info("#" * 80)
        info("")
        self.indexed_reflections = (self.reflections['id'] > -1)

        sel = flex.bool(len(self.reflections), False)
        lengths = 1/self.reflections['rlp'].norms()
        isel = (lengths >= self.d_min).iselection()
        sel.set_selected(isel, True)
        sel.set_selected(self.reflections['id'] > -1, False)
        self.unindexed_reflections = self.reflections.select(sel)

        reflections_for_refinement = self.reflections.select(
          self.indexed_reflections)
        try:
          refined_experiments, refined_reflections = self.refine(
            experiments, reflections_for_refinement)
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

        self.refined_reflections = refined_reflections.select(
          refined_reflections['id'] > -1)

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

    if not 'refined_experiments' in locals():
      raise Sorry("None of the experiments could refine.")

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import show_rotation_matrix_differences
      show_rotation_matrix_differences(
        self.refined_experiments.crystals(), out=info_handle)

    info("Final refined crystal models:")
    for i, crystal_model in enumerate(self.refined_experiments.crystals()):
      n_indexed = 0
      for i_expt in experiments.where(crystal=crystal_model):
        n_indexed += (self.reflections['id'] == i).count(True)
      info("model %i (%i reflections):" %(i+1, n_indexed))
      info(crystal_model)

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

  def experiment_list_for_crystal(self, crystal):
    experiments = ExperimentList()
    for imageset in self.imagesets:
      experiments.append(Experiment(imageset=imageset,
                                    beam=imageset.get_beam(),
                                    detector=imageset.get_detector(),
                                    goniometer=imageset.get_goniometer(),
                                    scan=imageset.get_scan(),
                                    crystal=crystal))
    return experiments

  def choose_best_orientation_matrix(self, candidate_orientation_matrices):
    from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList
    from logging import info
    import copy

    info('*' * 80)
    info('Selecting the best orientation matrix')
    info('*' * 80)

    from libtbx import group_args
    class candidate_info(group_args):
      pass
    candidates = []

    params = copy.deepcopy(self.all_params)

    for icm,cm in enumerate(candidate_orientation_matrices):
      sel = ((self.reflections['id'] == -1))
             #(1/self.reflections['rlp'].norms() > self.d_min))
      refl = self.reflections.select(sel)
      experiments = self.experiment_list_for_crystal(cm)
      self.index_reflections(experiments, refl)
      indexed = refl.select(refl['id'] >= 0)
      indexed = indexed.select(indexed.get_flags(indexed.flags.indexed))

      print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d initial outlier identification"%icm
      acceptance_flags = self.identify_outliers(params, experiments, indexed)
      #create a new "indexed" list with outliers thrown out:
      indexed = indexed.select(acceptance_flags)

      print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d refinement before outlier rejection"%icm
      R = e_refine(params = params, experiments=experiments, reflections=indexed, graph_verbose=False)
      ref_experiments = R.get_experiments()

      # try to improve the outcome with a second round of outlier rejection post-initial refinement:
      acceptance_flags = self.identify_outliers(params, ref_experiments, indexed)

      # insert a round of Nave-outlier rejection on top of the r.m.s.d. rejection
      nv0 = nave_parameters(params = params, experiments=ref_experiments, reflections=indexed, refinery=R, graph_verbose=False)
      crystal_model_nv0 = nv0()
      acceptance_flags_nv0 = nv0.nv_acceptance_flags
      indexed = indexed.select(acceptance_flags & acceptance_flags_nv0)

      print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d after positional and delta-psi outlier rejection"%icm
      R = e_refine(params = params, experiments=ref_experiments, reflections=indexed, graph_verbose=False)
      ref_experiments = R.get_experiments()

      nv = nave_parameters(params = params, experiments=ref_experiments, reflections=indexed, refinery=R, graph_verbose=False)
      crystal_model = nv()

      rmsd, _ = calc_2D_rmsd_and_displacements(R.predict_for_reflection_table(indexed))

      print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d done"%icm
      candidates.append(candidate_info(crystal = crystal_model,
                                       green_curve_area = nv.green_curve_area,
                                       ewald_proximal_volume = nv.ewald_proximal_volume(),
                                       n_indexed = len(indexed),
                                       rmsd = rmsd,
                                       indexed = indexed,
                                       experiments = ref_experiments))

    if len(candidates) == 0:
      raise Sorry("No suitable indexing solution found")

    print "**** ALL CANDIDATES:"
    for i,XX in enumerate(candidates):
      print "\n****Candidate %d"%i,XX
      cc = XX.crystal
      print "  half mosaicity %5.2f deg."%(cc._ML_half_mosaicity_deg)
      print "  domain size %.0f Ang."%(cc._ML_domain_size_ang)
    print "\n**** BEST CANDIDATE:"

    results = flex.double([c.rmsd for c in candidates])
    best = candidates[flex.min_index(results)]
    print best

    if best.rmsd > 1.5:
      raise Sorry ("RMSD too high, %f" %rmsd)

    if best.ewald_proximal_volume > 0.0015:
      raise Sorry ("Ewald proximity volume too high, %f"%best.ewald_proximal_volume)

    if len(candidates) > 1:
      for i in xrange(len(candidates)):
        if i == flex.min_index(results):
          continue
        if best.ewald_proximal_volume > candidates[i].ewald_proximal_volume:
          print "Couldn't figure out which candidate is best; picked the one with the best RMSD."

    best.indexed['entering'] = flex.bool(best.n_indexed, False)

    self._best_indexed = best.indexed
    return best.crystal, best.n_indexed

  def identify_outliers(self, params, experiments, indexed):
      print "$$$ stills_indexer::identify_outliers"
      refiner = e_refine(params, experiments, indexed, graph_verbose=False)

      RR = refiner.predict_for_reflection_table(indexed)

      px_sz = experiments[0].detector[0].get_pixel_size()

      class match: pass
      matches = []
      for item in RR:
        m = match()
        m.x_obs = item["xyzobs.px.value"][0]*px_sz[0]
        m.y_obs = item["xyzobs.px.value"][1]*px_sz[1]
        m.x_calc= item["xyzcal.px"][0]*px_sz[0]
        m.y_calc= item["xyzcal.px"][1]*px_sz[1]
        m.miller_index = item["miller_index"]
        matches.append(m)

      from rstbx.phil.phil_preferences import indexing_api_defs
      import iotbx.phil
      hardcoded_phil = iotbx.phil.parse(
      input_string=indexing_api_defs).extract()

      from rstbx.indexing_api.outlier_procedure import OutlierPlotPDF

      #comment this in if PDF graph is desired:
      #hardcoded_phil.indexing.outlier_detection.pdf = "outlier.pdf"
      # new code for outlier rejection inline here
      if hardcoded_phil.indexing.outlier_detection.pdf is not None:
        hardcoded_phil.__inject__("writer",OutlierPlotPDF(hardcoded_phil.indexing.outlier_detection.pdf))

      # execute Sauter and Poon (2010) algorithm
      from rstbx.indexing_api import outlier_detection
      od = outlier_detection.find_outliers_from_matches(
        matches,
        verbose=True,
        horizon_phil=hardcoded_phil)

      if hardcoded_phil.indexing.outlier_detection.pdf is not None:
        od.make_graphs(canvas=hardcoded_phil.writer.R.c,left_margin=0.5)
        hardcoded_phil.writer.R.c.showPage()
        hardcoded_phil.writer.R.c.save()

      return od.get_cache_status()

  def refine(self, experiments, reflections):

    sel = ((reflections['id'] >= -1))
    refl = reflections.select(sel)

    acceptance_flags = self.identify_outliers(self.all_params, experiments, refl)
    #create a new "indexed" list with outliers thrown out:
    refl = refl.select(acceptance_flags)

    print "$$$ stills_indexer::refine"
    R = e_refine(params = self.all_params, experiments=experiments, reflections=refl, graph_verbose=False)
    ref_experiments = R.get_experiments()

    # try to improve the outcome with a second round of outlier rejection post-initial refinement:
    acceptance_flags = self.identify_outliers(self.all_params, ref_experiments, refl)

    # insert a round of Nave-outlier rejection on top of the r.m.s.d. rejection
    nv0 = nave_parameters(params = self.all_params, experiments=ref_experiments, reflections=refl, refinery=R, graph_verbose=False)
    crystal_model_nv0 = nv0()
    acceptance_flags_nv0 = nv0.nv_acceptance_flags
    refl = refl.select(acceptance_flags & acceptance_flags_nv0)

    print "$$$ stills_indexer::refine after positional and delta-psi outlier rejection"
    refiner = e_refine(params = self.all_params, experiments=ref_experiments, reflections=refl, graph_verbose=False)

    matches = refiner.get_matches()
    xyzcal_mm = flex.vec3_double(len(refl))
    xyzcal_mm.set_selected(matches['iobs'], matches['xyzcal.mm'])
    refl['xyzcal.mm'] = xyzcal_mm
    refl.set_flags(matches['iobs'], refl.flags.used_in_refinement)
    refl['entering'] = flex.bool(len(refl), False)
    return refiner.get_experiments(), refl

""" Mixin class definitions that override the dials indexing class methods specific to stills """
class stills_indexer_known_orientation(indexer_known_orientation, stills_indexer):
  pass

class stills_indexer_real_space_grid_search(stills_indexer, indexer_real_space_grid_search):
  pass

class stills_indexer_fft3d(stills_indexer, indexer_fft3d):
  pass

class stills_indexer_fft1d(stills_indexer, indexer_fft1d):
  pass
