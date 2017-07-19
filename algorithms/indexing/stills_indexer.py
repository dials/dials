#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#

from __future__ import absolute_import, division
import logging
logger = logging.getLogger(__name__)

from dials.array_family import flex
from dials.util import log
debug_handle = log.debug_handle(logger)
info_handle = log.info_handle(logger)
import libtbx
from libtbx.utils import Sorry
from dials.algorithms.indexing.indexer import indexer_base
from dials.algorithms.indexing.known_orientation import indexer_known_orientation
from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search
from dials.algorithms.indexing.fft3d import indexer_fft3d
from dials.algorithms.indexing.fft1d import indexer_fft1d
from dxtbx.model.experiment_list import Experiment, ExperimentList
from dials.algorithms.indexing.nave_parameters import nave_parameters
import math
from dials.algorithms.indexing.indexer import master_params

def calc_2D_rmsd_and_displacements(reflections):

  displacements = flex.vec2_double(reflections['xyzobs.px.value'].parts()[0], reflections['xyzobs.px.value'].parts()[1]) - \
                  flex.vec2_double(reflections['xyzcal.px'].parts()[0], reflections['xyzcal.px'].parts()[1])
  rmsd = math.sqrt(flex.mean(displacements.dot(displacements)))

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
      reflections, experiments, verbosity=1, copy_experiments=True)

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
    # choose best orientation matrix, also don't use macrocycles) of refinement after indexing.
    # 2017 update: do accept multiple lattices per shot
    if self.params.refinement_protocol.n_macro_cycles > 1:
      raise Sorry("For stills, please set refinement_protocol.n_macro_cycles = 1")

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
          logger.info("Finish searching for more lattices: %i unindexed reflections remaining." %(
            min_reflections_for_indexing))
          break

      n_lattices_previous_cycle = len(experiments)

      # index multiple lattices per shot
      if len(experiments) == 0:
        experiments.extend(self.find_lattices())
        if len(experiments) == 0:
          raise Sorry("No suitable lattice could be found.")
      else:
        try:
          new = self.find_lattices()
          experiments.extend(new)
        except Sorry:
          logger.info("Indexing remaining reflections failed")

      self.index_reflections(experiments, self.reflections)

      if len(experiments) == n_lattices_previous_cycle:
        # no more lattices found
        break

      ### TODO verify things don't need to be re-indexed if a target is provided
      if False: #self.params.known_symmetry.space_group is not None:
        # now apply the space group symmetry only after the first indexing
        # need to make sure that the symmetrized orientation is similar to the P1 model
        target_space_group = self.target_symmetry_primitive.space_group()
        for i_cryst, cryst in enumerate(experiments.crystals()):
          if i_cryst >= n_lattices_previous_cycle:
            new_cryst, cb_op_to_primitive = self.apply_symmetry(
              cryst, target_space_group)
            if self.cb_op_primitive_inp is not None:
              new_cryst = new_cryst.change_basis(self.cb_op_primitive_inp)
              logger.info(new_cryst.get_space_group().info())
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

      # discard nearly overlapping lattices on the same shot
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

      self.indexed_reflections = (self.reflections['id'] > -1)
      if self.d_min is None:
        sel = self.reflections['id'] <= -1
      else:
        sel = flex.bool(len(self.reflections), False)
        lengths = 1/self.reflections['rlp'].norms()
        isel = (lengths >= self.d_min).iselection()
        sel.set_selected(isel, True)
        sel.set_selected(self.reflections['id'] > -1, False)
      self.unindexed_reflections = self.reflections.select(sel)

      reflections_for_refinement = self.reflections.select(
        self.indexed_reflections)

      if len(self.params.stills.isoforms) > 0:
        logger.info("")
        logger.info("#" * 80)
        logger.info("Starting refinement")
        logger.info("#" * 80)
        logger.info("")

        import copy
        isoform_experiments = ExperimentList()
        isoform_reflections = flex.reflection_table()
        # Note, changes to params after initial indexing. Cannot use tie to target when fixing the unit cell.
        self.all_params.refinement.reflections.outlier.algorithm = "null"
        self.all_params.refinement.parameterisation.crystal.fix = "cell"
        self.all_params.refinement.parameterisation.crystal.unit_cell.restraints.tie_to_target = []

        for expt_id, experiment in enumerate(experiments):
          reflections = reflections_for_refinement.select(reflections_for_refinement['id'] == expt_id)
          reflections['id'] = flex.int(len(reflections),0)
          refiners = []
          for isoform in self.params.stills.isoforms:
            iso_experiment = copy.deepcopy(experiment)
            crystal = iso_experiment.crystal
            if isoform.lookup_symbol != crystal.get_space_group().type().lookup_symbol():
              logger.info("Crystal isoform lookup_symbol %s does not match isoform %s lookup_symbol %s"%(crystal.get_space_group().type().lookup_symbol(), isoform.name, isoform.lookup_symbol))
              continue
            crystal.set_B(isoform.cell.fractionalization_matrix())

            logger.info("Refining isoform %s"%isoform.name)
            refiners.append(e_refine(params=self.all_params, experiments=ExperimentList([iso_experiment]), reflections=reflections, graph_verbose=False))

          if len(refiners) == 0:
            raise Sorry("No isoforms had a lookup symbol that matched")
          positional_rmsds = [math.sqrt(P.rmsds()[0] ** 2 + P.rmsds()[1] ** 2) for P in refiners]
          logger.info("Positional rmsds for all isoforms:" + str(positional_rmsds))
          minrmsd_mm = min(positional_rmsds)
          minindex = positional_rmsds.index(minrmsd_mm)
          logger.info("The smallest rmsd is %5.1f um from isoform %s" % (
            1000. * minrmsd_mm, self.params.stills.isoforms[minindex].name))
          if self.params.stills.isoforms[minindex].rmsd_target_mm is not None:
            logger.info("Asserting %f < %f"%(minrmsd_mm, self.params.stills.isoforms[minindex].rmsd_target_mm))
            assert minrmsd_mm < self.params.stills.isoforms[minindex].rmsd_target_mm
          logger.info("Acceptable rmsd for isoform %s." % (self.params.stills.isoforms[minindex].name))
          if len(self.params.stills.isoforms) == 2:
            logger.info("Rmsd gain over the other isoform %5.1f um." % (
            1000. * abs(positional_rmsds[0] - positional_rmsds[1])))
          R = refiners[minindex]
          # Now one last check to see if direct beam is out of bounds
          if self.params.stills.isoforms[minindex].beam_restraint is not None:
            from scitbx import matrix
            refined_beam = matrix.col(
              R.get_experiments()[0].detector[0].get_beam_centre_lab(experiments[0].beam.get_s0())[0:2])
            known_beam = matrix.col(self.params.stills.isoforms[minindex].beam_restraint)
            logger.info("Asserting difference in refined beam center and expected beam center %f < %f"%((refined_beam - known_beam).length(), self.params.stills.isoforms[
              minindex].rmsd_target_mm))
            assert (refined_beam - known_beam).length() < self.params.stills.isoforms[minindex].rmsd_target_mm
            # future--circle of confusion could be given as a separate length in mm instead of reusing rmsd_target

          experiment = R.get_experiments()[0]
          experiment.crystal.identified_isoform = self.params.stills.isoforms[minindex].name

          isoform_experiments.append(experiment)
          reflections['id'] = flex.int(len(reflections),expt_id)
          isoform_reflections.extend(reflections)
        experiments = isoform_experiments
        reflections_for_refinement = isoform_reflections

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
          logger.info("Refinement failed:")
          logger.info(s)
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

    # discard experiments with zero reflections after refinement
    id_set = set(self.refined_reflections['id'])
    if len(id_set) < len(self.refined_experiments):
      filtered_refined_reflections = flex.reflection_table()
      for i in xrange(len(self.refined_experiments)):
        if i not in id_set:
          del self.refined_experiments[i]
      for old, new in zip(sorted(id_set), range(len(id_set))):
        subset = self.refined_reflections.select(self.refined_reflections['id'] == old)
        subset['id'] = flex.int(len(subset), new)
        filtered_refined_reflections.extend(subset)
      self.refined_reflections = filtered_refined_reflections

    if len(self.refined_experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
           import show_rotation_matrix_differences
      show_rotation_matrix_differences(
        self.refined_experiments.crystals(), out=info_handle)

    logger.info("Final refined crystal models:")
    for i, crystal_model in enumerate(self.refined_experiments.crystals()):
      n_indexed = 0
      for i_expt in experiments.where(crystal=crystal_model):
        n_indexed += (self.reflections['id'] == i).count(True)
      logger.info("model %i (%i reflections):" %(i+1, n_indexed))
      logger.info(crystal_model)

    if 'xyzcal.mm' in self.refined_reflections: # won't be there if refine_all_candidates = False and no isoforms
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
    from dxtbx.model.experiment_list import Experiment, ExperimentList
    import copy

    logger.info('*' * 80)
    logger.info('Selecting the best orientation matrix')
    logger.info('*' * 80)

    from libtbx import group_args
    class candidate_info(group_args):
      pass
    candidates = []

    params = copy.deepcopy(self.all_params)

    n_cand = len(candidate_orientation_matrices)

    for icm,cm in enumerate(candidate_orientation_matrices):
      # Index reflections in P1
      sel = ((self.reflections['id'] == -1))
      refl = self.reflections.select(sel)
      experiments = self.experiment_list_for_crystal(cm)
      self.index_reflections(experiments, refl)
      indexed = refl.select(refl['id'] >= 0)
      indexed = indexed.select(indexed.get_flags(indexed.flags.indexed))

      # If target symmetry applied, try to apply it.  Then, apply the change of basis to the reflections
      # indexed in P1 to the target setting
      if self.params.known_symmetry.space_group is not None:
        target_space_group = self.target_symmetry_primitive.space_group()
        new_crystal, cb_op_to_primitive = self.apply_symmetry(cm, target_space_group)
        if new_crystal is None:
          print "Cannot convert to target symmetry, candidate %d/%d"%(icm, n_cand)
          continue
        new_crystal = new_crystal.change_basis(self.cb_op_primitive_inp)
        cm = candidate_orientation_matrices[icm] = new_crystal

        if not cb_op_to_primitive.is_identity_op():
          indexed['miller_index'] = cb_op_to_primitive.apply(indexed['miller_index'])
        if self.cb_op_primitive_inp is not None:
          indexed['miller_index'] = self.cb_op_primitive_inp.apply(indexed['miller_index'])

      if params.indexing.stills.refine_all_candidates:
        try:
          print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d/%d initial outlier identification"%(icm, n_cand)
          acceptance_flags = self.identify_outliers(params, experiments, indexed)
          #create a new "indexed" list with outliers thrown out:
          indexed = indexed.select(acceptance_flags)

          print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d/%d refinement before outlier rejection"%(icm, n_cand)
          R = e_refine(params = params, experiments=experiments, reflections=indexed, graph_verbose=False)
          ref_experiments = R.get_experiments()

          # try to improve the outcome with a second round of outlier rejection post-initial refinement:
          acceptance_flags = self.identify_outliers(params, ref_experiments, indexed)

          # insert a round of Nave-outlier rejection on top of the r.m.s.d. rejection
          nv0 = nave_parameters(params = params, experiments=ref_experiments, reflections=indexed, refinery=R, graph_verbose=False)
          crystal_model_nv0 = nv0()
          acceptance_flags_nv0 = nv0.nv_acceptance_flags
          indexed = indexed.select(acceptance_flags & acceptance_flags_nv0)

          print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d/%d after positional and delta-psi outlier rejection"%(icm, n_cand)
          R = e_refine(params = params, experiments=ref_experiments, reflections=indexed, graph_verbose=False)
          ref_experiments = R.get_experiments()

          nv = nave_parameters(params = params, experiments=ref_experiments, reflections=indexed, refinery=R, graph_verbose=False)
          crystal_model = nv()

          rmsd, _ = calc_2D_rmsd_and_displacements(R.predict_for_reflection_table(indexed))
        except Exception, e:
          print "Couldn't refine candiate %d/%d, %s"%(icm, n_cand, str(e))
        else:
          print "$$$ stills_indexer::choose_best_orientation_matrix, candidate %d/%d done"%(icm, n_cand)
          candidates.append(candidate_info(crystal = crystal_model,
                                           green_curve_area = nv.green_curve_area,
                                           ewald_proximal_volume = nv.ewald_proximal_volume(),
                                           n_indexed = len(indexed),
                                           rmsd = rmsd,
                                           indexed = indexed,
                                           experiments = ref_experiments))
      else:
        from dials.algorithms.refinement.prediction import ExperimentsPredictor
        ref_predictor = ExperimentsPredictor(experiments, force_stills=True,
                                             spherical_relp=params.refinement.parameterisation.spherical_relp_model)
        rmsd, _ = calc_2D_rmsd_and_displacements(ref_predictor(indexed))
        candidates.append(candidate_info(crystal = cm,
                                         n_indexed = len(indexed),
                                         rmsd = rmsd,
                                         indexed = indexed,
                                         experiments = experiments))
    if len(candidates) == 0:
      raise Sorry("No suitable indexing solution found")

    print "**** ALL CANDIDATES:"
    for i,XX in enumerate(candidates):
      print "\n****Candidate %d"%i,XX
      cc = XX.crystal
      if hasattr(cc, '_ML_half_mosaicity_deg'):
        print "  half mosaicity %5.2f deg."%(cc._ML_half_mosaicity_deg)
        print "  domain size %.0f Ang."%(cc._ML_domain_size_ang)
    print "\n**** BEST CANDIDATE:"

    results = flex.double([c.rmsd for c in candidates])
    best = candidates[flex.min_index(results)]
    print best

    if params.indexing.stills.refine_all_candidates:
      if best.rmsd > params.indexing.stills.rmsd_min_px:
        raise Sorry ("RMSD too high, %f" %rmsd)

      if best.ewald_proximal_volume > params.indexing.stills.ewald_proximal_volume_max:
        raise Sorry ("Ewald proximity volume too high, %f"%best.ewald_proximal_volume)

      if len(candidates) > 1:
        for i in xrange(len(candidates)):
          if i == flex.min_index(results):
            continue
          if best.ewald_proximal_volume > candidates[i].ewald_proximal_volume:
            print "Couldn't figure out which candidate is best; picked the one with the best RMSD."

    best.indexed['entering'] = flex.bool(best.n_indexed, False)

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
    acceptance_flags = self.identify_outliers(self.all_params, experiments, reflections)
    #create a new "reflections" list with outliers thrown out:
    reflections = reflections.select(acceptance_flags)

    R = e_refine(params = self.all_params, experiments=experiments, reflections=reflections, graph_verbose=False)
    ref_experiments = R.get_experiments()

    # try to improve the outcome with a second round of outlier rejection post-initial refinement:
    acceptance_flags = self.identify_outliers(self.all_params, ref_experiments, reflections)

    # insert a round of Nave-outlier rejection on top of the r.m.s.d. rejection
    nv0 = nave_parameters(params = self.all_params, experiments=ref_experiments, reflections=reflections, refinery=R, graph_verbose=False)
    crystal_model_nv0 = nv0()
    acceptance_flags_nv0 = nv0.nv_acceptance_flags
    reflections = reflections.select(acceptance_flags & acceptance_flags_nv0)

    R = e_refine(params = self.all_params, experiments=ref_experiments, reflections=reflections, graph_verbose=False)
    ref_experiments = R.get_experiments()

    nv = nave_parameters(params = self.all_params, experiments=ref_experiments, reflections=reflections, refinery=R, graph_verbose=False)
    nv()
    rmsd, _ = calc_2D_rmsd_and_displacements(R.predict_for_reflection_table(reflections))

    matches = R.get_matches()
    xyzcal_mm = flex.vec3_double(len(reflections))
    xyzcal_mm.set_selected(matches['iobs'], matches['xyzcal.mm'])
    reflections['xyzcal.mm'] = xyzcal_mm
    reflections.set_flags(matches['iobs'], reflections.flags.used_in_refinement)
    reflections['entering'] = flex.bool(len(reflections), False)
    return ref_experiments, reflections

""" Mixin class definitions that override the dials indexing class methods specific to stills """
class stills_indexer_known_orientation(indexer_known_orientation, stills_indexer):
  pass

class stills_indexer_real_space_grid_search(stills_indexer, indexer_real_space_grid_search):
  pass

class stills_indexer_fft3d(stills_indexer, indexer_fft3d):
  pass

class stills_indexer_fft1d(stills_indexer, indexer_fft1d):
  pass
