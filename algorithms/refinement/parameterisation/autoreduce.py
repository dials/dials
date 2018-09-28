from __future__ import absolute_import, division
import logging
logger = logging.getLogger(__name__)

from dials_refinement_helpers_ext import pgnmn_iter as pgnmn
from dials_refinement_helpers_ext import ucnmn_iter as ucnmn
from dials_refinement_helpers_ext import mnmn_iter as mnmn

from scitbx.array_family import flex

# PHIL
phil_str = '''
      min_nref_per_parameter = 5
        .help = "the smallest number of reflections per parameter for a"
                "model parameterisation below which the parameterisation will"
                "not be made in full, but the action described below will be"
                "triggered."
        .type = int(value_min=1)

      action = *fail fix remove
        .help = "action to take if there are too few reflections across the"
                "experiments related to a particular model parameterisation."
                "If fail, an exception will be raised and refinement will not"
                "proceed. If fix, refinement will continue but with the"
                "parameters relating to that model remaining fixed at their"
                "initial values. If remove, parameters relating to that model"
                "will be fixed, and in addition all reflections related to"
                "that parameterisation will be removed. This will therefore"
                "remove these reflections from other parameterisations of the"
                "global model too. For example, if a crystal model could not"
                "be parameterised it will be excised completely and not"
                "contribute to the joint refinement of the detector and beam."
                "In the fix mode, reflections emanating from that crystal will"
                "still form residuals and will contribute to detector and beam"
                "refinement."
        .type = choice

      detector_reduce = False
        .type = bool
        .help = "Special case designed for detector metrology refinement"
                "(particularly of the CSPAD). See detector_reduce_list for"
                "details."
        .expert_level = 1

      detector_reduce_list = Dist Tau2 Tau3
        .type = strings
        .help = "Partial names to match to detector parameters to try fixing."
                "If there are still not"
                "enough parameters for refinement after fixing these, then"
                "fail. This is to ensure that metrology refinement never"
                "completes if it is not able to refine some panels. The default"
                "is to try fixing the distance as well as Tau2 and Tau3"
                "rotations of detector panel, leaving the in-plane shifts and"
                "the rotation around the detector normal for refinement."
                "groups only."
        .expert_level = 1
'''

# Parameter auto reduction options
def model_nparam_minus_nref(options, p, reflections):
  cutoff = options.auto_reduction.min_nref_per_parameter * p.num_free()

  #Replaced Python code
  '''
  exp_ids = p.get_experiment_ids()
  # Do we have enough reflections to support this parameterisation?
  nparam = p.num_free()
  cutoff = options.auto_reduction.min_nref_per_parameter * nparam
  isel = flex.size_t()
  for exp_id in exp_ids:
    isel.extend((reflections['id'] == exp_id).iselection())
  nref = len(isel)
  return nref - cutoff
  '''
  return mnmn(reflections["id"],p.get_experiment_ids()).result - cutoff

def unit_cell_nparam_minus_nref(options, p, reflections):
  '''Special version of model_nparam_minus_nref for crystal unit cell
  parameterisations. In some cases certain parameters of a unit cell
  parameterisation may affect only some subset of the total number of
  reflections. For example, for an orthorhombic cell the g_param_0 parameter
  has no effect on predictions in the plane (0,k,l). Here, take the number
  of affected reflections for each parameter into account.'''

  F_dbdp=flex.mat3_double( p.get_ds_dp() )
  min_nref = options.auto_reduction.min_nref_per_parameter
  # if no free parameters, do as model_nparam_minus_nref
  if len(F_dbdp) == 0:
    exp_ids = p.get_experiment_ids()
    isel = flex.size_t()
    for exp_id in exp_ids:
      isel.extend((reflections['id'] == exp_id).iselection())
    return len(isel)

  #Replaced Python code
  '''
  exp_ids = p.get_experiment_ids()
  isel = flex.size_t()
  for exp_id in exp_ids:
    isel.extend((reflections['id'] == exp_id).iselection())
  ref = reflections.select(isel)
  h = ref['miller_index'].as_vec3_double()
  dB_dp = p.get_ds_dp()
  # if no free parameters, do as model_nparam_minus_nref
  if len(dB_dp) == 0: return len(isel)
  nref_each_param = []
  min_nref = options.auto_reduction.min_nref_per_parameter
  for der in dB_dp:
    der_mat = flex.mat3_double(len(h), der.elems)
    tst = (der_mat * h).norms()
    nref_each_param.append((tst > 0.0).count(True))

  return min([nref - min_nref for nref in nref_each_param])
  '''
  return ucnmn(reflections["id"], reflections["miller_index"], p.get_experiment_ids(), F_dbdp).result - min_nref

def panel_gp_nparam_minus_nref(options, p, pnl_ids, group, reflections, verbose=False):
  """
  :param p: ModelParameterisation; parameters in model
  :param pnl_ids: panel IDs
  :param group: group ID
  :panel reflections: flex table of reflections
  :panel verbose:
  :return: returns surplus {int}
  """
  exp_ids = p.get_experiment_ids() #Experiments parameterised by this ModelParameterisation
  # Do we have enough reflections to support this parameterisation?
  gp_params = [gp == group for gp in p.get_param_panel_groups()] #select the group ids for each param that matches the arg 'group'
  fixlist = p.get_fixed() # Get the fixed parameters; list says yes or no over all
  free_gp_params = [a and not b for a,b in zip(gp_params, fixlist)] #Free params is the total less the fixed
  nparam = free_gp_params.count(True)
  cutoff = options.auto_reduction.min_nref_per_parameter * nparam
  isel = flex.size_t()
  #Use Boost.Python extension module to replace below code
  surplus = pgnmn(reflections["id"], reflections["panel"], pnl_ids, exp_ids, cutoff).result

  #Replaced Python code
  '''
  for exp_id in exp_ids:
    sub_expID = (reflections['id'] == exp_id).iselection()
    sub_panels_expID = reflections['panel'].select(sub_expID)
    for pnl in pnl_ids:
      isel.extend(sub_expID.select(sub_panels_expID == pnl))
  nref = len(isel)
  surplus = nref - cutoff
  '''
  if surplus < 0 and verbose:
    logger.warning('{0} reflections on panels {1} with a cutoff of {2}'.format(nref, pnl_ids, cutoff))
  return surplus

def weak_parameterisation_search(options, beam_params, xl_ori_params, xl_uc_params,
  det_params, gon_params, reflections):
  weak = None
  nref_deficit = 0
  panels = None
  pnl_gp = None
  name = None
  for i, p in enumerate(beam_params):
    net_nref = model_nparam_minus_nref(options, p, reflections)
    if net_nref < nref_deficit:
      nref_deficit = net_nref
      weak = p
      name = 'Beam{0}'.format(i + 1)
  for i, p in enumerate(xl_ori_params):
    net_nref = model_nparam_minus_nref(options, p, reflections)
    if net_nref < nref_deficit:
      nref_deficit = net_nref
      weak = p
      name = 'Crystal{0} orientation'.format(i + 1)
  for i, p in enumerate(xl_uc_params):
    net_nref = unit_cell_nparam_minus_nref(options, p, reflections)
    if net_nref < nref_deficit:
      nref_deficit = net_nref
      weak = p
      name = 'Crystal{0} unit cell'.format(i + 1)
  for i, p in enumerate(det_params):
    try:
      pnl_groups = p.get_panel_ids_by_group()
      for igp, gp in enumerate(pnl_groups):
        net_nref = panel_gp_nparam_minus_nref(options, p, gp, igp, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          panels = gp
          pnl_gp = igp
          name = 'Detector{0}PanelGroup{1}'.format(i + 1, pnl_gp + 1)
    except Exception:
      net_nref = model_nparam_minus_nref(options, p, reflections)
      if net_nref < nref_deficit:
        nref_deficit = net_nref
        weak = p
        panels = None
        pnl_gp = None
        name = 'Detector{0}'.format(i + 1)
  for i, p in enumerate(gon_params):
    net_nref = model_nparam_minus_nref(options, p, reflections)
    if net_nref < nref_deficit:
      nref_deficit = net_nref
      weak = p
      name = 'Goniometer{0}'.format(i + 1)
  return {'parameterisation':weak,
          'panels':panels,
          'panel_group_id':pnl_gp,
          'name':name}

def autoreduce(options, det_params, beam_params, xl_ori_params, xl_uc_params,
    gon_params, reflections):

  # In the scan-varying case we can't calculate dB_dp before composing the
  # model, so revert to the original function
  if options.scan_varying:
    _unit_cell_nparam_minus_nref = model_nparam_minus_nref
  else:
    _unit_cell_nparam_minus_nref = unit_cell_nparam_minus_nref

  # As a special case for detector metrology, try reducing the number of
  # detector parameters if there are too few for some panel group. If this is
  # unsuccessful, fail outright.
  if options.auto_reduction.detector_reduce:
    reduce_list = options.auto_reduction.detector_reduce_list
    for i, dp in enumerate(det_params):
      to_fix = flex.bool(dp.get_fixed())
      try: # test for hierarchical detector parameterisation
        pnl_groups = dp.get_panel_ids_by_group()
        for igp, gp in enumerate(pnl_groups):
          surplus = panel_gp_nparam_minus_nref(options, dp, gp, igp, reflections, verbose=True)
          if surplus < 0:
            msg = ('Require {0} more reflections to parameterise Detector{1} '
                   'panel group {2}').format(-1*surplus, i + 1, igp + 1)
            logger.warning(msg + '\nAttempting reduction of non-essential parameters')
            names = cls._filter_parameter_names(dp)
            prefix = 'Group{0}'.format(igp + 1)
            reduce_this_group = [prefix + e for e in reduce_list]
            to_fix |= flex.bool(string_sel(reduce_this_group, names))
            # try again, and fail if still unsuccessful
            surplus = panel_gp_nparam_minus_nref(options, dp, gp, igp, reflections, verbose=True)
            if surplus < 0:
              msg = msg.format(-1*surplus, i + 1, igp + 1)
              raise Sorry(msg + '\nFailing.')
      except AttributeError:
        if model_nparam_minus_nref(options, dp, reflections) < 0:
          mdl = 'Detector{0}'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)
      dp.set_fixed(to_fix)

  if options.auto_reduction.action == 'fail':
    failmsg = 'Too few reflections to parameterise {0}'
    failmsg += '\nTry modifying refinement.parameterisation.auto_reduction options'
    for i, bp in enumerate(beam_params):
      if model_nparam_minus_nref(options, bp, reflections) < 0:
        mdl = 'Beam{0}'.format(i + 1)
        msg = failmsg.format(mdl)
        raise Sorry(msg)

    for i, xlo in enumerate(xl_ori_params):
      if model_nparam_minus_nref(options, xlo, reflections) < 0:
        mdl = 'Crystal{0} orientation'.format(i + 1)
        msg = failmsg.format(mdl)
        raise Sorry(msg)

    for i, xluc in enumerate(xl_uc_params):
      if _unit_cell_nparam_minus_nref(options, xluc, reflections) < 0:
        mdl = 'Crystal{0} unit cell'.format(i + 1)
        msg = failmsg.format(mdl)
        raise Sorry(msg)

    for i, dp in enumerate(det_params):
      try: # test for hierarchical detector parameterisation
        pnl_groups = dp.get_panel_ids_by_group()
        for igp, gp in enumerate(pnl_groups):
          if panel_gp_nparam_minus_nref(options, dp, gp, igp, reflections) < 0:
            msg = 'Too few reflections to parameterise Detector{0} panel group {1}'
            msg = msg.format(i + 1, igp + 1)
            msg += '\nTry modifying refinement.parameterisation.auto_reduction options'
            raise Sorry(msg)
      except AttributeError:
        if model_nparam_minus_nref(options, dp, reflections) < 0:
          mdl = 'Detector{0}'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

    for i, gonp in enumerate(gon_params):
      if model_nparam_minus_nref(options, gonp, reflections) < 0:
        mdl = 'Goniometer{0}'.format(i + 1)
        msg = failmsg.format(mdl)
        raise Sorry(msg)

  elif options.auto_reduction.action == 'fix':
    warnmsg = 'Too few reflections to parameterise {0}'
    tmp = []
    for i, bp in enumerate(beam_params):
      if model_nparam_minus_nref(options, bp, reflections) >= 0:
        tmp.append(bp)
      else:
        mdl = 'Beam{0}'.format(i + 1)
        msg = warnmsg.format(mdl)
        logger.warning(msg)
    beam_params = tmp

    tmp = []
    for i, xlo in enumerate(xl_ori_params):
      if model_nparam_minus_nref(options, xlo, reflections) >= 0:
        tmp.append(xlo)
      else:
        mdl = 'Crystal{0} orientation'.format(i + 1)
        msg = warnmsg.format(mdl)
        logger.warning(msg)
    xl_ori_params = tmp

    tmp = []
    for i, xluc in enumerate(xl_uc_params):
      if _unit_cell_nparam_minus_nref(options, xluc, reflections) >= 0:
        tmp.append(xluc)
      else:
        mdl = 'Crystal{0} unit cell'.format(i + 1)
        msg = warnmsg.format(mdl)
        logger.warning(msg)
    xl_uc_params = tmp

    tmp = []
    for i, dp in enumerate(det_params):
      fixlist = dp.get_fixed()
      try: # test for hierarchical detector parameterisation
        pnl_groups = dp.get_panel_ids_by_group()
        for igp, gp in enumerate(pnl_groups):
          if panel_gp_nparam_minus_nref(options, dp, gp, igp, reflections) < 0:
            msg = 'Too few reflections to parameterise Detector{0}PanelGroup{1}'
            msg = msg.format(i + 1, igp + 1)
            logger.warning(msg)
            gp_params = [gp == igp for gp in dp.get_param_panel_groups()]
            for j, val in enumerate(gp_params):
              if val: fixlist[j] = True
        dp.set_fixed(fixlist)
        if dp.num_free() > 0:
          tmp.append(dp)
        else:
          msg = 'No parameters remain free for Detector{0}'.format(i + 1)
          logger.warning(msg)
      except AttributeError:
        if model_nparam_minus_nref(options, dp, reflections) >= 0:
          tmp.append(dp)
        else:
          mdl = 'Detector{0}'.format(i + 1)
          msg = warnmsg.format(mdl)
          logger.warning(msg)
    det_params = tmp

    tmp = []
    for i, gonp in enumerate(gon_params):
      if model_nparam_minus_nref(options, gonp, reflections) >= 0:
        tmp.append(gonp)
      else:
        mdl = 'Goniometer{0}'.format(i + 1)
        msg = warnmsg.format(mdl)
        logger.warning(msg)
    gon_params = tmp

  elif options.auto_reduction.action == 'remove':
    # if there is only one experiment, it should be multi-panel for remove to make sense
    if len(experiments) == 1:
      if not det_params[-1].is_multi_state():
        raise Sorry("For single experiment, single panel refinement "
          "auto_reduction.action=remove cannot be used as it could only "
          "remove all reflections from refinement")
    warnmsg = 'Too few reflections to parameterise {0}'
    warnmsg += '\nAssociated reflections will be removed from the Reflection Manager'
    while True:
      dat = weak_parameterisation_search(options, beam_params, xl_ori_params,
          xl_uc_params, det_params, gon_params, reflections)
      if dat['parameterisation'] is None: break
      exp_ids = dat['parameterisation'].get_experiment_ids()
      if dat['panels'] is not None:
        msg = warnmsg.format(dat['name'])
        fixlist = dat['parameterisation'].get_fixed()
        pnl_gps = dat['parameterisation'].get_param_panel_groups()
        for i, gp in enumerate(pnl_gps):
          if gp == dat['panel_group_id']: fixlist[i] = True
        dat['parameterisation'].set_fixed(fixlist)
        # identify observations on this panel group from associated experiments
        obs = reflection_manager.get_obs()
        isel=flex.size_t()
        for exp_id in exp_ids:
          subsel = (obs['id'] == exp_id).iselection()
          panels_this_exp = obs['panel'].select(subsel)
          for pnl in dat['panels']:
            isel.extend(subsel.select(panels_this_exp == pnl))
      else:
        msg = warnmsg.format(dat['name'])
        fixlist = [True] * dat['parameterisation'].num_total()
        dat['parameterisation'].set_fixed(fixlist)
        # identify observations from the associated experiments
        obs = reflection_manager.get_obs()
        isel=flex.size_t()
        for exp_id in exp_ids:
          isel.extend((obs['id'] == exp_id).iselection())
      # Now remove the selected reflections
      sel = flex.bool(len(obs), True)
      sel.set_selected(isel, False)
      reflection_manager.filter_obs(sel)
      reflections = reflection_manager.get_matches()
      logger.warning(msg)

    # Strip out parameterisations with zero free parameters
    beam_params = [p for p in beam_params if p.num_free() > 0]
    xl_ori_params = [p for p in xl_ori_params if p.num_free() > 0]
    xl_uc_params = [p for p in xl_uc_params if p.num_free() > 0]
    det_params = [p for p in det_params if p.num_free() > 0]
    gon_params = [p for p in gon_params if p.num_free() > 0]

  return det_params, beam_params, xl_ori_params, xl_uc_params, gon_params
