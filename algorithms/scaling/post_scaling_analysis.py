import logging
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.scaling_library import create_Ih_table

logger = logging.getLogger('dials')

def exclude_on_batch_rmerge(reflections, experiments, rmerge_cutoff=2.0):
  """Analyse R-merge per batch and set the user_excluded_in_scaling flag if
    R-merge is above rmerge_cutoff."""

  batchno = 0
  target_Ih_Table = create_Ih_table(experiments, reflections)
  new_reflections = []

  for i, (refl, exp) in enumerate(zip(reflections, experiments)):
    outlier_sel = refl.get_flags(refl.flags.bad_for_scaling, all=False)
    if outlier_sel.count(True) != refl.size():
      Ihtable = IhTable([(refl, None)], exp.crystal.get_space_group())
      data = Ihtable.blocked_data_list[0]
      data.set_Ih_values_to_target(target_Ih_Table.blocked_data_list[0])
      refl['Ih_values'] = flex.double(refl.size(), 0.0)
      refl['Ih_values'].set_selected(data.nonzero_weights, data.Ih_values)
      bad_sf_sel = refl['inverse_scale_factor'] <= 0.0
      sel3 = outlier_sel | bad_sf_sel
      if not exp.scan:
        # Assume no scan, therefore all data is in one batch
        data = refl.select(~sel3)
        if data:
          Rmerge = flex.sum(flex.abs((data['intensity']/data['inverse_scale_factor'])
            - data['Ih_values']))/ flex.sum((data['intensity']/data['inverse_scale_factor']))
        else:
          Rmerge = 0.0
        if Rmerge > rmerge_cutoff or Rmerge <= 0.0:
          refl.set_flags(flex.bool(refl.size(), True), refl.flags.user_excluded_in_scaling)
          logger.info(
            'Batch no %s has an Rmerge of %s, this data will be excluded for calculating '
            'merging statistics, further scaling and mtz output if applicable'
            % (batchno, Rmerge))
        batchno += 1
      else:
        z = refl['xyzobs.px.value'].parts()[2]
        refl_to_delete_sel = flex.bool(refl.size(), False)
        refl_to_delete = False
        i0, i1 = exp.scan.get_image_range()
        image_vals = range(i0, i1+1)
        for i in image_vals:
          if i == i0:
            sel1 = float(i - 1) - 0.3 < z #add some tolerance - necessary?
            sel2 = float(i) > z
          elif i == i1:
            sel1 = float(i - 1) < z
            sel2 = float(i) + 0.3 > z  #add some tolerance - necessary?
          else:
            sel1 = float(i - 1) < z
            sel2 = float(i) > z
          sel = sel1 & sel2
          sel4 = sel & ~sel3
          data = refl.select(sel4)
          if data:
            Rmerge = flex.sum(flex.abs((data['intensity']/data['inverse_scale_factor'])
              - data['Ih_values']))/ flex.sum((data['intensity']/data['inverse_scale_factor']))
          else:
            Rmerge = 0.0
          if Rmerge > rmerge_cutoff or Rmerge <= 0.0:
            refl_to_delete = True
            refl_to_delete_sel.set_selected(sel4.iselection(), True)
            logger.info(
              'Batch no %s has an Rmerge of %s, this data will be excluded for calculating '
              'merging statistics, further scaling and mtz output if applicable'
              % (batchno, Rmerge))
          batchno += 1
        if refl_to_delete:
          refl.set_flags(refl_to_delete_sel, refl.flags.user_excluded_in_scaling)
          logger.info("%s reflections excluded from dataset with experiment identifier %s" % (
            refl_to_delete_sel.count(True), exp.identifier))
    new_reflections.append(refl)

  return new_reflections

def exclude_on_image_scale(reflections, experiments, scale_cutoff=0.2):
  """Analyse the scale per image and set the user_excluded_in_scaling flag if
    the scale is below scale_cutoff."""
  new_reflections = []
  for i, (reflection, experiment) in enumerate(zip(reflections, experiments)):
    if experiment.scaling_model.id_ == 'KB':
      s = experiment.scaling_model.components['scale'].parameters[0]
      if s < scale_cutoff:
        logger.info('The presence of a bad image was detected for dataset %s'
          'Excluding %s reflections'
          % (experiment.identifier, reflection.size()))
        reflection.set_flags(flex.bool(reflection.size(), True),
          reflection.flags.user_excluded_in_scaling)
    elif experiment.scaling_model.id_ == 'physical':
      scale_values = flex.double([])
      config = experiment.scaling_model.configdict
      i0, i1 = experiment.scan.get_image_range()
      image_vals = range(i0, i1+1)
      scale_comp = experiment.scaling_model.components['scale']
      for i in image_vals:
        n = (i - 0.5) * config['s_norm_fac']
        s = scale_comp._smoother.value_weight(n, scale_comp.value)[0]
        scale_values.append(s)
      sel = scale_values < scale_cutoff
      bad_images = flex.int(image_vals).select(sel)
      if bad_images:
        logger.info('The presence of potentially bad images was detected in'
          'dataset %s. List of potential bad image numbers: \n%s' % (
            experiment.identifier, list(bad_images)))
        refl_to_delete_sel = flex.bool(reflection.size(), False)
        z = reflection['xyzobs.px.value'].parts()[2]
        for image in bad_images:
          sel1 = float(image - 1) < z
          sel2 = float(image) > z
          sel = sel1 & sel2
          logger.info('Excluding %s reflections from image number %s' % (
            sel.count(True), image))
          refl_to_delete_sel.set_selected(sel.iselection(), True)
        reflection.set_flags(refl_to_delete_sel,
          reflection.flags.user_excluded_in_scaling)
    else:
      logger.info(
        'Exclude on image scale not implemented for scaling model type %s'
        'Skipping dataset %s' % (experiment.scaling_model.id_,
        experiment.identifier))
    new_reflections.append(reflection)

  return new_reflections
