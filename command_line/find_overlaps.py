from __future__ import division
#!/usr/bin/env python
#
# find_overlaps.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

import math
import iotbx.phil
from libtbx.phil import command_line
from dials.array_family import flex
from dxtbx.model.experiment.experiment_list import ExperimentList


master_phil_scope = iotbx.phil.parse("""
sigma_divergence = None
  .type = float
sigma_mosaicity = None
  .type = float
n_sigma = 3
  .type = float(value_min=0)
d_min = None
  .type = float(value_min=0.0)
max_overlap_fraction = 0.0
  .type = float(value_min=0.0)
max_overlap_pixels = 0
  .type = int(value_min=0)
nproc = 1
  .type = int(value_min=1)
save_overlaps = True
  .type = bool
""")


deg_to_rad = math.pi / 180


def run(args):
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
    args=args, custom_processor="collect_remaining")
  working_phil.show()
  params = working_phil.extract()

  from dials.util.command_line import Importer
  from libtbx import easy_pickle

  json_files = [arg for arg in args if arg.endswith('.json')]
  args = [arg for arg in args if arg not in json_files]
  experiments = ExperimentList()
  for json_file in json_files:
    importer = Importer([json_file], check_format=False)
    assert importer.experiments is not None
    assert len(importer.experiments) == 1
    experiments.extend(importer.experiments)
  reflections = []
  for arg in args:
    reflections.append(easy_pickle.load(arg))
  assert len(reflections) == len(experiments)

  reflection_table = flex.reflection_table()
  for i_lattice, ref_table in enumerate(reflections):
    ref_table['id'] = flex.int(len(ref_table), i_lattice)
    reflection_table.extend(ref_table)

  for i_lattice, expt in enumerate(experiments):
    if params.sigma_divergence is not None:
      expt.beam.set_sigma_divergence(params.sigma_divergence, deg=True)
    if params.sigma_mosaicity is not None:
      expt.crystal.set_mosaicity(params.sigma_mosaicity, deg=True)
    #print expt.beam.get_sigma_divergence(deg=True)
    #print expt.crystal.get_mosaicity(deg=True)

  overlaps = find_overlaps(
    experiments, reflection_table, n_sigma=params.n_sigma, nproc=params.nproc,
    d_min=params.d_min,
    max_overlap_fraction=params.max_overlap_fraction,
    max_overlap_pixels=params.max_overlap_pixels)
  if params.save_overlaps:
    overlaps.overlapping_reflections.as_pickle('overlaps.pickle')
  if 0:
    overlaps.plot_histograms()

  return overlaps


class find_overlaps(object):
  def __init__(self, experiments, reflections, n_sigma=3, nproc=1, d_min=None,
               max_overlap_fraction=0.0, max_overlap_pixels=0):
    from dials.algorithms import shoebox
    from dials.algorithms.shoebox import MaskCode

    reflection_table = self._prepare_reflections(experiments, reflections)
    if d_min is not None:
      S = reflection_table['s1'] - experiments[0].beam.get_s0()
      d_spacings = 1/S.norms()
      print (d_spacings > d_min).count(True), (d_spacings > d_min).count(False)
      reflection_table = reflection_table.select(d_spacings > d_min)
    reflection_table = self._predict_for_reflection_table(
      experiments, reflection_table)
    # XXX not sure why zero length s1 vectors come out of the prediction code
    reflection_table = reflection_table.select(
      reflection_table['s1'].norms() > 0)

    sel = flex.bool(len(reflection_table), False)
    frame_numbers = reflection_table['xyzcal.px'].parts()[2]
    for i_lattice, expt in enumerate(experiments):
      sel.set_selected(
        (reflection_table['id'] == i_lattice) &
        (frame_numbers >= expt.scan.get_array_range()[0]) &
        (frame_numbers < expt.scan.get_array_range()[1]), True)
    reflection_table = reflection_table.select(sel)

    reflection_table = self._compute_shoebox_mask(
      experiments, reflection_table, n_sigma=n_sigma)

    # Find overlapping reflections
    overlaps = shoebox.find_overlapping(reflection_table['bbox'])

    fg_code = MaskCode.Valid | MaskCode.Foreground
    # Loop through all edges
    overlapping = set()
    overlapping_bbox = set()
    self._n_overlapping_pixels = flex.size_t()
    self._fraction_overlapping_pixels = flex.double()
    #print len(list(overlaps.edges()))

    def func(args):
      overlapping = set()
      overlapping_bbox = set()
      n_overlapping_pixels = flex.size_t()
      fraction_overlapping_pixels = flex.double()
      edges = args[0]
      for e in edges:
        v1, v2 = overlaps[e]
        o1 = reflection_table[v1]
        o2 = reflection_table[v2]
        coords1 = o1['shoebox'].coords().select(
          (o1['shoebox'].mask == fg_code).as_1d())
        coords2 = o2['shoebox'].coords().select(
          (o2['shoebox'].mask == fg_code).as_1d())
        intersection = set(coords1).intersection(set(coords2))
        if len(intersection):
          n_overlapping_pixels.append(len(intersection))
          fraction_overlapping_pixels.append(len(intersection)/(
          len(coords1)+len(coords2)-len(intersection)))
          #print "Overlapping pixels: %i/%i (%.1f%%)" %(
            #len(intersection), (len(coords1) + len(coords2)),
            #len(intersection)/(len(coords1) + len(coords2)) * 100)
          if (n_overlapping_pixels[-1] > max_overlap_pixels or
              fraction_overlapping_pixels[-1] > max_overlap_fraction):
            overlapping.add(v1)
            overlapping.add(v2)
        overlapping_bbox.add(v1)
        overlapping_bbox.add(v2)
      return (overlapping, overlapping_bbox,
              n_overlapping_pixels, fraction_overlapping_pixels)

    edges = list(overlaps.edges())
    args = []
    n_per_proc = int(math.ceil(len(edges) / nproc))
    for i in range(nproc):
      args.append((edges[i*n_per_proc:min((i+1)*n_per_proc, len(edges))],
                   ))
    from libtbx import easy_mp
    results = easy_mp.parallel_map(func, args, processes=nproc)
    for r in results:
      overlapping.update(r[0])
      overlapping_bbox.update(r[1])
      self._n_overlapping_pixels.extend(r[2])
      self._fraction_overlapping_pixels.extend(r[3])

    print "Number of overlaps: %i/%i (%.1f%%)" %(
      len(overlapping), len(reflection_table),
      100*len(overlapping)/len(reflection_table))

    print "Number of overlapping bboxes: %i/%i (%.1f%%)" %(
      len(overlapping_bbox), len(reflection_table),
      100*len(overlapping_bbox)/len(reflection_table))

    miller_indices = reflection_table['miller_index']
    overlapping_isel = flex.size_t(list(overlapping))
    overlapping_sel = flex.bool(len(reflection_table), False)
    overlapping_sel.set_selected(overlapping_isel, True)

    self.overlap_selection = overlapping_sel
    self.reflections = reflection_table
    self.overlapping_reflections = self.reflections.select(overlapping_sel)

  def plot_histograms(self):
    from matplotlib import pyplot
    hist = flex.histogram(self._n_overlapping_pixels.as_double(), n_slots=50)
    #pyplot.bar(hist.slot_centers(), hist.slots(), width=hist.slot_width())
    pyplot.plot(hist.slot_centers(), hist.slots())
    pyplot.yscale('log')
    pyplot.xlabel('Number of overlapping pixels')
    pyplot.ylabel('Frequency')
    #pyplot.show()

    hist = flex.histogram(self._fraction_overlapping_pixels, n_slots=50)
    #pyplot.bar(hist.slot_centers(), hist.slots(), width=hist.slot_width())
    pyplot.plot(hist.slot_centers(), hist.slots())
    pyplot.yscale('log')
    pyplot.xlabel('Fraction of overlapping pixels')
    pyplot.ylabel('Frequency')
    #pyplot.show()

  def _prepare_reflections(self, experiments, reflections):
    reflection_table = flex.reflection_table()
    for i_lattice, expt in enumerate(experiments):
      ref_table = reflections.select(reflections['id'] == i_lattice)

      if not ref_table.has_key('xyzobs.mm.value'):

        xyzobs_px_value_orig = ref_table['xyzobs.px.value'].deep_copy()

        # map centroids pixels to mm
        from dials.algorithms.indexing.indexer import indexer_base
        ref_table['xyzobs.mm.value'] = flex.vec3_double(len(ref_table))
        ref_table['xyzobs.mm.variance'] = flex.vec3_double(len(ref_table))
        # some centroids are unobserved (0.0,0.0,0.0) so use xyzcal.px values instead
        ref_table['xyzobs.px.value'] = ref_table['xyzcal.px']
        ref_table['xyzobs.px.variance'] = flex.vec3_double(len(ref_table), (1,1,1))
        ref_table = indexer_base.map_spots_pixel_to_mm_rad(
          ref_table, expt.detector, expt.scan)

        # reset these to the original values
        ref_table['xyzobs.px.value'] = xyzobs_px_value_orig

      # compute s1 vectors
      if not ref_table.has_key('s1'):
        ref_table['s1'] = flex.vec3_double(len(ref_table))
        panel_numbers = flex.size_t(ref_table['panel'])
        reciprocal_space_points = flex.vec3_double()
        for i_panel in range(len(expt.detector)):
          sel = (panel_numbers == i_panel)
          isel = sel.iselection()
          spots_panel = ref_table.select(panel_numbers == i_panel)
          x, y, rot_angle = spots_panel['xyzobs.mm.value'].parts()
          s1 = expt.detector[i_panel].get_lab_coord(flex.vec2_double(x,y))
          s1 = s1/s1.norms() * (1/expt.beam.get_wavelength())
          ref_table['s1'].set_selected(isel, s1)

      reflection_table.extend(ref_table)

    return reflection_table

  def _predict_for_reflection_table(self, experiments, reflections):
    assert reflections.has_key('s1')
    reflection_table = flex.reflection_table()

    # calculate entering flags
    from dials.algorithms.refinement.reflection_manager import calculate_entering_flags
    entering = calculate_entering_flags(reflections, experiments)
    reflections['entering'] = entering

    for i_lattice, expt in enumerate(experiments):
      ref_table = reflections.select(reflections['id'] == i_lattice)
      miller_indices = ref_table['miller_index']

      # predict reflections
      from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
      predict = ScanStaticReflectionPredictor(expt)
      ub_array = flex.mat3_double(miller_indices.size(), expt.crystal.get_A())
      predict.for_reflection_table(ref_table, expt.crystal.get_A())

      # re-calculate phi ourselves since the predicted phi is in the range 0 -> 2pi
      for i, refl in enumerate(ref_table):
        phi_deg = expt.scan.get_angle_from_array_index(refl['xyzcal.px'][2])
        phi_rad = deg_to_rad * phi_deg
        pred = refl['xyzcal.mm']
        ref_table['xyzcal.mm'][i] = (pred[0], pred[1], phi_rad)

      reflection_table.extend(ref_table)

    return reflection_table

  def _compute_shoebox_mask(self, experiments, reflections, n_sigma=3):
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    from dials.algorithms.profile_model.profile_model import ProfileModel
    assert reflections.has_key('xyzcal.mm')
    assert reflections.has_key('id')
    reflection_table = reflections.copy()

    # Create the profile model
    profile_model = ProfileModelList()
    for i_lattice, expt in enumerate(experiments):
      profile_model.append(ProfileModel(
        n_sigma,
        expt.beam.get_sigma_divergence(deg=False),
        expt.crystal.get_mosaicity(deg=False)))

    # Compute the bounding boxes
    reflection_table.compute_bbox(experiments, profile_model, sigma_b_multiplier=1.0)

    # Allocate the shoeboxes
    reflection_table['shoebox'] = flex.shoebox(
      reflection_table['panel'],
      reflection_table['bbox'])
    reflection_table['shoebox'].allocate_with_value(MaskCode.Valid)

    # Create the function to mask the shoebox profiles
    reflection_table.compute_mask(experiments, profile_model)

    # Return reflection table
    return reflection_table


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
