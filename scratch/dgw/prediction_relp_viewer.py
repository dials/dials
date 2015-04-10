# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
from dials.command_line.reciprocal_lattice_viewer import *

class PredRelpViewer(ReciprocalLatticeViewer):

  # modified from dials.algorithms.indexing.indexer to map predictions
  # from pixels to mm, rather than observations
  @staticmethod
  def map_spots_pixel_to_mm_rad(spots, detector, scan):
    from dials.algorithms.centroid import centroid_px_to_mm_panel
    ## ideally don't copy, but have separate spot attributes for mm and pixel
    import copy
    spots = copy.deepcopy(spots)
    spots['xyzcal.mm'] = flex.vec3_double(len(spots))
    #spots['xyzobs.mm.variance'] = flex.vec3_double(len(spots))
    dummy_variance = flex.vec3_double(len(spots))
    panel_numbers = flex.size_t(spots['panel'])
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      isel = sel.iselection()
      spots_panel = spots.select(panel_numbers == i_panel)
      centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
        detector[i_panel], scan,
        spots_panel['xyzcal.px'],
        dummy_variance,
        flex.vec3_double(len(spots_panel), (1,1,1)))
      spots['xyzcal.mm'].set_selected(sel, centroid_position)
      #spots['xyzobs.mm.variance'].set_selected(sel, centroid_variance)
    return spots

  @staticmethod
  def map_centroids_to_reciprocal_space(spots_mm, detector,
    beam, goniometer):

    if 'imageset_id' not in spots_mm:
      spots_mm['imageset_id'] = spots_mm['id']
    if 's1' not in spots_mm: spots_mm['s1'] = flex.vec3_double(len(spots_mm))
    spots_mm['rlp'] = flex.vec3_double(len(spots_mm))
    panel_numbers = flex.size_t(spots_mm['panel'])
    for i_panel in range(len(detector)):
      sel = (panel_numbers == i_panel)
      spots_panel = spots_mm.select(panel_numbers == i_panel)
      x, y, rot_angle = spots_panel['xyzcal.mm'].parts()
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

  # override this from the base class to map predictions rather than
  # observations to reciprocal space
  def map_points_to_reciprocal_space(self):
    goniometer = copy.deepcopy(self.goniometer)
    if goniometer is not None and self.settings.reverse_phi:
      goniometer.set_rotation_axis([-i for i in goniometer.get_rotation_axis()])


    # set parameters to latest values
    self.dp.set_param_vals(self.detector_parameters)
    self.bp.set_param_vals(self.beam_parameters)

    reflections = self.map_spots_pixel_to_mm_rad(
      self.reflections, self.detector, self.scan)
    self.map_centroids_to_reciprocal_space(reflections, self.detector,
      self.beam, goniometer)

    d_spacings = 1/reflections['rlp'].norms()
    if self.settings.d_min is not None:
      reflections = reflections.select(d_spacings > self.settings.d_min)
    else:
      self.settings.d_min = flex.min(d_spacings)
      self.settings_panel.d_min_ctrl.SetValue(self.settings.d_min)
    points = reflections['rlp'] * 100
    self.viewer.set_points(points)

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections

  parser = OptionParser(
    phil=master_phil,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)[0]

  if len(datablocks) == 0:
    if len(experiments) > 0:
      imagesets = experiments.imagesets()
    else:
      parser.print_help()
      return
  elif len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()

  if len(imagesets) > 1:
    raise Sorry("Only one ImageSet can be processed at a time")
  imageset = imagesets[0]

  import wxtbx.app
  a = wxtbx.app.CCTBXApp(0)
  a.settings = params
  f = PredRelpViewer(
    None, -1, "Prediction reciprocal lattice viewer", size=(1024,768))
  f.load_models(imageset, reflections)
  f.Show()
  a.SetTopWindow(f)
  #a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
