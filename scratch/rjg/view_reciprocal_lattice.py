from __future__ import division
from gltbx import wx_viewer
import wx
from gltbx.gl import *
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex

class MyGLWindow(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(MyGLWindow, self).__init__(*args, **kwds)
    self.points = flex.vec3_double()
    self.flag_show_minimum_covering_sphere = False
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)

  def set_points(self, points):
    self.points = points
    self.points_display_list = None
    self.draw_points()
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)
    self.fit_into_viewport()


class MyApp(wx_viewer.App):

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = MyGLWindow(self.frame, size=(600,600), orthographic=False)
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

  def load_models(self, imageset, reflections):
    from dials.algorithms.indexing import indexer
    reflections = indexer.indexer_base.map_spots_pixel_to_mm_rad(
      reflections, imageset.get_detector(), imageset.get_scan())
    indexer.indexer_base.map_centroids_to_reciprocal_space(
      reflections, imageset.get_detector(), imageset.get_beam(),
      imageset.get_goniometer())
    points = reflections['rlp'] * 100
    self.view_objects.set_points(points)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_reflections

  parser = OptionParser(
    phil=None,
    read_datablocks=True,
    read_reflections=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)[0]
  assert len(datablocks) == 1
  imageset = datablocks[0].extract_imagesets()[0]

  a = MyApp(title="Reciprocal space viewer")
  a.load_models(imageset, reflections)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
