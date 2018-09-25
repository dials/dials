# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# LIBTBX_SET_DISPATCHER_NAME dev.dials.diffraction_geometry_viewer

from __future__ import absolute_import, division, print_function

from dials.command_line.reciprocal_lattice_viewer import *

class DGVWindow(RLVWindow):

  def draw_ewald_sphere(self):

    if self.beam_vector is None: return
    from scitbx import matrix
    s0 = matrix.col(self.beam_vector)
    ewald_sphere = gltbx.util.WireSphere
    x = (-100*s0).elems
    r = 100*s0.length()
    grid = 200
    glPushMatrix()
    glTranslated(*(x))
    ewald_sphere(radius=r, slices=grid, stacks=grid)
    glPopMatrix()

  def DrawGL(self):
    super(DGVWindow, self).DrawGL()
    self.draw_ewald_sphere()

class DiffractionGeometryViewer(ReciprocalLatticeViewer):

  def create_viewer_panel (self) :
    self.viewer = DGVWindow(settings=self.settings, parent=self, size=(800,600),
      #orthographic=True
      )


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] experiments.json reflections.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) == 0 or len(reflections) == 0:
    parser.print_help()
    exit(0)

  reflections = reflections[0]

  imagesets = experiments.imagesets()

  import wxtbx.app
  a = wxtbx.app.CCTBXApp(0)
  a.settings = params
  f = DiffractionGeometryViewer(
    None, -1, "Diffraction Geometry viewer", size=(1024,768))
  f.load_models(imagesets, reflections)
  f.Show()
  a.SetTopWindow(f)
  #a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
