from __future__ import division

def run(args):
  from matplotlib import pyplot
  from dials.util.command_line import Importer
  from scitbx.array_family import flex
  importer = Importer(args)
  #datablocks = importer.datablocks
  #assert len(datablocks) == 1
  ##imagesets = datablocks[0].extract_imagesets()
  observed_xy = flex.vec2_double()
  for reflection_list in importer.reflections:
    for refl in reflection_list:
      centroid = refl['xyzobs.px.value']
      if centroid != (0,0,0):
        x, y = centroid[:2]
        observed_xy.append((x,y))
  obs_x, obs_y = observed_xy.parts()
  fig = pyplot.figure()
  pyplot.scatter(obs_x, obs_y, marker='.', c='black', lw=0, s=0.1)
  pyplot.axes().set_aspect('equal')
  pyplot.xlim((0, pyplot.xlim()[1]))
  pyplot.ylim((0, pyplot.ylim()[1]))
  pyplot.savefig("virtual_powder.png",
                 size_inches=(10,10), dpi=600)
  #pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
