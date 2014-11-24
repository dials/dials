def analyse(reflections, detector, beam):
  from dials.array_family import flex

  p = reflections['panel']
  x, y = reflections['xyzobs.px.value'].parts()[:2]
  resolutions = flex.double(len(reflections), 0.0)

  # FIXME move this calculation to C++

  for j, r in enumerate(reflections):
    d = detector[p[j]].get_resolution_at_pixel(beam.get_s0(), (x[j], y[j]))
    resolutions[j] = d
  return len(resolutions) - (resolutions < 4).count(True) - \
    (resolutions > 40).count(True)

def work(filename, cl=[]):
  from dials.command_line.find_spots import phil_scope as params
  from dxtbx.datablock import DataBlockFactory
  from dials.array_family import flex
  interp = params.command_line_argument_interpreter()
  for cla in cl:
    params = params.fetch(interp.process(cla))
  datablock = DataBlockFactory.from_filenames([filename])[0]
  reflections = flex.reflection_table.from_observations(
    datablock, params.extract())
  detector = datablock.unique_detectors()[0]
  beam = datablock.unique_beams()[0]
  return analyse(reflections, detector, beam)

if __name__ == '__main__':
  import sys
  filename = sys.argv[1]
  params = sys.argv[2:]
  print work(filename, params)
