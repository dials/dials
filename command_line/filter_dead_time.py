# LIBTBX_SET_DISPATCHER_NAME dev.dials.filter_dead_time

from __future__ import absolute_import, division, print_function

import iotbx.phil

help_message = '''\
'''

phil_scope = iotbx.phil.parse('''\
dead_time = 0
  .help = "Detector dead time in ms, assumed to be at the end of the exposure time."
  .type = float(value_min=0)
reject_fraction = 0
  .help = "Reject reflections which overlap by more than the given fraction"
          "with the dead region of the image."
  .type = float(value_min=0, value_max=1)
output {
  reflections = filtered.pickle
    .type = path
}
''', process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = '%s [options] experiments.json integrated.pickle' %(
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
    exit()

  experiment = experiments[0]
  reflections = reflections[0]

  sel = reflections.get_flags(reflections.flags.integrated)
  reflections = reflections.select(sel)

  goniometer = experiment.goniometer
  beam = experiment.beam

  m2 = goniometer.get_rotation_axis()
  s0 = beam.get_s0()

  from dials.array_family import flex
  phi1 = flex.double()
  phi2 = flex.double()

  phi_range = reflections.compute_phi_range(
    goniometer.get_rotation_axis(),
    beam.get_s0(),
    experiment.profile.sigma_m(deg=False),
    experiment.profile.n_sigma())
  phi1, phi2 = phi_range.parts()

  scan = experiment.scan
  exposure_time = scan.get_exposure_times()[0]
  assert scan.get_exposure_times().all_eq(exposure_time)
  phi_start, phi_width = scan.get_oscillation(deg=False)
  phi_range_dead = phi_width * (params.dead_time/1000) / exposure_time

  sel_good = flex.bool(len(reflections), True)

  start, end = scan.get_array_range()
  for i in range(start, end):
    phi_dead_start = phi_start + (i+1) * phi_width - phi_range_dead
    phi_dead_end = phi_dead_start + phi_range_dead

    left = phi1.deep_copy()
    left.set_selected(left < phi_dead_start, phi_dead_start)

    right = phi2.deep_copy()
    right.set_selected(right > phi_dead_end, phi_dead_end)

    overlap = (right - left)/(phi2-phi1)

    sel = overlap > params.reject_fraction

    sel_good.set_selected(sel, False)
    print('Rejecting %i reflections from image %i' %(sel.count(True), i))

  print('Keeping %i reflections (rejected %i)' %(
    sel_good.count(True), sel_good.count(False)))

  from libtbx import easy_pickle
  easy_pickle.dump(params.output.reflections, reflections.select(sel_good))



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
