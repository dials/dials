from __future__ import division
def run(args):
  from dials.util.command_line import Importer
  assert len(args) % 2 == 0

  experiment_filenames = args[:len(args)//2]
  reflection_filenames = args[len(args)//2:]

  from dxtbx.model.experiment.experiment_list import ExperimentList
  from dials.array_family import flex

  experiment_list = ExperimentList()
  reflection_table = None
  i_cryst = -1


  detector = None
  beam = None
  scan = None
  goniometer = None

  for i in range(len(experiment_filenames)):
    print experiment_filenames[i], reflection_filenames[i]
    importer = Importer([experiment_filenames[i], reflection_filenames[i]])
    assert importer.experiments is not None
    assert importer.reflections is not None
    assert len(importer.reflections) == 1
    experiments = importer.experiments
    n_expts = len(experiments)
    reflections = importer.reflections[0]
    assert len(experiments.detectors()) == 1
    assert len(experiments.beams()) == 1
    assert len(experiments.goniometers()) == 1
    assert len(experiments.scans()) == 1
    if detector is None:
      detector = experiments.detectors()[0]
    else:
      try:
        assert detector == experiments.detectors()[0]
      except AssertionError:
        print detector, experiments.detectors()[0]

    assert (reflections['id'] < 0).count(True) == 0
    updated_ref_ids = flex.int(len(reflections), -1)
    if scan is None:
      scan = experiments.scans()[0]
    else:
      try:
        assert scan == experiments.scans()[0]
      except AssertionError:
        print scan, experiments.scans()[0]

    if goniometer is None:
      goniometer = experiments.goniometers()[0]
    else:
      try:
        assert goniometer == experiments.goniometers()[0]
      except AssertionError:
        print goniometer, experiments.goniometers()[0]

    if beam is None:
      beam = experiments.beams()[0]
    else:
      try:
        assert beam == experiments.beams()[0]
      except AssertionError:
        print beam, experiments.beams()[0]

    for i in range(n_expts):
      experiments[i].detector = detector
      experiments[i].beam = beam
      experiments[i].goniometer = goniometer
      #experiments[i].scan = scan
      i_cryst += 1
      print i_cryst
      sel = (reflections['id'] == i)
      assert sel.count(True) > 0
      print sel.count(True)
      updated_ref_ids.set_selected(sel, i_cryst)
      print (updated_ref_ids == i_cryst).count(True)
    reflections['id'] = updated_ref_ids
    reflections = reflections.select(reflections['id'] > -1)
    #assert (reflections['id'] < 0).count(True) == 0
    experiment_list.extend(experiments)
    if reflection_table is None:
      reflection_table = reflections
    else:
      reflection_table.extend(reflections)

  print
  print len(experiments)
  for i in range(len(experiments)):
    print (reflection_table['id'] == i).count(True)

  assert (reflection_table['id'] > -1).count(True) == len(reflection_table)

  from dials.model.serialize import dump as dials_dump
  from dxtbx.serialize import dump as dxtbx_dump
  dxtbx_dump.experiment_list(experiment_list, 'combined_experiments.json')
  dials_dump.reflections(reflection_table, 'combined_reflections.pickle')

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
