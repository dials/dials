from dials.util.export_mtz import export_mtz

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 4:
    raise RuntimeError, '%s integrated.pickle experiments.json hklout.mtz'

  import cPickle as pickle
  from dials.array_family import flex
  from dials.model.experiment.experiment_list import ExperimentListFactory

  integrated_data = pickle.load(open(sys.argv[1], 'rb'))
  experiment_list = ExperimentListFactory.from_json_file(sys.argv[2])

  m = export_mtz(integrated_data, experiment_list, sys.argv[3])
  m.show_summary()
