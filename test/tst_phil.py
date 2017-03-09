
from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")
    self.datablock_path = join(self.path, "datablock.json")
    self.experiments_path = join(self.path, "experiments.json")
    self.reflections_path = join(self.path, "integrated.pickle")

  def run(self):
    from dials.phil import parse

    phil_scope = parse('''
      input {
        reflections = %s
          .type = reflection_table
        datablock = %s
          .type = datablock
        experiments = %s
          .type = experiment_list
      }
    ''' % (self.reflections_path,
           self.datablock_path,
           self.experiments_path))

    params = phil_scope.extract()
    from dials.array_family import flex
    from dxtbx.datablock import DataBlock
    from dxtbx.model.experiment_list import ExperimentList
    assert(params.input.reflections.filename == self.reflections_path)
    assert(params.input.datablock.filename == self.datablock_path)
    assert(params.input.experiments.filename == self.experiments_path)
    assert(isinstance(params.input.reflections.data, flex.reflection_table))
    assert(len(params.input.datablock.data) == 1)
    assert(isinstance(params.input.datablock.data[0], DataBlock))
    assert(isinstance(params.input.experiments.data, ExperimentList))
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
