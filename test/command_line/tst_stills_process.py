from __future__ import division
from dials.array_family import flex # import dependency


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "image_examples/SACLA_MPCCD_Cheetah")

  def run(self):
    self.test_sacla_h5()

  def test_sacla_h5(self):
    from os.path import join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    assert exists(join(self.path, 'run266702-0-subset.h5'))

    f = open("process.phil", 'w')
    f.write("""
      indexing {
        known_symmetry {
          space_group = P43212
          unit_cell = 78.9 78.9 38.1 90 90 90
        }
        method = real_space_grid_search
        refinement_protocol.d_min_start = 2.2
      }

      spotfinder {
        filter.min_spot_size = 2
      }

      refinement {
        parameterisation {
          beam.fix = all
          detector.fix_list = Dist,Tau1
          auto_reduction {
            action = fix
            min_nref_per_parameter = 1
          }
          crystal {
            unit_cell {
              restraints {
                tie_to_target {
                  values = 78.9,78.9,38.1,90,90,90
                  sigmas = 1,1,1,0,0,0
                  apply_to_all = True
                }
              }
            }
          }
        }
      }
      integration {
        integrator = stills
        profile.fitting = False
        background {
          algorithm = simple
          simple {
            model.algorithm = linear2d
            outlier.algorithm = tukey
          }
        }
      }
      profile {
        gaussian_rs {
          min_spots.overall = 0
        }
      }
      """)
    f.close()

    # Call dials.stills_process
    result = easy_run.fully_buffered([
      'dials.stills_process',
      join(self.path, 'run266702-0-subset.h5'),
      'process.phil',
    ]).raise_if_errors()
    result.show_stdout()

    import cPickle as pickle
    for result, n_refls in zip(["idx-run266702-0-subset_00000_integrated.pickle",
                                "idx-run266702-0-subset_00001_integrated.pickle"],
                                [range(109,114), range(80,85)]): # large ranges to handle platform-specific differences
      table = pickle.load(open(result, 'rb'))
      assert len(table) in n_refls, len(table)
      assert 'id' in table
      assert (table['id'] == 0).count(False) == 0
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
