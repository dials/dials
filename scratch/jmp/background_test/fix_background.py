def generate_reflections():
  from dials.algorithms.simulation.generate_test_reflections import main
  from dials.algorithms.simulation.generate_test_reflections import \
    master_phil
  from libtbx.phil import command_line
  cmd = command_line.argument_interpreter(master_params = master_phil)
  working_phil = cmd.process_and_fetch(args = ["""
    nrefl = 1000
    shoebox_size {
      x = 10
      y = 10
      z = 10
    }
    spot_size {
      x = 1
      y = 1
      z = 1
    }
    spot_offset {
      x = -0.5
      y = -0.5
      z = -0.5
    }
    mask_nsigma = 3.0
    counts = 0
    background = 10
    pixel_mask = all *static precise
    background_method = *xds mosflm
    integration_methpd = *xds mosflm
    output {
      over = None
      under = None
      all = all_refl.pickle
    }
    rotation {
      axis {
        x = 0
        y = 0
        z = 1
      }
      angle = 0
    }

    """])
  main(working_phil.extract())
  print 'OK'


if __name__ == '__main__':
  generate_reflections()
