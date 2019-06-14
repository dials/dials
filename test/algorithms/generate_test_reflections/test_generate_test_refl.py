from __future__ import absolute_import, division, print_function


def test_generate_test_reflections(tmpdir):
    tmpdir.chdir()
    from dials.algorithms.simulation.generate_test_reflections import main
    from dials.algorithms.simulation.generate_test_reflections import master_phil
    from libtbx.phil import command_line

    cmd = command_line.argument_interpreter(master_phil=master_phil)
    working_phil = cmd.process_and_fetch(
        args=[
            """nrefl = 10
shoebox_size {
  x = 20
  y = 20
  z = 20
}
spot_size {
  x = 1
  y = 3
  z = 1
}
spot_offset {
  x = -0.5
  y = -0.5
  z = -0.5
}
mask_nsigma = 3.0
counts = 10000
background = 0
background_a = 10
background_b = 0.1
background_c = -0.1
background_d = 0
pixel_mask = *all static precise
background_method = *xds mosflm
integration_methpd = *xds mosflm
output {
  over = None
  under = None
  all = all_refl.refl
}
rotation {
  axis {
    x = 0
    y = 0
    z = 1
  }
  angle = 45
}

"""
        ]
    )
    main(working_phil.extract())
    assert (tmpdir / "all_refl.refl").check()
