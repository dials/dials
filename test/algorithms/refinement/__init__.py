from __future__ import absolute_import, division, print_function

geometry_phil = """geometry.parameters
    .help = "Set up an experimental model for refinement test case"
{
    random_seed = 2
        .type = int
    beam
    {
        wavelength
        {
            random = True
                .type = bool
            range = 0.8 1.2
                .type = floats(size=2)
            value = None
                .type = float
        }
        direction
        {
            method = *inclination close_to exactly
                .type = choice
            inclination
                .help = "inclination in the X-Z plane"
            {
                random = True
                    .type = bool
                angle = 0.5
                    .type = float
            }
            close_to
            {
                direction = None
                    .type = floats(size=3)
                sd = 0.5
                    .type = float
            }
            exactly = None
                .type = floats(size=3)
        }
    }
    crystal
    {
        space_group_symbol = "P 1"
            .type = str
        a
        {
            length
            {
                random = True
                    .type = bool
                range = 10 20
                    .type = floats(size=2)
                value = None
                    .type = float
            }
            direction
            {
                method = *close_to exactly
                    .type = choice
                close_to
                {
                        direction = 1., 0., 0.
                            .type = floats(size=3)
                        sd = 0.5
                            .type = float
                }
                exactly
                {
                        direction = 1., 0., 0.
                            .type = floats(size=3)
                }
            }


        }
        b
        {
            length
            {
                random = True
                    .type = bool
                range = 10 20
                    .type = floats(size=2)
                value = None
                    .type = float
            }
            direction
            {
                method = *close_to exactly
                    .type = choice
                close_to
                {
                        direction = 0., 1., 0.
                            .type = floats(size=3)
                        sd = 0.5
                            .type = float
                }
                exactly
                {
                        direction = 0., 1., 0.
                            .type = floats(size=3)
                }
            }


        }
        c
        {
            length
            {
                random = True
                    .type = bool
                range = 10 20
                    .type = floats(size=2)
                value = None
                    .type = float
            }
            direction
            {
                method = *close_to exactly
                    .type = choice
                close_to
                {
                        direction = 0., 0., 1.
                            .type = floats(size=3)
                        sd = 0.5
                            .type = float
                }
                exactly
                {
                        direction = 0., 0., 1.
                            .type = floats(size=3)
                }
            }


        }
    }
    goniometer
    {
        axis = 1. 0. 0.
            .type = floats(size=3)
    }
    detector
    {
        directions
        {
            method = *close_to exactly
                .type=choice
            close_to
            {
                dir1 = 1., 0., 0.
                    .type = floats(size=3)
                norm = 0., 0., 1.
                    .type = floats(size=3)
                sd = 0.5
                    .type = float
            }
            exactly
            {
                dir1 = 1., 0., 0.
                    .type = floats(size=3)
                norm = 0., 0., 1.
                    .type = floats(size=3)
            }
        }
        centre
        {
            method = *close_to exactly
                .type=choice
            close_to
            {
                value = 0., 0., 200.
                    .type = floats(size=3)
                sd = 0.5
                    .type = float
            }
            exactly
            {
                value = 0., 0., 200.
                    .type = floats(size=3)
            }
        }
        npx_fast = 1475
            .type = int
        npx_slow = 1679
            .type = int
        pix_size = 0.172
            .type = float
    }
}"""


minimiser_phil = """minimiser.parameters
    .help = "Set up an minimiser for refinement test case"
{
    engine = SimpleLBFGS LBFGScurvs *GaussNewton
        .type = choice
    logfile = tst_orientation_refinement.log
        .type = path
}"""
