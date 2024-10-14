from dials.algorithms.beam_position import MidpointMethodSolver
from dials.algorithms.beam_position.plot import Figure


def compute_beam_position(image, params, index=None):

    method_x, method_y = resolve_projection_methods(params)

    if method_x == 'midpoint':
        solver_x = MidpointMethodSolver(image, params, axis='x')

    if method_y == 'midpoint':
        solver_y = MidpointMethodSolver(image, params, axis='y')

    x = solver_x.find_beam_position()
    y = solver_y.find_beam_position()

    fig = Figure('fig.png')

    fig.plot_main(image, params, beam_position=(x, y))
    solver_x.plot(fig)
    solver_y.plot(fig)

    print('Solved xy', x, y)
    fig.save_and_close()

    # if params.bad_pixel_threshold:
    #    image[np.where(image > params.bad_pixel_threshold)] = 0

    # If plot make plotting skeleton

    # Compute and plot along x
    # Compute and plot along y

    # Plot along the central region with excluded points

    # Save figure and return


def resolve_projection_methods(params):
    """
    Resolves which method to use along x and y direction.
    """

    p = params.projection

    # Only method_x supplied
    if len(p.method_x) == 1 and len(p.method_y) > 1:
        p.method_y = p.method_x
    # Only method_y supplied
    elif len(p.method_y) == 1 and len(p.method_x) > 1:
        p.method_x = p.method_y

    # Check if both options are set
    if not (len(p.method_x) == 1 and len(p.method_y) == 1):
        msg = "Projection method_x and method_y require a single option"
        raise ValueError(msg)

    options = ['midpoint', 'maximum', 'inversion']

    if not (p.method_x[0] in options):
        msg = "Wrong projection method along the x: '%s' \n" % p.method_x[0]
        msg += "Correct options: 'midpoint', 'maximum', or 'inversion'"
        raise ValueError(msg)

    if not (p.method_x[0] in options):
        msg = "Wrong projection method along the x: '%s' \n" % p.method_x[0]
        msg += "Correct options: 'midpoint', 'maximum', or 'inversion'"
        raise ValueError(msg)

    return p.method_x[0], p.method_y[0]
