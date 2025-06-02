from __future__ import annotations

from dials.algorithms.beam_position import (
    InversionMethodSolver,
    MaximumMethodSolver,
    MidpointMethodSolver,
)
from dials.algorithms.beam_position.plot import Figure


def compute_beam_position(image, params, image_index=None, imageset_index=None):
    method_x, method_y = resolve_projection_methods(params)

    if method_x == "midpoint":
        solver_x = MidpointMethodSolver(image, params, axis="x")
    elif method_x == "maximum":
        solver_x = MaximumMethodSolver(image, params, axis="x")
    elif method_x == "inversion":
        solver_x = InversionMethodSolver(image, params, axis="x")

    if method_y == "midpoint":
        solver_y = MidpointMethodSolver(image, params, axis="y")
    elif method_y == "maximum":
        solver_y = MaximumMethodSolver(image, params, axis="y")
    elif method_y == "inversion":
        solver_y = InversionMethodSolver(image, params, axis="y")

    x = solver_x.find_beam_position()
    y = solver_y.find_beam_position()

    if image_index is not None:
        image_str = "_img_%05d" % image_index
    else:
        image_str = ""

    if imageset_index is not None:
        imageset_str = "_imageset_%05d" % imageset_index
    else:
        imageset_str = ""

    if params.projection.plot:
        fig = Figure(f"beam_position{imageset_str}{image_str}.png")

        fig.plot_main(image, params, beam_position=(x, y))
        solver_x.plot(fig)
        solver_y.plot(fig)

        fig.save_and_close()

    if not params.projection.per_image:
        print(f"{x:.2f}, {y:.2f}")
    return x, y


def resolve_projection_methods(params):
    """
    Resolves which method to use along the x and the y direction.
    """

    options = ["midpoint", "maximum", "inversion"]

    p = params.projection
    resolved_method_x = None
    resolved_method_y = None

    if params.method[0] in options and (len(params.method) == 1):
        resolved_method_x = params.method[0]
        resolved_method_y = params.method[0]

    # Parameters `method_x` and `method_y` overwrite the `method` parameter
    # Only `method_x` supplied
    # Note that by default `method_x` and `method_y` hold a list of all options
    if len(p.method_x) == 1 and len(p.method_y) > 1:
        resolved_method_x = p.method_x[0]

        # If no `method` is set, assume same method along x and y
        if params.method == "default":
            resolved_method_y = p.method_x[0]

    # Only `method_y` supplied
    elif len(p.method_y) == 1 and len(p.method_x) > 1:
        resolved_method_y = p.method_y[0]

        # If no `method` is set, assume same method along x and y
        if params.method == "default":
            resolved_method_x = p.method_y[0]

    if resolved_method_x not in options:
        msg = f"Wrong projection method along the x: {resolved_method_x} \n"
        msg += "Correct options: 'midpoint', 'maximum', or 'inversion'"
        raise ValueError(msg)

    if resolved_method_y not in options:
        msg = f"Wrong projection method along the y: {resolved_method_y} \n"
        msg += "Correct options: 'midpoint', 'maximum', or 'inversion'"
        raise ValueError(msg)

    return resolved_method_x, resolved_method_y
