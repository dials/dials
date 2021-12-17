from matplotlib import ticker


@ticker.FuncFormatter
def resolution_formatter(d_star_sq, pos):
    # d = 1 / d_star_sq ** 0.5 if d_star_sq > 0 else None
    return f"{1 / d_star_sq ** 0.5:.2f}" if d_star_sq > 0 else ""
