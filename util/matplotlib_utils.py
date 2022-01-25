from __future__ import annotations

from matplotlib import ticker


@ticker.FuncFormatter
def resolution_formatter(d_star_sq, pos):
    return f"{1 / d_star_sq ** 0.5:.2f}" if d_star_sq > 0 else ""
