from __future__ import annotations

from dxtbx.model.tof_helpers import frame_to_tof_interpolator
from scitbx.array_family import flex


def centroid_px_to_mm(detector, scan, position, variance, sd_error):
    """Convenience function to calculate centroid in mm/rad from px"""

    # Get the pixel to millimeter function
    assert len(detector) == 1
    if scan is not None and scan.has_property("time_of_flight"):
        return tof_centroid_px_to_mm_panel(
            detector[0], scan, position, variance, sd_error
        )
    return centroid_px_to_mm_panel(detector[0], scan, position, variance, sd_error)


def centroid_px_to_mm_panel(panel, scan, position, variance, sd_error):
    """Convenience function to calculate centroid in mm/rad from px"""
    # Get the pixel to millimeter function

    if scan is not None and scan.has_property("time_of_flight"):
        return tof_centroid_px_to_mm_panel(panel, scan, position, variance, sd_error)

    pixel_size = panel.get_pixel_size()
    if scan is None or not scan.has_property("oscillation"):
        oscillation = (0, 0)
    else:
        oscillation = scan.get_oscillation(deg=False)
    scale = pixel_size + (oscillation[1],)
    scale2 = tuple(s * s for s in scale)

    if isinstance(position, tuple):
        # Convert Pixel coordinate into mm/rad
        x, y, z = position
        xy_mm = panel.pixel_to_millimeter((x, y))

        if scan is None or not scan.has_property("oscillation"):
            z_rad = 0
        else:
            z_rad = scan.get_angle_from_array_index(z, deg=False)

        # Set the position, variance and squared width in mm/rad
        # N.B assuming locally flat pixel to millimeter transform
        # for variance calculation.
        position_mm = xy_mm + (z_rad,)
        variance_mm = [var * s for var, s in zip(variance, scale2)]
        sd_error_mm = [sde * s for sde, s in zip(sd_error, scale2)]

    else:
        # Convert Pixel coordinate into mm/rad
        x, y, z = position.parts()
        xy_mm = panel.pixel_to_millimeter(flex.vec2_double(x, y))

        if scan is None or not scan.has_property("oscillation"):
            z_rad = flex.double(z.size(), 0)
        else:
            z_rad = scan.get_angle_from_array_index(z, deg=False)

        # Set the position, variance and squared width in mm/rad
        # N.B assuming locally flat pixel to millimeter transform
        # for variance calculation.
        x_mm, y_mm = xy_mm.parts()
        position_mm = flex.vec3_double(x_mm, y_mm, z_rad)
        v0, v1, v2 = variance.parts()
        variance_mm = flex.vec3_double(v0 * scale2[0], v1 * scale2[1], v2 * scale2[2])
        s0, s1, s2 = sd_error.parts()
        sd_error_mm = flex.vec3_double(s0 * scale2[0], s1 * scale2[1], s2 * scale2[2])

    # Return the stuff in mm/rad
    return position_mm, variance_mm, sd_error_mm


def tof_centroid_px_to_mm_panel(panel, scan, position, variance, sd_error):
    """Convenience function to calculate centroid in mm/tof from px"""
    assert scan is not None and scan.has_property("time_of_flight")

    pixel_size = panel.get_pixel_size()
    tof = scan.get_property("time_of_flight")  # (usec)
    frames = list(range(len(tof)))
    frame_to_tof = frame_to_tof_interpolator(frames, tof)

    if isinstance(position, tuple):
        # Convert Pixel coordinate into mm/rad
        x, y, z = position
        xy_mm = panel.pixel_to_millimeter((x, y))
        z_tof = flex.double(frame_to_tof(z))
        scale = pixel_size + (z_tof,)
        scale2 = tuple(s * s for s in scale)

        # Set the position, variance and squared width in mm/tof
        # N.B assuming locally flat pixel to millimeter transform
        # for variance calculation.
        position_mm = xy_mm + (z_tof,)
        variance_mm = [var * s for var, s in zip(variance, scale2)]
        sd_error_mm = [sde * s for sde, s in zip(sd_error, scale2)]

    else:
        # Convert Pixel coordinate into mm/tof
        x, y, z = position.parts()
        xy_mm = panel.pixel_to_millimeter(flex.vec2_double(x, y))
        z_tof = flex.double(frame_to_tof(z))
        scale = tuple(s * s for s in pixel_size)

        # Set the position, variance and squared width in mm/tof
        # N.B assuming locally flat pixel to millimeter transform
        # for variance calculation.
        x_mm, y_mm = xy_mm.parts()
        position_mm = flex.vec3_double(x_mm, y_mm, z_tof)
        v0, v1, v2 = variance.parts()
        variance_mm = flex.vec3_double(v0 * scale[0], v1 * scale[1], v2 * z_tof * z_tof)
        s0, s1, s2 = sd_error.parts()
        sd_error_mm = flex.vec3_double(s0 * scale[0], s1 * scale[1], s2 * z_tof * z_tof)

    return position_mm, variance_mm, sd_error_mm
