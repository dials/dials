from scitbx.array_family import flex


def calculate_reflection_mask_roi(beam, detector, goniometer,
                                  diffracted_beam_vector,
                                  rotation_angle,
                                  delta_divergence,
                                  delta_mosaicity):

    from dials.geometry import XdsCoordinateSystem
    from dials.geometry.transform import FromXdsToDetector
    from dials.geometry.transform import FromXdsE3ToPhi
    from math import floor, ceil
 
    # Create the coordinate system for the reflection
    xds_coordinate_system = XdsCoordinateSystem(beam.direction,
                                                diffracted_beam_vector,
                                                goniometer.rotation_axis,
                                                rotation_angle)
    
    # Create the transformer from the xds coordinate system to the detector
    from_xds_to_detector = FromXdsToDetector(xds_coordinate_system,
                                             diffracted_beam_vector,
                                             detector)

    # Create the transformer from Xds E3 to rotation angle                                             
    from_xds_e3_to_phi = FromXdsE3ToPhi(xds_coordinate_system.zeta, 
                                        rotation_angle)
                                                                                     
    # Find the detector coordinates at the following xds coordinates:
    #   (-delta_d/2, -delta_d/2, 0) 
    #   (+delta_d/2, -delta_d/2, 0)
    #   (-delta_d/2, +delta_d/2, 0) 
    #   (+delta_d/2, +delta_d/2, 0)
    point = delta_divergence / 2.0
    x1, y1 = from_xds_to_detector.apply((-point, -point, 0.0))
    x2, y2 = from_xds_to_detector.apply((+point, -point, 0.0))
    x3, y3 = from_xds_to_detector.apply((-point, +point, 0.0))
    x4, y4 = from_xds_to_detector.apply((+point, +point, 0.0))

    # Get the image volume z coordinates (zero based) at the following XDS e3 
    # coordinates: -delta_m/2, +delta_m/2
    z1 = goniometer.get_frame_from_angle(
                from_xds_e3_to_phi.apply(
                    -delta_mosaicity / 2.0)) - goniometer.starting_frame
    z2 = goniometer.get_frame_from_angle(
                from_xds_e3_to_phi.apply(
                    +delta_mosaicity / 2.0)) - goniometer.starting_frame

    # Return the roi in the following form:
    # (minx, maxx, miny, maxy, minz, maxz)
    # Min's are rounded down to the nearest integer, Max's are rounded up
    minx = int(floor(min([x1, x2, x3, x4])))
    maxx = int(ceil (max([x1, x2, x3, x4])))
    miny = int(floor(min([y1, y2, y3, y4])))
    maxy = int(ceil (max([y1, y2, y3, y4])))
    minz = int(floor(min([z1, z2])))
    maxz = int(ceil (max([z1, z2])))
    return (minx, maxx, miny, maxy, minz, maxz)    
    
def extract_reflection_profiles(paths):
    
    from scitbx import matrix

    import numpy

    from matplotlib import pylab, cm
    from mpl_toolkits.mplot3d import axes3d

    from dials.equipment import Beam
    from dials.equipment import Detector
    from dials.equipment import Goniometer
    from dials.geometry import XdsCoordinateSystem
    from dials.geometry.transform import FromDetectorToBeamVector
    from dials.geometry.transform import FromBeamVectorToXds
    from dials.integration import ReflectionMask
    from dials.io.pycbf_extra import search_for_image_volume

    
    # Set a load of parameters from the GXPARM file
    z0 = 1
    phi0 = 1.0
    dphi = 1.0
    s0 = matrix.col((0.013142, 0.002200, 1.450476))
    dx0 = 244.836136
    dy0 = 320.338531
    dx = 487
    dy = 619
    px = 0.172000
    py = 0.172000
    d1 = matrix.col((0.000000, -0.939693, -0.342020))
    d2 = matrix.col((1.000000,  0.000000,  0.000000))
    d3 = matrix.col((0.000000, -0.342020,  0.939693))
    m2 = matrix.col((0.999975, -0.001289, -0.006968))
    f  = 122.124901
    l  = 0.689400
    sigma_d = 0.034
    sigma_m = 0.082

    # Set some values of a test reflection
    s1 = matrix.col((-0.0175199028171, -0.24775259748, 1.42844630118))
    phi_dash = 5.0
    phi = 5.83575672475
    sx = 117.588714455
    sy = 311.621428845
    sz = z0 + (phi - phi0) / dphi

    # Load the CBF image volume
    image_volume = search_for_image_volume(paths)
    image_volume = flex.int(image_volume)
    image_volume_copy = image_volume.deep_copy()

    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    # Create the beam, goniometer and detector
    beam = Beam(s0, l)
    gonio = Goniometer(m2, phi0, dphi, z0)
    detector = Detector(d1, d2, d3, (dx0, dy0), (px, py), (dx, dy), f)

    #roi = calculate_reflection_mask_roi(beam, detector, gonio, s1, phi,
    #                              10 * sigma_d,
    #                              10 * sigma_m)
    
    
    #print roi
    
    from dials.integration import ReflectionMaskRoi
    
    reflection_mask_roi = ReflectionMaskRoi(beam, detector, gonio, 10 * sigma_d, 10 * sigma_m)
    
    roi = reflection_mask_roi.calculate(flex.vec3_double(1, s1), flex.double(1, phi))
    
    print roi[0]
    
    #print roi
    
    
    #rcs = XdsCoordinateSystem(s0, s1, m2, phi)

    #roi_size = (1, 4, 4)
    #x0 = int(sx) - roi_size[2]
    #x1 = int(sx) + roi_size[2]
    #y0 = int(sy) - roi_size[1]
    #y1 = int(sy) + roi_size[1]
    #z0 = int(sz) - int(gonio.starting_frame) - roi_size[1]
    #z1 = int(sz) - int(gonio.starting_frame) + roi_size[1]    

    
    reflection_mask = ReflectionMask(image_volume.all())
    #reflection_mask = ReflectionMask(image_volume.all(), roi_size=roi_size)
    #reflection_mask.create(flex.vec3_double(1, (sx, sy, sz-z0)))
    xyz = (sx, sy, sz - z0)
    
    reflection_mask.create(flex.vec3_double(1, xyz), roi)

    mask = reflection_mask.mask
    
    from matplotlib import pylab, cm
    mask = mask.as_numpy_array()

    for i in range(0, 10):
        pylab.imshow(mask[i,:,:], cmap=cm.Greys_r)
        pylab.show()
    
if __name__ == '__main__':

    import sys
    # Extract the reflection profiles
    extract_reflection_profiles(sys.argv[1])    

    
    
    
