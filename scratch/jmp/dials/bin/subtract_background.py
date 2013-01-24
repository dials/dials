

def extract_reflection_profiles(paths):
    
    from scitbx import matrix

    import numpy

    from matplotlib import pylab, cm
    from mpl_toolkits.mplot3d import axes3d

    from dials.equipment import Beam
    from dials.equipment import Detector
    from dials.equipment import Goniometer
    from dials.integration import ReflectionMaskRoi
    from dials.integration import ReflectionMask
    from dials.integration import BackgroundIntensity
    from dials.integration import SubtractBackground
    from dials.io.pycbf_extra import search_for_image_volume
    from scitbx.array_family import flex
    
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
    n_sigma = 10

    # Set some values of a test reflection
    s1 = matrix.col((-0.0175199028171, -0.24775259748, 1.42844630118))
    phi_dash = 5.0
    phi = 5.83575672475
    sx = 117.588714455
    sy = 311.621428845
    sz = z0 + (phi - phi0) / dphi

    # Load the CBF image volume
    image_volume = search_for_image_volume(paths)
    
    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    # Create the beam, goniometer and detector
    beam = Beam(s0, l)
    gonio = Goniometer(m2, phi0, dphi, z0)
    detector = Detector(d1, d2, d3, (dx0, dy0), (px, py), (dx, dy), f)
    
    # Create the reflection mask regions of interest
    reflection_mask_roi = ReflectionMaskRoi(
                            beam, detector, gonio, 
                            n_sigma * sigma_d, 
                            n_sigma * sigma_m)
    region_of_interest = reflection_mask_roi.calculate(
                            flex.vec3_double(1, s1), 
                            flex.double(1, phi))
    
    reflection_mask = ReflectionMask(image_volume.shape)
    xyz = (sx, sy, sz-gonio.starting_frame)
    reflection_mask.create(flex.vec3_double(1, xyz), region_of_interest)
    
        
    #subtract_background = SubtractBackground(image_volume, roi_size, 0.1, 0.1)
    #subtract_background.subtract(flex.vec3_double(1, (sx, sy, sz-gonio.starting_frame)))
      

if __name__ == '__main__':

    import sys
    # Extract the reflection profiles
    extract_reflection_profiles("/home/upc86896/Projects/data/300k/ximg2700*.cbf")    

    
    
    
