
def find_cctbx_include_dirs():
    import os
    
    try:
        build = os.environ['LIBTBX_BUILD']
        source = os.environ['LIBTBX_SOURCE']
    except:
        print ('Can\'t find CCTBX source and build directories\n'
               'Try setting the following environment variables\n'
               ' - LIBTBX_BUILD\n'
               ' - LIBTBX_SOURCE\n')
        raise
        
        
    return [build + "/include", source + "/cctbx_project"]

def setup_package():
    from distutils.core import setup, Extension

    CCTBX_INCLUDE_DIRS = find_cctbx_include_dirs()

    setup(name='dials',
          version='0.1.0',
          description='Dials Integration',
          packages=['dials', 
                    'dials.util',
                    'dials.equipment',
                    'dials.geometry',
                    'dials.geometry.transform'],
                    
          ext_modules=[
            Extension('util_ext', [
                'dials/util/boost_python/util_ext.cc'],
                include_dirs=CCTBX_INCLUDE_DIRS,
                libraries=['boost_python']),
            Extension('equipment_ext', [
                'dials/equipment/boost_python/beam.cc',
                'dials/equipment/boost_python/detector.cc',
                'dials/equipment/boost_python/goniometer.cc',
                'dials/equipment/boost_python/equipment_ext.cc'],
                depends = [
                    'dials/equipment/beam.h',
                    'dials/equipment/detector.h',
                    'dials/equipment/goniometer.h'],
                include_dirs=CCTBX_INCLUDE_DIRS,
                libraries=['boost_python']),
            Extension('geometry_ext', [
                'dials/geometry/boost_python/detector_coordinate_system.cc',
                'dials/geometry/boost_python/reciprocal_lattice_coordinate_system.cc',
                'dials/geometry/boost_python/geometry_ext.cc'],
                depends = [
                    'dials/geometry/detector_coordinate_system.h',
                    'dials/geometry/reciprocal_lattice_coordinate_system.h'],
                include_dirs=CCTBX_INCLUDE_DIRS,
                libraries=['boost_python']),
            Extension('transform_ext', [
                'dials/geometry/transform/boost_python/from_detector_to_beam_vector.cc',
                'dials/geometry/transform/boost_python/from_beam_vector_to_detector.cc',
                'dials/geometry/transform/boost_python/from_hkl_to_beam_vector.cc',
                'dials/geometry/transform/boost_python/from_hkl_to_detector.cc',
                'dials/geometry/transform/boost_python/transform_ext.cc'],
                depends = [
                    'dials/geometry/transform/from_detector_to_beam_vector.h',
                    'dials/geometry/transform/from_beam_vector_to_detector.h',
                    'dials/geometry/transform/from_hkl_to_beam_vector.h',
                    'dials/geometry/transform/from_hkl_to_detector.h'],
                include_dirs=CCTBX_INCLUDE_DIRS,
                libraries=['boost_python'])]
     )

if __name__ == '__main__':
    setup_package()
