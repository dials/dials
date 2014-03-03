from __future__ import division

class Test:

  def __init__(self):
    pass

  def run(self):
    self.tst_singleton()
    self.tst_get_interfaces()
    self.tst_get_extensions()
    self.tst_dict_access()
    self.tst_init_ext()

  def tst_singleton(self):
    from dials.framework.registry import Registry
    r1 = Registry()
    r2 = Registry()
    assert(r1 is r2)
    print 'OK'

  def tst_get_interfaces(self):

    from dials.framework.registry import Registry
    registry = Registry()

    assert(len(registry) == 4)
    interfaces = list(registry.interfaces())
    assert(len(list(interfaces)) == 4)

    from dials.interfaces import SpotFinderThresholdIface
    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntensityIface

    assert(CentroidIface in interfaces)
    assert(SpotFinderThresholdIface in interfaces)
    assert(BackgroundIface in interfaces)
    assert(IntensityIface in interfaces)

    print 'OK'

  def tst_get_extensions(self):

    from dials.framework.registry import Registry
    registry = Registry()
    extensions = list(registry.all_extensions())

    from dials.extensions.simple_centroid_ext import SimpleCentroidExt
    from dials.extensions.null_background_ext import NullBackgroundExt
    from dials.extensions.fable_background_ext import FableBackgroundExt
    from dials.extensions.xds_background_ext import XdsBackgroundExt
    from dials.extensions.summation_3d_integration_ext import Summation3dIntegrationExt
    from dials.extensions.profile_fitting_rs_integration_ext import ProfileFittingRSIntegrationExt

    assert(SimpleCentroidExt in extensions)
    assert(NullBackgroundExt in extensions)
    assert(FableBackgroundExt in extensions)
    assert(XdsBackgroundExt in extensions)
    assert(Summation3dIntegrationExt in extensions)
    assert(ProfileFittingRSIntegrationExt in extensions)

    print 'OK'

  def tst_dict_access(self):

    from dials.framework.registry import Registry
    registry = Registry()

    centroid_iface = registry.interface('integration.centroid')
    background_iface = registry.interface('integration.background')
    integration_iface = registry.interface('integration.intensity')

    from dials.interfaces import IntensityIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface

    assert(integration_iface == IntensityIface)
    assert(background_iface == BackgroundIface)
    assert(centroid_iface == CentroidIface)

    centroid_ext = registry['integration.centroid']
    background_ext = registry['integration.background']
    integration_ext = registry['integration.intensity']

    from dials.extensions.simple_centroid_ext import SimpleCentroidExt
    from dials.extensions.null_background_ext import NullBackgroundExt
    from dials.extensions.fable_background_ext import FableBackgroundExt
    from dials.extensions.xds_background_ext import XdsBackgroundExt
    from dials.extensions.summation_3d_integration_ext import Summation3dIntegrationExt
    from dials.extensions.profile_fitting_rs_integration_ext import ProfileFittingRSIntegrationExt

    assert(integration_ext == Summation3dIntegrationExt)
    assert(background_ext == XdsBackgroundExt)
    assert(centroid_ext == SimpleCentroidExt)
    print 'OK'

  def tst_init_ext(self):

    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntensityIface
    from dials.framework.registry import init_ext
    algorithm1 = init_ext("integration.centroid", None)
    algorithm2 = init_ext("integration.background", None)
    algorithm3 = init_ext("integration.intensity", None)

    assert(isinstance(algorithm1, CentroidIface))
    assert(isinstance(algorithm2, BackgroundIface))
    assert(isinstance(algorithm3, IntensityIface))

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
