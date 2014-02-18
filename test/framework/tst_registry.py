from __future__ import division

class Test:

  def __init__(self):
    pass

  def run(self):
    self.tst_singleton()
    self.tst_get_interfaces()
    self.tst_get_extensions()
    self.tst_dict_access()

  def tst_singleton(self):
    from dials.framework.registry import Registry
    r1 = Registry()
    r2 = Registry()
    assert(r1 is r2)
    print 'OK'

  def tst_get_interfaces(self):

    from dials.framework.registry import Registry
    registry = Registry()

    assert(len(registry) == 3)
    interfaces = list(registry.interfaces())
    assert(len(list(interfaces)) == 3)

    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntegrationIface

    assert(CentroidIface in interfaces)
    assert(BackgroundIface in interfaces)
    assert(IntegrationIface in interfaces)

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

    centroid_iface = registry['centroid']
    background_iface = registry['background']
    integration_iface = registry['integration']

    from dials.interfaces import IntegrationIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface

    assert(integration_iface == IntegrationIface)
    assert(background_iface == BackgroundIface)
    assert(centroid_iface == CentroidIface)

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
