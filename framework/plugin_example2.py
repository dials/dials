
def singleton(cls):
    instance = cls()
    instance.__call__ = lambda: instance
    return instance


class Factory(object):

    def __init__(self, plugins):
        self._plugins = plugins

    def create(self, name, *args, **kwargs):
        return self._plugins[name](*args, **kwargs)


@singleton
class Registry(object):

    def __init__(self):
        from collections import defaultdict
        self._plugins = defaultdict(dict)

    def append(self, mount, plugin):
        self._plugins[mount.name][plugin.name] = plugin

    def factory(self, mount_name):
        return Factory(self._plugins[mount_name])


class Extension(type):

    def __init__(self, name, bases, attrs):
        super(Extension, self).__init__(name, bases, attrs)



class Integrator(object):
    __metaclass__ = Extension

    name = 'integrator'


class XdsPlugin(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "xds"

    config = '''
        param1 = 10
          .type = int
          .help = "This is a parameter"
        param2 = 20.5
          .type = float
          .help = "This is another parameter"
    '''

    def __init__(self, sweep, crystal, params):
        print "Init XDS:"
        print "  ", sweep, crystal
        print "  ", params.integration.xds.param1
        print "  ", params.integration.xds.param2
        print "  ", params.lookup.gain

    def integrate(self):
        print "Call XDS"


class MosflmPlugin(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "mosflm"

    def __init__(self, sweep, crystal, params):
        print "Init Mosflm"

    def integrate(self):
        print "Call Mosflm"


class params:
    pass

params.integration = params()
params.integration.algorithm = "xds"
params.integration.xds = params()
params.integration.xds.param1 = 10
params.integration.xds.param2 = 20.5
params.lookup = params()
params.lookup.gain = "A Gain Map"

factory = Registry.factory('integrator')

integrator = factory.create(params.integration.algorithm, None, None, params)
integrator.integrate()

integrator = factory.create('mosflm', None, None, params)
integrator.integrate()
