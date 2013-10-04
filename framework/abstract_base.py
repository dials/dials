


from abc import ABCMeta, abstractmethod

class Factory(object):

    def __init__(self, plugins):
        self._plugins = plugins

    def create(self, name, *args, **kwargs):
        return self._plugins[name](*args, **kwargs)

def singleton(cls):
    instance = cls()
    instance.__call__ = lambda: instance
    return instance

@singleton
class Registry:

    def __init__(self):
        from collections import defaultdict
        self._mount_points = dict()

    def append(self, mount):
        self._mount_points[mount] = mount

    def factory(self, mount_name):
        mount = self._mount_points[mount_name]
        plugins = dict()
        for sc in  mount.__subclasses__():
            plugins[sc.name] = sc
        return Factory(plugins)



class Extension(ABCMeta):

    def __init__(self, name, bases, attrs):
        super(Extension, self).__init__(name, bases, attrs)
        if not hasattr(self, 'registered'):
            self.registered = True
            Registry.append(self)


class Integrator(object):

    __metaclass__ = Extension

    @abstractmethod
    def integrate(self):
        pass


class XdsIntegrator(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "xds"

    def __init__(self, sweep, crystal, params):
        super(XdsIntegrator, self).__init__()
        print "Init XDS:"

    def integrate(self):
        print "Integrate"



class MosflmIntegrator(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "mosflm"

    def __init__(self, sweep, crystal, params):
        super(MosflmIntegrator, self).__init__()
        print "Init Mosflm:"

    def integrate(self):
        print "Integrate"


class params:
    pass

params.integration = params()
params.integration.algorithm = "xds"
params.integration.xds = params()
params.integration.xds.param1 = 10
params.integration.xds.param2 = 20.5
params.lookup = params()
params.lookup.gain = "A Gain Map"

factory = Registry.factory(Integrator)

integrator = factory.create(params.integration.algorithm, None, None, params)
print integrator.integrate()


integrator = factory.create('mosflm', None, None, params)
print integrator.integrate()


class HighLevelInterface(object):

    __metaclass__ = ABCMeta

    def __init__(self, maxiter=20):
        self._prepared = False
        self._done = False
        self._finished = False
        self._maxiter = maxiter

    def __call__(self):
        niter = [0, 0, 0]
        while not self.finished:
            niter[1:] = [0, 0]
            while not self.done:
                niter[2] = 0
                while not self.prepared:
                    self.prepare()
                    niter[2] += 1
                    if niter[2] >= self._maxiter:
                        raise RuntimeError('niter > maxiter')
                self.do()
                niter[1] += 1
                if niter[1] >= self._maxiter:
                    raise RuntimeError('niter > maxiter')
            self.finish()
            niter[0] += 1
            if niter[0] >= self._maxiter:
                raise RuntimeError('niter > maxiter')

    @property
    def prepared(self):
        return self._prepared

    @property
    def done(self):
        return self._done

    @property
    def finished(self):
        return self._finished

    @prepared.setter
    def prepared(self, value):
        self._prepared = value
        if not self.prepared:
            self.done = False

    @done.setter
    def done(self, value):
        self._done = value
        if not self.done:
            self._finished = False

    @finished.setter
    def finished(self, value):
        self._finished = value

    @abstractmethod
    def prepare(self):
        pass

    @abstractmethod
    def do(self):
        pass

    @abstractmethod
    def finish(self):
        pass


class HighLevelIntegrator(HighLevelInterface):

    def __init__(self):
        super(HighLevelIntegrator, self).__init__()

    def prepare(self):
        self.prepared = True
        print "Prepare:"
        print " Predict reflections"
        print " Extract shoeboxes"

    def do(self):
        self.done = True
        print "Do:"
        print " Calculate centroids"
        print " Calculate background"
        print " Calculate intensities"

    def finish(self):
        self.finished = True
        print "Finished:"
        print " Perform corrections"


integrator = HighLevelIntegrator()
integrator()
