from dials.framework import toplevel

class Integrator(toplevel.Interface):

  def __init__(self, sweep, crystal, params):

    from dials.framework.plugin import Registry
    from dials.interfaces.standard_interfaces import IntegratorBoundingBox
    from dials.interfaces.standard_interfaces import IntegratorBackground
    from dials.interfaces.standard_interfaces import IntegratorCentroid
    from dials.interfaces.standard_interfaces import IntegratorIntensity

    # Save some stuff
    self.sweep = sweep
    self.crystal = crystal

    # Configure the bounding box algorithm
    bbox_algorithm = Registry.factory(IntegratorBoundingBox)(
        params.integrator.bounding_box.algorithm,
        sweep, crystal, params)

    # Configure the background algorithm
    self.compute_background = Registry.factory(IntegratorBackground)(
        params.integrator.background.algorithm,
        sweep, crystal, params)

    # Configure the centroid algorithm
    self.compute_centroid = Registry.factory(IntegratorCentroid)(
        params.integrator.centroid.algorithm,
        sweep, crystal, params)

    # Configure the intensity algorithm
    self.compute_intensity = Registry.factory(IntegratorIntensity)(
        params.integrator.intensity.algorithm,
        sweep, crystal, params)

  def prepare(self):
    self.predictions = self.predict()
    self.shoeboxes = self.extract_shoeboxes(self.predictions)
    self.prepared = True

  def process(self):

    self.compute_background(self.shoeboxes, self.predictions)
    self.observations = flex.observation(
        self.shoeboxes.panel(),
        self.compute_centroid(self.shoeboxes, self.predictions),
        self.compute_intensity(self.shoeboxes, self.predictions))

    self.reflections.compute_background()
    self.reflections.compute_centroid()
    self.reflections.compute_intensity()

    self.processed = True

  def finish(self):
    self.finished = True


class IntegratorFactory(object):

  @staticmethod
  def create(sweep, crystal, params):
    return IntegratorInterface(sweep, crystal, params)
