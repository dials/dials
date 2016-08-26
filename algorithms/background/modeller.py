
from __future__ import division
from dials_algorithms_background_modeller_ext import *



class BackgroundModellerExecutor(object):

  def __init__(self):
    self.result = None

  def process(self, image_volume, experiments, reflections):
    from dials.algorithms.integration.processor import job
    from logging import info

    # Write some output
    info(" Background modelling; job: %d; frames: %d -> %d; # Reflections: %d" % (
      job.index,
      image_volume.frame0(),
      image_volume.frame1(),
      len(reflections)))

    # Compute the shoebox mask
    reflections.compute_mask(
      experiments  = experiments,
      image_volume = image_volume)

    # Compute the sum, sum^2 and the number of contributing pixels
    return MultiPanelBackgroundStatistics(image_volume)

  def accumulate(self, index, data):
    if self.result is None:
      self.result = data
    else:
      self.result += data

  def finalize_model(self):
    from logging import info
    info("Finalizing model")
    return self.result

#     mean = self.result.get(0).mean()
#     var = self.result.get(0).variance()
#     disp = self.result.get(0).dispersion()

#     return (mean, var, disp)


class BackgroundModeller(object):
  '''
  A class to help with background modelling

  '''

  def __init__(self,
               experiments,
               reflections,
               params):
    '''
    Initialize the modeller

    :param experiments: The experiment list
    :param reflections: The reflections to process
    :param params: The parameters to use

    '''
    # Check all reflections have same imageset and get it
    imageset = experiments[0].imageset
    for expr in experiments:
      assert expr.imageset == imageset, "All experiments must share and imageset"

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    self.model = None

  def compute(self):
    '''
    Integrate the data

    '''
    from dials.algorithms.integration.image_integrator import ProcessorImage
    from dials.util.command_line import heading
    from logging import info, debug

    # Init the report
    self.profile_model_report = None
    self.integration_report = None

    # Create summary format
    fmt = (
      ' Processing the following experiments:\n'
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      ' Imagesets:   %d\n'
    )

    # Print the summary
    info(fmt % (
      len(self.experiments),
      len(self.experiments.beams()),
      len(self.experiments.detectors()),
      len(self.experiments.goniometers()),
      len(self.experiments.scans()),
      len(self.experiments.crystals()),
      len(self.experiments.imagesets())))

    # Print a heading
    info("=" * 80)
    info("")
    info(heading("Modelling background"))
    info("")

    # Compute some reflection properties
    self.reflections.compute_zeta_multi(self.experiments)
    self.reflections.compute_d(self.experiments)
    self.reflections.compute_bbox(self.experiments)

    # Construvt the image integrator processor
    processor = ProcessorImage(
      self.experiments,
      self.reflections,
      self.params)
    processor.executor = BackgroundModellerExecutor()

    # Do the processing
    _, time_info = processor.process()

    # Compute the model
    self.model = processor.executor.finalize_model()

    # Print the time info
    info(str(time_info))
    info("")

    # Return the reflections
    return self.model
