
from __future__ import division
from dials_algorithms_background_modeller_ext import *


class FinalizeModel(object):

  def __init__(self, experiments, kernel_size=10, niter=100):
    from dials.algorithms.background.gmodel import PolarTransform

    # Set some parameters
    self.kernel_size = (10, 1)
    self.niter = niter

    # Check the input
    assert len(experiments) == 1
    experiment = experiments[0]
    assert len(experiment.detector) == 1

    # Create the transform object
    self.transform = PolarTransform(
      experiment.beam,
      experiment.detector[0],
      experiment.goniometer)

  def finalize(self, data, mask):
    from logging import info
    from dials.algorithms.image.filter import median_filter
    from dials.algorithms.image.fill_holes import diffusion_fill
    from dials.algorithms.image.fill_holes import simple_fill
    from dials.array_family import flex

    # Print some image properties
    sub_data = data.as_1d().select(mask.as_1d())
    info('Raw image statistics:')
    info('  min:  %d' % int(flex.min(sub_data)))
    info('  max:  %d' % int(flex.max(sub_data)))
    info('  mean: %d' % int(flex.mean(sub_data)))

    # Transform to polar
    info('Transforming image data to polar grid')
    result = self.transform.to_polar(data, mask)
    data = result.data()
    mask = result.mask()
    sub_data = data.as_1d().select(mask.as_1d())
    info('Polar image statistics:')
    info('  min:  %d' % int(flex.min(sub_data)))
    info('  max:  %d' % int(flex.max(sub_data)))
    info('  mean: %d' % int(flex.mean(sub_data)))

    from matplotlib import pylab
    pylab.imshow(mask.as_numpy_array())
    pylab.show()
    pylab.imshow(data.as_numpy_array(), vmax=200)
    pylab.show()


    # Filter the image to remove noise
    info('Applying median filter')
    data = median_filter(data, mask, self.kernel_size)
    sub_data = data.as_1d().select(mask.as_1d())
    info('Median polar image statistics:')
    info('  min:  %d' % int(flex.min(sub_data)))
    info('  max:  %d' % int(flex.max(sub_data)))
    info('  mean: %d' % int(flex.mean(sub_data)))

    # Fill any remaining holes
    info("Filling holes")
    data = simple_fill(data, mask)
    data = diffusion_fill(data, mask, self.niter)
    mask = flex.bool(data.accessor(), True)
    sub_data = data.as_1d().select(mask.as_1d())
    info('Filled polar image statistics:')
    info('  min:  %d' % int(flex.min(sub_data)))
    info('  max:  %d' % int(flex.max(sub_data)))
    info('  mean: %d' % int(flex.mean(sub_data)))

    # Transform back
    info('Transforming image data from polar grid')
    result = self.transform.from_polar(data, mask)
    data = result.data()
    mask = result.mask()
    sub_data = data.as_1d().select(mask.as_1d())
    info('Final image statistics:')
    info('  min:  %d' % int(flex.min(sub_data)))
    info('  max:  %d' % int(flex.max(sub_data)))
    info('  mean: %d' % int(flex.mean(sub_data)))

    from matplotlib import pylab
    pylab.imshow(data.as_numpy_array())
    pylab.show()
    # Return the result
    return data



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
