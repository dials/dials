
from __future__ import division

def read_experiments_and_reflections():
  from collections import namedtuple
  from dials.array_family import flex
  import cPickle as pickle
  from dxtbx.model.experiment_list import ExperimentListFactory
  import libtbx.load_env
  from os.path import join
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'SKIP: dials_regression not configured'
    exit(0)

  directory = join(dials_regression, "integration_test_data", "shoeboxes")
  experiments_filename = join(directory, "integrated_experiments.json")
  reflections_filename = join(directory, "shoeboxes_0_0.pickle")
  reference_filename = join(directory, "reference_profiles.pickle")

  experiments = ExperimentListFactory.from_json_file(experiments_filename)
  reflections = flex.reflection_table.from_pickle(reflections_filename)
  reference = pickle.load(open(reference_filename))

  Data = namedtuple("Data", ["experiments", "reflections", "reference"])
  return Data(
    experiments=experiments,
    reflections=reflections,
    reference=reference)


data = read_experiments_and_reflections()


def tst_gaussianrs_mask_calculator():

  from dials.algorithms.integration.parallel_integrator import MaskCalculator
  from dials.algorithms.integration.parallel_integrator import MultiCrystalMaskCalculator

  algorithm = MultiCrystalMaskCalculator()
  for e in data.experiments:
    alg = MaskCalculator(
      e.beam,
      e.detector,
      e.goniometer,
      e.scan,
      e.profile.delta_b(deg=False),
      e.profile.delta_m(deg=False))
    algorithm.append(alg)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  for r in reflections:
    algorithm(r, False)

  print 'OK'


def tst_simple_background_calculator():

  from dials.algorithms.integration.parallel_integrator import SimpleBackgroundCalculatorFactory

  algorithm = SimpleBackgroundCalculatorFactory.create(data.experiments)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  count = 0
  for r in reflections:
    try:
      algorithm(r)
    except Exception:
      count += 1
  assert len(reflections) == 15193, len(reflections)
  assert count == 333, count

  print 'OK'


def tst_glm_background_calculator():

  from dials.algorithms.integration.parallel_integrator import GLMBackgroundCalculatorFactory

  algorithm = GLMBackgroundCalculatorFactory.create(data.experiments)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  count = 0
  for r in reflections:
    try:
      algorithm(r)
    except Exception:
      count += 1
  assert len(reflections) == 15193, len(reflections)
  assert count == 333, count

  print 'OK'


def tst_gmodel_background_calculator():
  pass


class IntensityCalculatorFactory(object):

  @classmethod
  def create(Class,
             detector_space=False,
             deconvolution=False):
    from dials.algorithms.integration.parallel_integrator import GaussianRSIntensityCalculatorFactory
    from dials.algorithms.integration.parallel_integrator import GaussianRSReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import GaussianRSMultiCrystalReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import  ReferenceProfileData
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.array_family import flex
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

    reference  = data.reference[0]
    experiments = data.experiments

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9

    data_spec = GaussianRSMultiCrystalReferenceProfileData()
    for e in experiments:

      sampler = CircleSampler(
        e.detector[0].get_image_size(),
        e.scan.get_array_range(),
        num_scan_points)


      spec = TransformSpec(
        e.beam,
        e.detector,
        e.goniometer,
        e.scan,
        e.profile.sigma_b(deg=False),
        e.profile.sigma_m(deg=False),
        e.profile.n_sigma() * 1.5,
        5)

      temp = reference

      reference = ReferenceProfileData()
      for d, m in temp:
        reference.append(d, m)

      spec = GaussianRSReferenceProfileData(reference, sampler, spec)

      data_spec.append(spec)

    return GaussianRSIntensityCalculatorFactory.create(
      data_spec,
      detector_space,
      deconvolution)


def tst_gaussianrs_reciprocal_space_intensity_calculator():

  algorithm = IntensityCalculatorFactory.create(
    detector_space = False,
    deconvolution = False)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  count = 0
  nint = 0
  for r in reflections:
    try:
      algorithm(r, [])
    except Exception, e:
      count += 1

  assert len(reflections) == 15193, len(reflections)
  assert count == 5295, count

  print 'OK'



def tst_gaussianrs_detector_space_intensity_calculator():

  algorithm = IntensityCalculatorFactory.create(
    detector_space = True,
    deconvolution = False)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  count = 0
  nint = 0
  for r in reflections:
    try:
      algorithm(r, [])
      partiality_old = r.get('partiality_old')
      partiality_new = r.get('partiality')
    except Exception, e:
      count += 1
      continue

    assert partiality_old < 1.0 and partiality_old >= 0, partiality_old
    assert partiality_new < 1.0 and partiality_new >= 0, partiality_new


  assert len(reflections) == 15193, len(reflections)
  assert count == 4801, count


  print 'OK'


def tst_gaussianrs_detector_space_with_deconvolution_intensity_calculator():

  algorithm = IntensityCalculatorFactory.create(
    detector_space = True,
    deconvolution = True)

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  count = 0
  nint = 0
  for r in reflections:
    try:
      algorithm(r, [])
      partiality_old = r.get('partiality_old')
      partiality_new = r.get('partiality')
    except Exception, e:
      count += 1
      continue

    assert partiality_old < 1.0 and partiality_old >= 0, partiality_old
    assert partiality_new < 1.0 and partiality_new >= 0, partiality_new


  assert len(reflections) == 15193, len(reflections)
  assert count == 4801, count

  print 'OK'

def tst_gaussianrs_detector_space_with_deconvolution_intensity_calculator2():
  from scitbx import matrix

  from dials.array_family import flex
  reflections = flex.reflection_table_to_list_of_reflections(data.reflections)

  R = None
  for r in reflections:
    if r.get('partiality') > 0.9 and r.get("intensity.sum.value") > 100:
      R = r
      break
  assert R is not None

  s1 = R.get("s1")
  px = R.get("xyzcal.px")
  mm = R.get("xyzcal.mm")
  sbox = R.get("shoebox")

  px1 = (px[0]-3, px[1]-3, px[2])
  px2 = (px[0]+3, px[1]+3, px[2])
  mm1 = (mm[0]-3*0.172, mm[1]-3*0.172, mm[2])
  mm2 = (mm[0]+3*0.172, mm[1]+3*0.172, mm[2])

  s11 = matrix.col(data.experiments[0].detector[0].get_pixel_lab_coord(px1[0:2])).normalize()
  s12 = matrix.col(data.experiments[0].detector[0].get_pixel_lab_coord(px2[0:2])).normalize()

  R1 = R.copy()
  R2 = R.copy()

  R1.set_vec3_double("xyzcal.px", px1)
  R1.set_vec3_double("xyzcal.mm", mm1)
  R1.set_vec3_double("s1", s11)

  R2.set_vec3_double("xyzcal.px", px2)
  R2.set_vec3_double("xyzcal.mm", mm2)
  R2.set_vec3_double("s1", s12)

  compute_intensity= IntensityCalculatorFactory.create(
    detector_space = True,
    deconvolution = False)

  compute_intensity(R, [])
  R.set_double("intensity.prf.value.old", R.get("intensity.prf.value"))
  R.set_double("intensity.prf.variance.old", R.get("intensity.prf.variance"))

  R_no_deconvolution = R.copy()
  #print "Partiality", R.get("partiality")
  #print "Partiality Old", R.get("partiality_old")
  #print "Intensity", R.get("intensity.prf.value")
  #print "Intensity Old", R.get("intensity.prf.value.old")


  from dials.algorithms.integration.parallel_integrator import MaskCalculator
  from dials.algorithms.integration.parallel_integrator import MultiCrystalMaskCalculator

  mask_calculator = MultiCrystalMaskCalculator()
  for e in data.experiments:
    alg = MaskCalculator(
      e.beam,
      e.detector,
      e.goniometer,
      e.scan,
      e.profile.delta_b(deg=False),
      e.profile.delta_m(deg=False))
    mask_calculator.append(alg)

  mask_calculator(R1, True)
  mask_calculator(R2, True)

  #print R.get("shoebox").mask.as_numpy_array()[0,:,:]

  compute_intensity= IntensityCalculatorFactory.create(
    detector_space = True,
    deconvolution = True)

  compute_intensity(R, [R1, R2])
  #compute_intensity(R, [R1, R2])

  R_deconvolution = R.copy()

  partiality1 = R_no_deconvolution.get("partiality")
  partiality2 = R_deconvolution.get("partiality")
  intensity = R_deconvolution.get("intensity.prf.value")
  variance = R_deconvolution.get("intensity.prf.variance")

  assert partiality1 <= partiality2, (partiality1, partiality2)

  assert abs(intensity - 179.04238328) < 1e-7, intensity
  assert abs(variance - 181.396028623) < 1e-7, variance

  #print "Partiality", R.get("partiality")
  #print "Partiality Old", R.get("partiality_old")
  #print "Intensity", R.get("intensity.prf.value")
  #print "Intensity Old", R.get("intensity.prf.value.old")
  #print "Variance", R.get("intensity.prf.variance")
  #print "Variance Old", R.get("intensity.prf.variance.old")

  #print R.get("shoebox").mask.all()
  #sbox = R.get("shoebox")
  #print sbox.mask.count(5) + sbox.mask.count(37) + sbox.mask.count(51)
  print 'OK'


def tst_gaussianrs_profile_data_pickling():
    from dials.algorithms.integration.parallel_integrator import GaussianRSReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import GaussianRSMultiCrystalReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import  ReferenceProfileData
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.array_family import flex
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

    reference  = data.reference[0]
    experiments = data.experiments

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9

    data_spec = GaussianRSMultiCrystalReferenceProfileData()
    for e in experiments:

      sampler = CircleSampler(
        e.detector[0].get_image_size(),
        e.scan.get_array_range(),
        num_scan_points)


      spec = TransformSpec(
        e.beam,
        e.detector,
        e.goniometer,
        e.scan,
        e.profile.sigma_b(deg=False),
        e.profile.sigma_m(deg=False),
        e.profile.n_sigma() * 1.5,
        5)

      temp = reference

      reference = ReferenceProfileData()
      for d, m in temp:
        reference.append(d, m)

      spec = GaussianRSReferenceProfileData(reference, sampler, spec)

      data_spec.append(spec)

    import cPickle as pickle

    s = pickle.dumps(data_spec)

    data_spec2 = pickle.loads(s)

    print 'OK'


def tst_gaussianrs_intensity_calculator():
  tst_gaussianrs_reciprocal_space_intensity_calculator()
  tst_gaussianrs_detector_space_intensity_calculator()
  tst_gaussianrs_detector_space_with_deconvolution_intensity_calculator()
  tst_gaussianrs_detector_space_with_deconvolution_intensity_calculator2()
  tst_gaussianrs_profile_data_pickling()



if __name__ == '__main__':

  tst_gaussianrs_mask_calculator()
  tst_simple_background_calculator()
  tst_glm_background_calculator()
  tst_gmodel_background_calculator()
  tst_gaussianrs_intensity_calculator()
