"""Simulate a rotation dataset with a smoothly-varying beam position for
refinement testing. Script based on tst_nanoBragg_basic.py"""

from __future__ import division
from scitbx.array_family import flex
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import miller
from scitbx import matrix
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.beam import Beam, BeamFactory
from dxtbx.model.scan import ScanFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model import Crystal
from dxtbx.model.experiment_list import Experiment
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.experiment_list import ExperimentListDumper

pdb_lines = """HEADER TEST
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""

class Simulation(object):

  def __init__(self, override_fdp=None):

    # Set up detector
    distance = 100
    pixel_size = 0.1
    image_size = (1000, 1000)
    beam_centre_mm = (pixel_size * image_size[0] / 2,
                      pixel_size * image_size[1] / 2)
    self.detector = DetectorFactory().simple(
          'CCD', distance, beam_centre_mm, '+x', '-y', (pixel_size, pixel_size),
          image_size)

    # Set up beam
    self.beam = BeamFactory().simple(wavelength=1)

    # Set up scan
    sweep_width = 90.0
    osc_start = 0.0
    image_width = 0.2
    oscillation = (osc_start, image_width)
    from math import ceil
    nframes = int(ceil(sweep_width / image_width))
    image_range = (1, nframes)
    exposure_times = 0.0
    epochs = [0] * nframes
    self.scan=ScanFactory().make_scan(image_range, exposure_times, oscillation,
        epochs, deg=True)

    # Set up goniometer
    self.goniometer = GoniometerFactory.known_axis(
        self.detector[0].get_fast_axis())

    # Set up simulated structure factors
    self.sfall = self.fcalc_from_pdb(resolution=1.6,algorithm="direct",
        override_fdp=override_fdp)

    # Set up crystal
    self.crystal = Crystal(real_space_a = (50,0,0),
                           real_space_b = (0,60,0),
                           real_space_c = (0,0,70),
                           space_group_symbol="P1")
    axis = matrix.col(elems=(-0.14480368275412925,
                             -0.6202131724405818,
                             -0.7709523423610766))
    self.crystal.set_U(axis.axis_and_angle_as_r3_rotation_matrix(
        angle=0.625126343998969))

    return

  def dump_experiments(self, file_name):
    e = Experiment(beam=self.beam, detector=self.detector, crystal=self.crystal,
        scan=self.scan, goniometer=self.goniometer)
    el = ExperimentList()
    el.append(e)
    dump = ExperimentListDumper(el)
    dump.as_json(file_name)
    return

  def fcalc_from_pdb(self, resolution, algorithm=None, override_fdp=None):
    from iotbx import pdb
    pdb_inp = pdb.input(source_info=None,lines = pdb_lines)
    xray_structure = pdb_inp.xray_structure_simple()
    wavelength=self.beam.get_wavelength()
    #
    # take a detour to calculate anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import sasaki, henke
      #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = override_fdp if override_fdp is not None else expected_henke.fdp()

    # how do we do bulk solvent?
    primitive_xray_structure = xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    fcalc = P1_primitive_xray_structure.structure_factors(
      d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
    return fcalc.amplitudes()

  def set_varying_beam(self, along="fast", npixels_drift=5):

    assert along in ['fast', 'slow', 'both']
    num_scan_points = self.scan.get_num_images() + 1
    s0 = matrix.col(self.beam.get_s0())
    beam_centre_px = self.detector[0].get_beam_centre_px(s0)
    if along == 'fast':
      start_beam_centre = (beam_centre_px[0] - npixels_drift/2,
                           beam_centre_px[1])
      end_beam_centre = (beam_centre_px[0] + npixels_drift/2,
                         beam_centre_px[1])
    elif along == 'slow':
      start_beam_centre = (beam_centre_px[0],
                           beam_centre_px[1] - npixels_drift/2)
      end_beam_centre = (beam_centre_px[0],
                         beam_centre_px[1] + npixels_drift/2)
    elif along == 'both':
      from math import sqrt
      offset = sqrt(2.0) * npixels_drift / 4.0
      start_beam_centre = (beam_centre_px[0] - offset,
                           beam_centre_px[1] - offset)
      end_beam_centre = (beam_centre_px[0] + offset,
                         beam_centre_px[1] + offset)

    start_lab = matrix.col(self.detector[0].get_pixel_lab_coord(
        start_beam_centre))
    end_lab = matrix.col(self.detector[0].get_pixel_lab_coord(end_beam_centre))
    axis = start_lab.cross(end_lab).normalize()
    full_angle = start_lab.angle(end_lab)
    angle_step = full_angle / self.scan.get_num_images()
    angles = [e * angle_step for e in range(num_scan_points)]

    start_s0 = start_lab.normalize() * s0.length()
    s0_list = [start_s0.rotate_around_origin(axis=axis, angle=e)
               for e in angles]

    self.beam.set_s0_at_scan_points(s0_list)

    return

  def generate_image(self, image_no=1):

    # Set the beam for this image
    beam = Beam(self.beam)
    if self.beam.get_num_scan_points() == self.scan.get_num_images() + 1:
      # get s0 for the midpoint of this image
      s0_start = matrix.col(self.beam.get_s0_at_scan_point(image_no - 1))
      s0_end = matrix.col(self.beam.get_s0_at_scan_point(image_no))
      s0_mid = (s0_start + s0_end).normalize() * s0_start.length()
      beam.set_s0(s0_mid)

    print "Generating image {0}".format(image_no)
    phi_deg = self.scan.get_angle_from_image_index(image_no, deg=True)
    print "Rotation angle={0} degrees".format(phi_deg)
    print "DIALS beam centre will be", self.detector[0].get_beam_centre(
      beam.get_s0())

    # Construct simulation
    SIM = nanoBragg(self.detector, beam, verbose=0)
    SIM.Ncells_abc=(10,10,10)
    # set different random number seed for noise generation for each image
    SIM.seed = image_no
    SIM.oversample=1
    SIM.progress_meter=False # only really useful for long runs
    # SIM.default_F=100 # this will become F000, marking the beam center
    # use crystal structure to initialize Fhkl array
    SIM.Fhkl=self.sfall

    # This does not 'stick': "ERROR: cannot initialize without a cell"
    #SIM.Amatrix = self.crystal.get_A()

    # WORKAROUND: Instead, use nanoBragg to set the A matrix by missets, then
    # update the dxtbx crystal to keep track
    SIM.missets_deg= (10,20,30)
    # Apparently this needs the transpose
    self.crystal.set_A(matrix.sqr(SIM.Amatrix).transpose())

    SIM.xtal_shape=shapetype.Tophat # fastest option, least realistic
    SIM.flux=1e12 # photons/s
    SIM.beamsize_mm=0.1 # assumes round beam
    SIM.exposure_s=0.1

    # Set rotation image parameters
    SIM.phi_deg = phi_deg # spindle phi rotation start angle (deg)
    SIM.osc_deg = self.scan.get_oscillation()[1] # phi rotation range (deg)
    SIM.phisteps = 100 # 1000 would be better, for those with patience

    # Set low detector calibration noise of 1% (instead of default 3%)
    SIM.detector_calibration_noise_pct=1

    # Now actually burn up some CPU
    SIM.add_nanoBragg_spots()

    # Simulated crystal is 1000 unit cells (50 nm wide). Amplify spot
    # signal to simulate physical crystal of 4000x larger: 200 um (64e9 x the
    # volume)
    SIM.raw_pixels *= 64e9

    # Write out the noise-free image with a pedestal matching that reported in
    # the header
    SIM.raw_pixels += SIM.adc_offset_adu
    fileout = "intimage_{0:03d}.img".format(image_no)
    SIM.to_smv_format(fileout=fileout, intfile_scale=1)
    SIM.raw_pixels -= SIM.adc_offset_adu

    # Add scatter from rough approximation to water: interpolation points for
    # sin(theta/lambda) vs structure factor
    bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),
        (0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),
        (0.345,4.36),(0.436,3.77),(0.5,3.17)])
    SIM.Fbg_vs_stol = bg
    SIM.amorphous_sample_thick_mm = 0.1
    SIM.amorphous_density_gcm3 = 1
    SIM.amorphous_molecular_weight_Da = 18
    SIM.add_background()

    # Add scatter from rough approximation to air
    bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),
        (0.5,4.22)])
    SIM.Fbg_vs_stol = bg
    SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
    SIM.amorphous_density_gcm3 = 1.2e-3
    SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
    SIM.add_background()

    # detector PSF params for apply_psf(), called within add_noise()
    SIM.detector_psf_kernel_radius_pixels=5;
    SIM.detector_psf_fwhm_mm=0.08;
    SIM.detector_psf_type=shapetype.Fiber

    # Add the various noise sources
    SIM.add_noise()

    # Write out image with noise
    fileout = "noiseimage_{0:03d}.img".format(image_no)
    SIM.to_smv_format(fileout=fileout, intfile_scale=1)

    SIM.free_all()

  def generate_all_images(self):
    for i in range(self.scan.get_num_images()):
      self.generate_image(i + 1)

def check_anomalous_sign(mtz_file_name, override_fdp=None):
  '''Test that the anomalous differences in a processed dataset are the right
  way round, so we haven't entered the 'inverted world' during the simulation.
  Using an MTZ from Aimless for this rather than data straight out of
  dials.integrate, mainly to use the simple any_reflection_file reader and to
  avoid having to decide between I and IPR etc.'''

  # Extract miller array of intensities from the MTZ
  from iotbx.reflection_file_reader import any_reflection_file
  hkl_in = any_reflection_file(file_name=mtz_file_name)
  mtz_miller_arrays = hkl_in.as_miller_arrays()
  # assume I(+),SIGI(+),I(-),SIGI(-) is in the 2nd array
  mtz_intens = mtz_miller_arrays[1]

  assert mtz_intens.is_xray_intensity_array()
  mtz_ampl = mtz_intens.f_sq_as_f()
  print "MTZ anomalous signal = {0:.2f}%".format(
      mtz_ampl.anomalous_signal()*100)
  mtz_anom_diff = mtz_ampl.anomalous_differences()

  # Take the largest 10% of anomalous differences
  #abs_diff = flex.abs(mtz_anom_diff.data())
  #cutoff = flex.sorted(abs_diff)[int(len(abs_diff) * 0.9)]
  #mtz_anom_diff = mtz_anom_diff.select(abs_diff > cutoff)

  # Simulation miller array of anomalous differences
  sim=Simulation(override_fdp=override_fdp)
  print "Simulation anomalous signal = {0:.2f}%".format(
      sim.sfall.anomalous_signal()*100)
  sim_anom_diff = sim.sfall.anomalous_differences()

  # CCs of intensities and amplitudes
  sim_intens = sim.sfall.f_as_f_sq()
  mtz_intens_comm, sim_intens_comm = mtz_intens.common_sets(other=sim_intens)
  cc = flex.linear_correlation(
      mtz_intens_comm.data(), sim_intens_comm.data()).coefficient()
  print "CC between intensities from MTZ and simulation = {0:.2f}%".format(
      cc * 100)
  mtz_ampl_comm, sim_ampl_comm = mtz_ampl.common_sets(other=sim.sfall)
  cc = flex.linear_correlation(
      mtz_ampl_comm.data(), sim_ampl_comm.data()).coefficient()
  print "CC between amplitudes from MTZ and simulation = {0:.2f}%".format(
      cc * 100)

  # Write ideal SFs to an MTZ for comparison with viewhkl
  mtz_dataset = sim_intens_comm.as_mtz_dataset(column_root_label="I")
  mtz_dataset.add_miller_array(sim_ampl_comm, column_root_label="F")
  mtz_dataset.mtz_object().write("sim_data.mtz")

  # Get common sets of anomalous differences
  mtz_dano_comm, sim_dano_comm = mtz_anom_diff.common_sets(other=sim_anom_diff)

  cc = flex.linear_correlation(
      mtz_dano_comm.data(), sim_dano_comm.data()).coefficient()
  msg = ("CC between anomalous differences of common "
         "reflections = {0:.2f}%").format(cc *100)
  print msg

  # Compare signs
  mtz_pos = mtz_dano_comm.data() >= 0.0
  sim_pos = sim_dano_comm.data() >= 0.0
  n = len(mtz_pos)
  nsame = (mtz_pos==sim_pos).count(True)
  print "Anomalous differences with the same sign: {0} ({1:0.2f}%)".format(
    nsame, 100.*nsame/n)

def process_data(prefix='int'):
  assert prefix in ['int', 'noise']
  from libtbx import easy_run
  cmd = 'dials.import {0}image_*.img'.format(prefix)
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.find_spots datablock.json nproc=8'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.index datablock.json strong.pickle'
  result = easy_run.fully_buffered(command=cmd)
  # now patch experiments.json to ensure the crystal model is indexed with the
  # right basis
  el = ExperimentListFactory.from_json_file('experiments.json')
  ref = ExperimentListFactory.from_json_file('sim_experiments.json')
  el[0].crystal = ref[0].crystal
  dump = ExperimentListDumper(el)
  dump.as_json('mod_experiments.json')
  cmd = 'dials.index mod_experiments.json strong.pickle'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.refine experiments.json indexed.pickle'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.refine refined_experiments.json refined.pickle'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.integrate refined_experiments.json refined.pickle nproc=8'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'dials.export integrated_experiments.json integrated.pickle'
  result = easy_run.fully_buffered(command=cmd)
  cmd = 'aimless hklin integrated.mtz hklout {0}_scaled.mtz'.format(prefix)
  stdin = ['resolution 2.0', 'anomalous on']
  result = easy_run.fully_buffered(command=cmd, stdin_lines=stdin)
  with open('{0}_aimless.log'.format(prefix), 'w') as f:
    f.writelines(result.stdout_lines)

def run_no_drift():
  sim=Simulation()
  sim.generate_all_images()
  #sim.generate_image(1)
  sim.dump_experiments('sim_experiments.json')
  process_data(prefix='int')
  process_data(prefix='noise')
  print "Checking processed data from noise-free images"
  check_anomalous_sign("int_scaled.mtz")
  print
  print "Checking processed data from noisy images"
  check_anomalous_sign("noise_scaled.mtz")

def run_drift_along_fast():
  sim=Simulation()
  sim.set_varying_beam()
  sim.generate_all_images()

def run_drift_along_slow():
  sim=Simulation()
  sim.set_varying_beam(along='slow')
  sim.generate_all_images()

def run_drift_along_diag():
  sim=Simulation()
  sim.set_varying_beam(along='both')
  sim.generate_all_images()

if __name__=="__main__":

  run_no_drift()
  #run_drift_along_fast()

