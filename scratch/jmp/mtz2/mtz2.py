


from __future__ import division
import h5py


schema_url = 'https://github.com/nexusformat/definitions/blob/master/contributed_definitions/NXmx.nxdl.xml'


class NXobject(object):

  def __init__(self, handle):
    self._handle = handle

  @classmethod
  def factory(Class, handle, index=None):
    if index is None:
      name = Class.name
    else:
      name = '%s%d' % (Class.name, index)
    if Class.name in handle:
      obj = handle[name]
      assert(obj.attrs['NX_class'] == 'NX%s' % Class.name)
    else:
      obj = handle.create_group(name)
      obj.attrs['NX_class'] = 'NX%s' % Class.name
    return Class(obj)

  def __setitem__(self, name, value):
    self._handle.create_dataset(name, data=value)

  def __getitem__(self, name):
    return self._handle[name]


class NXobject_list(object):

  def __init__(self, handle):
    self._handle = handle

  def __getitem__(self, index):
    return self.Class.factory(self._handle, index=index)



class NXcollection(NXobject):

  name = 'collection'

class NXtransformations(NXobject):

  name = 'transformations'

class NXattenuator(NXobject):

  name = 'attenuator'

class NXdetector_module(NXobject):

  name = 'detector_module'

class NXdetector(NXobject):

  name = 'detector'

  @property
  def transformations(self):
    return NXtransformations.factory(self._handle)

  @property
  def collection(self):
    return NXcollection.factory(self._handle)

  @property
  def detector_module(self):
    return NXdetector_module.factory(self._handle)

class NXinstrument(NXobject):

  name = 'instrument'

  @property
  def attenuator(self):
    return NXattenuator.factory(self._handle)

  @property
  def detector(self):
    return NXdetector.factory(self._handle)

class NXbeam(NXobject):

  name = 'beam'

  @property
  def transformations(self):
    return NXtransformations.factory(self._handle)

class NXsample(NXobject):

  name = 'sample'

  @property
  def beam(self):
    return NXbeam.factory(self._handle)

class NXnote(NXobject):

  name = 'note'

class NXprocess(NXobject):

  name = 'process'


class NXprocess_list(NXobject_list):

  Class = NXprocess


class NXuser(NXobject):

  name = 'user'

class NXdata(NXobject):

  name = 'data'

class NXexperiment(NXobject):

  name = 'experiment'

class NXexperiment_list(NXobject_list):

  Class = NXexperiment


class NXdiffraction(NXobject):

  name = 'diffraction'

  @property
  def experiments(self):
    return NXexperiment_list(self._handle)



class NXentry(NXobject):

  name = 'entry'

  @property
  def data(self):
    return NXdata.factory(self._handle)

  @property
  def instrument(self):
    return NXinstrument.factory(self._handle)

  @property
  def sample(self):
    return NXsample.factory(self._handle)

  @property
  def process(self):
    return NXprocess_list(self._handle)

  @property
  def user(self):
    return NXuser.factory(self._handle)

  @property
  def diffraction(self):
    return NXdiffraction.factory(self._handle)

  @classmethod
  def factory(Class, handle):
    if 'entry' in handle:
      entry = handle['entry']
      assert(entry.attrs['NX_class'] == 'NXentry')
      assert(entry['definition'].value == 'NXmx')
      assert(entry['definition'].attrs['version'] >= 1)
    else:
      entry = handle.create_group('entry')
      entry.attrs['NX_class'] = 'NXentry'
      definition = entry.create_dataset('definition', data='NXmx')
      definition.attrs['version'] = 1
      definition.attrs['URL'] = schema_url
    return Class(entry)

class NXmx(object):

  def __init__(self, handle):
    self._handle = handle

  @property
  def entry(self):
    return NXentry.factory(self._handle)


class File(NXmx):

  def __init__(self, filename, mode='a'):
    super(File, self).__init__(h5py.File(filename, mode))

  def flush(self):
    self._handle.flush()

  def close(self):
    self._handle.close()



def test_writing():

  from time import strftime

  # Open the file
  outfile = File('test_file.mtz2', 'w')

  # Get the entry
  entry = outfile.entry

  # Start time and end time
  entry['start_time'] = strftime('%Y-%m-%dT%H:%M:%S')
  entry['end_time'] = strftime('%Y-%m-%dT%H:%M:%S')

  # Program info
  entry['program_name'] = 'dials.do_something'
  entry['program_name'].attrs['version'] = 1
  entry['program_name'].attrs['configuration'] = 'parameter=something'

  # Set some user data
  user = entry.user
  user['name'] = 'James Parkhurst'
  user['email'] = 'james.parkhurst@diamond.ac.uk'

  # Set some processing information (each program should add itself)
  process = entry.process[0]
  process['program'] = 'dials.first_program'
  process['version'] = 1
  process['date'] = strftime('%Y-%m-%dT%H:%M:%S')

  process = entry.process[1]
  process['program'] = 'dials.second_program'
  process['version'] = 1
  process['date'] = strftime('%Y-%m-%dT%H:%M:%S')

  # Flush the file
  outfile.flush()


def test_export_dials():

  from dials.array_family import flex
  print 'Creating dummy reflection table...'
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index(100)
  table['id'] = flex.size_t(100)
  table['intensity.sum.value'] = flex.double(100)
  table['intensity.sum.variance'] = flex.double(100)
  table['intensity.prf.value'] = flex.double(100)
  table['intensity.prf.variance'] = flex.double(100)
  table['lp'] = flex.double(100)
  table['panel'] = flex.size_t(100)
  table['bbox'] = flex.int6(100)
  table['xyzcal.px'] = flex.vec3_double(100)
  table['xyzcal.mm'] = flex.vec3_double(100)
  table['xyzobs.px.value'] = flex.vec3_double(100)
  table['xyzobs.px.variance'] = flex.vec3_double(100)
  table['xyzobs.mm.value'] = flex.vec3_double(100)
  table['xyzobs.mm.variance'] = flex.vec3_double(100)
  table['partiality'] = flex.double(100)
  table['d'] = flex.double(100)
  table['s1'] = flex.vec3_double(100)
  table['rlp'] = flex.vec3_double(100)
  table['background.mean'] = flex.double(100)
  table['entering'] = flex.bool(100)
  table['flags'] = flex.size_t(100)
  table['profile.correlation'] = flex.double(100)

  # Open the file
  outfile = File('test_file.mtz2', 'w')

  # Get the entry
  entry = outfile.entry

  print 'Writing reflection table stuff...'
  # Save some processed data
  diffraction = entry.diffraction

  # Set the experiments
  experiment = diffraction.experiments[0]
  experiment['beam'] = 0
  experiment['detector'] = 0
  experiment['goniometer'] = 0
  experiment['scan'] = 0
  experiment['crystal'] =0

  # Get columns into array
  col1, col2, col3 = zip(*list(table['miller_index']))
  col4 = table['id']
  col5 = table['intensity.sum.value']
  col6 = table['intensity.sum.variance']
  col7 = table['intensity.prf.value']
  col8 = table['intensity.prf.variance']
  col9 = table['lp']
  col10 = table['panel']
  col11, col12, col13, col14, col15, col16 = table['bbox'].parts()
  col17, col18, col19 = table['xyzcal.px'].parts()
  col20, col21, col22 = table['xyzcal.mm'].parts()
  col23, col24, col25 = table['xyzobs.px.value'].parts()
  col26, col27, col28 = table['xyzobs.px.variance'].parts()
  col29, col30, col31 = table['xyzobs.mm.value'].parts()
  col32, col33, col34 = table['xyzobs.mm.variance'].parts()
  col35 = table['partiality']
  col36 = table['d']
  col37 = table['background.mean']
  col38 = table['entering']
  col39 = table['flags']
  col40 = table['profile.correlation']

  # Some data
  diffraction['h'] = col1
  diffraction['k'] = col2
  diffraction['l'] = col3
  diffraction['id'] = col4
  diffraction['intensity_sum_value'] = col5
  diffraction['intensity_sum_variance'] = col6
  diffraction['intensity_prf_value'] = col7
  diffraction['intensity_prf_variance'] = col8
  diffraction['lp'] = col9
  diffraction['detector_module'] = col10
  diffraction['x0'] = col11
  diffraction['x1'] = col12
  diffraction['y0'] = col13
  diffraction['y1'] = col14
  diffraction['z0'] = col15
  diffraction['z1'] = col16
  diffraction['predicted_px_x'] = col17
  diffraction['predicted_px_y'] = col18
  diffraction['predicted_frame'] = col19
  diffraction['predicted_mm_x'] = col20
  diffraction['predicted_mm_y'] = col21
  diffraction['predicted_phi'] = col22
  diffraction['observed_px_x_value'] = col23
  diffraction['observed_px_x_variance'] = col24
  diffraction['observed_px_y_value'] = col25
  diffraction['observed_px_y_variance'] = col26
  diffraction['observed_frame_value'] = col27
  diffraction['observed_frame_variance'] = col28
  diffraction['observed_mm_x_value'] = col29
  diffraction['observed_mm_x_variance'] = col30
  diffraction['observed_mm_y_value'] = col31
  diffraction['observed_mm_y_variance'] = col32
  diffraction['observed_phi_value'] = col33
  diffraction['observed_phi_variance'] = col34
  diffraction['partiality'] = col35
  diffraction['d'] = col36
  diffraction['background.mean'] = col37
  diffraction['entering'] = col38
  diffraction['flags'] = col39
  diffraction['profile_correlation'] = col40

  # Flush the file
  outfile.flush()

if __name__ == '__main__':

  test_export_dials()
