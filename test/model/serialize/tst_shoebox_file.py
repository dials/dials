
from __future__ import division
from dials.model.serialize import ShoeboxFileExporter
from dials.model.serialize import ShoeboxFileImporter

class Test(object):

  def __init__(self):
    from dials.array_family import flex
    from random import randint

    self.filename = "shoebox.dat"
    self.nframes = 5
    self.npanels = 2
    self.width = 1000
    self.height = 1000
    self.nrefl = 1000
    self.panels = flex.size_t(self.nrefl)
    self.bboxes = flex.int6(self.nrefl)

    for i in range(self.nrefl):
      self.panels[i] = randint(0, self.npanels-1)
      x0 = randint(0, self.width-10)
      x1 = randint(x0+1, x0 + 10)
      y0 = randint(0, self.height-10)
      y1 = randint(y0+1, y0 + 10)
      z0 = randint(0, self.nframes-2)
      z1 = randint(z0+1, self.nframes)
      self.bboxes[i] = (x0, x1, y0, y1, z0, z1)

    self.images = []
    a = 0
    b = a + self.width * self.height
    for i in range(self.nframes):
        image = flex.int(range(a, b))
        image.reshape(flex.grid(self.height, self.width))
        self.images.append(image)
        a = b
        b = a + self.width * self.height

  def run(self):
    
    self.export_data()
    self.import_data()
    self.check_data()

  def export_data(self):
    from os.path import isfile

    # Create the exporter
    exporter = ShoeboxFileExporter(
      self.filename, 
      self.panels, 
      self.bboxes, 
      self.nframes, 
      self.npanels)

    # Add all the images
    for image in self.images:
      exporter.next(image)
      
    # Check the file exists
    assert(isfile(self.filename))

    # Writing test passed
    print 'OK'

  def import_data(self):
    
    from dials.array_family import flex

    # Create the lookup
    gain = tuple([flex.double(flex.grid(self.height, self.width), 1.0) 
                  for i in range(self.npanels)])

    dark = tuple([flex.double(flex.grid(self.height, self.width), 0.0) 
                  for i in range(self.npanels)])
    
    mask = tuple([flex.bool(flex.grid(self.height, self.width), True) 
                  for i in range(self.npanels)])
    
    # Create the importer
    importer = ShoeboxFileImporter(self.filename, gain, dark, mask)

    # Run some tests
    self.test_magic_number(importer)
    self.test_read_data(importer)
    self.test_data_values(importer)

  def test_magic_number(self, importer):
    pass

  def test_read_data(self, importer):
    pass

  def test_data_values(self, importer):
    pass




if __name__ == '__main__':
  test = Test()
  test.run()
