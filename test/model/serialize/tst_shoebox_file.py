
from __future__ import division
from dials.model.serialize import ShoeboxFileExporter
from dials.model.serialize import ShoeboxFileImporter

class Test(object):

  def __init__(self):
    from dials.array_family import flex
    from random import randint

    self.filename = "shoebox.dat"
    self.nframes = 5
    self.npanels = 1
    self.width = 1000
    self.height = 1000
    self.nrefl = 1000
    self.panels = flex.size_t(self.nrefl)
    self.bboxes = flex.int6(self.nrefl)
    self.size = flex.size_t(self.nrefl)

    self.images = []
    a = 0
    b = a + self.width * self.height
    for i in range(self.nframes * self.npanels):
        image = flex.int_range(a, b)
        image.reshape(flex.grid(self.height, self.width))
        self.images.append(image)
        a = b
        b = a + self.width * self.height

    self.expected = []
    for i in range(self.nrefl):
      self.panels[i] = randint(0, self.npanels-1)
      x0 = randint(0, self.width-10)
      x1 = randint(x0+1, x0 + 10)
      y0 = randint(0, self.height-10)
      y1 = randint(y0+1, y0 + 10)
      z0 = randint(0, self.nframes-2)
      z1 = randint(z0+1, self.nframes)
      self.bboxes[i] = (x0, x1, y0, y1, z0, z1)
      self.size[i] = (x1 - x0) * (y1 - y0) * (z1 - z0)
      data = flex.int(flex.grid(z1-z0, y1-y0, x1-x0))
      for z in range(z0, z1):
        for y in range(y0, y1):
          for x in range(x0, x1):
            data[z-z0,y-y0,x-x0] = self.images[self.panels[i] + z * self.npanels][y,x]
      self.expected.append(data)

  def run(self):
    
    self.export_data()
    self.import_data()

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
    i = 0
    for frame in range(self.nframes):
      for panel in range(self.npanels):
        f, p  = exporter.next(self.images[i])
        assert(f == frame)
        assert(p == panel)
        i += 1

    # Check the file exists
    assert(isfile(self.filename))
    assert(exporter.finished())

    # Writing test passed
    print 'OK'

  def import_data(self):
    
    from dials.array_family import flex

    # Test the magic numbers
    self.test_contents()

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
    self.test_read_data(importer)

  def test_contents(self):
   
    import struct

    # Open the file
    infile = open(self.filename, "rb")

    # Read the magic numbers
    magic = struct.unpack('@I', infile.read(4))[0]
    version = struct.unpack('@I', infile.read(4))[0]
    assert(magic == 58008)
    assert(version == 1)

    # Read the header
    begin = struct.unpack('@I', infile.read(4))[0]
    assert(begin == 2)
    nrefl = struct.unpack('@I', infile.read(4))[0]
    assert(nrefl == self.nrefl)
   
    # Check the panels
    for i in range(nrefl):
      p = struct.unpack('@I', infile.read(4))[0]
      assert(p == self.panels[i])

    # Check the bboxes
    for i in range(nrefl):
      for j in range(6):
        b = struct.unpack('@i', infile.read(4))[0]
        assert(b == self.bboxes[i][j])

    # Read the end of the header
    end = struct.unpack('@I', infile.read(4))[0]
    assert(end == 3)

    # Test the start of the data
    begin = struct.unpack('@I', infile.read(4))[0]
    assert(begin == 4)
    
    # Check all the shoeboxes are laid out properly
    for i in range(nrefl):
      shoebox = struct.unpack('@I', infile.read(4))[0]
      assert(shoebox == 6)
      infile.seek(self.size[i] * 4, 1)

    # Check end of data
    end = struct.unpack('@I', infile.read(4))[0]
    assert(end == 5)

    # Test passed
    print 'OK'

  def test_read_data(self, importer):
      
    from dials.array_family import flex    

    # Check size of importer
    assert(len(importer) == self.nrefl)

    # Read by index
    for i in range(len(importer)):
      sbox = importer[i]
      self.check(sbox, i)

    # Read by range
    sbox = importer.select(0, len(importer))
    assert(len(sbox) == len(importer))
    self.check_range(sbox, 0, len(importer))

    # Read by indices
    indices = flex.random_size_t(100, self.nrefl)
    sbox = importer.select(indices)
    assert(len(sbox) == 100)
    self.check_many(sbox, indices)

    # Read by iterator
    for i, sbox in enumerate(importer):
      self.check(sbox, i)

    # Test passed
    print 'OK'

  def check(self, sbox, i):
    sbox2 = self.expected[i]
    assert(len(sbox.data) == len(sbox2))
    assert(len(sbox.data.all()) == len(sbox2.all()))
    for i in range(len(sbox.data.all())):
      assert(sbox.data.all()[i] == sbox2.all()[i])

    EPS = 1e-7
    for i in range(len(sbox.data)):
      v1 = sbox.data[i]
      v2 = sbox2[i]
      assert(abs(v1 - v2) < EPS)

  def check_range(self, sbox, i0, i1):
    for i in range(i0, i1):
      self.check(sbox[i-i0], i)

  def check_many(self, sbox, indices):
    for i in range(len(sbox)):
      self.check(sbox[i], indices[i])


if __name__ == '__main__':
  test = Test()
  test.run()
