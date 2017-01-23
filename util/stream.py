from __future__ import absolute_import, division

import logging
logger = logging.getLogger(__name__)


class ZMQStream():
  '''
  A class to listen to a zeromq stream

  '''
  def __init__(self, host, port=9999):
    """
    create stream listener object

    """
    self.host = host
    self.port = port

    # Connect to the stream
    self.connect()

  def connect(self):
    """
    Enable stream, connect to zmq host

    :return: receiver object

    """
    import zmq

    # Set the URL
    url = "tcp://{0}:{1}".format(self.host, self.port)
    logger.info("Connecting to ZMQ stream at %s" % url)

    # Get the zmq context
    context = zmq.Context()

    # Create the reciever
    self.receiver = context.socket(zmq.PULL)
    self.receiver.connect(url)
    return self.receiver

  def receive(self):
    """
    Receive and return zmq frames

    """
    logger.info("Waiting for frame")
    frames = self.receiver.recv_multipart(copy = False)
    return frames

  def close(self):
    """
    Close and disable stream

    """
    logger.info("Closing stream")
    return self.receiver.close()


class Result(object):
  '''
  A class to represent a result

  '''
  def __init__(self):
    pass

  def is_header(self):
    '''
    Return not a header

    '''
    return False

  def is_image(self):
    '''
    Return not an image

    '''
    return False

  def is_endofseries(self):
    '''
    Return not an endofseries

    '''
    return False


class Header(Result):
  '''
  A class representing a header

  '''

  def __init__(self, frames, directory, image_template):
    '''
    Create a header object

    '''
    import json

    super(Header, self).__init__()

    # Load the header
    head = json.loads(frames[0].bytes)

    # Get info
    if head['header_detail'] == 'basic':
      conf = json.loads(frames[1].bytes)
    elif head['header_detail'] == 'all':
      conf = json.loads(frames[1].bytes)
      # TODO use flatfield, pixel mask and count rate
      # flat_head = json.loads(frames[2].bytes)
      # flat_data = frames[3].bytes
      # mask_head = json.loads(frames[5].bytes)
      # mask_data = frames[6].bytes
      # rate_head = json.loads(frames[7].bytes)
      # rate_data = frames[8].bytes
    else:
      raise RuntimeError('Need header_detail to be "all" or "basic"')
    print conf

    # Save the header data
    self.header = {
      'htype' : 'eiger-stream',
      'directory' : directory,
      'image_template' : image_template,
      'configuration' : conf
    }

  def is_header(self):
    '''
    This is a header object

    '''
    return True

  def as_imageset(self, filename):
    '''
    Return as an imageset

    '''
    from dxtbx.format.FormatEigerStream import FormatEigerStream
    from dxtbx.imageset import MultiFileReader, ImageSweep

    # Create the reader
    reader = MultiFileReader(FormatEigerStream, [filename])

    # Create the sweep
    return ImageSweep(reader)


class Image(Result):
  '''
  A class to represent an image

  '''

  def __init__(self, frames, header):
    '''
    Create the image object

    '''
    import json

    super(Image, self).__init__()

    # Load stuff
    head = json.loads(frames[0].bytes)
    info = json.loads(frames[1].bytes)
    data = frames[2].bytes
    time = json.loads(frames[3].bytes)

    # The image number and data
    self.count = head['frame']
    self.data = data
    self.info = info

    # The dimensions
    shape = info['shape']
    dtype = info['type']
    denco = info['encoding']
    dsize = info['size']

    # The timing info
    self.start_time = time['start_time']
    self.stop_time = time['stop_time']
    self.real_time = time['real_time']

    # Check the sizes
    assert shape[1] == header.header['configuration']['x_pixels_in_detector']
    assert shape[0] == header.header['configuration']['y_pixels_in_detector']

    # # Check the image data hash
    # assert hashlib.md5(data).hexdigest() == head['hash']
    # if header.compression == 'lz4':
    #   assert denco == 'lz4<'
    # elif header.compression == 'bslz4':
    #   assert denco == 'bs32-lz4<'
    # else:
    #   raise RuntimeError('Unknown compression')

  def is_image(self):
    '''
    Return that the object is an image

    '''
    return True


class EndOfSeries(Result):
  '''
  Class to represent the end of a series

  '''

  def is_endofseries(self):
    '''
    Get the end of series

    '''
    return True


class Decoder(object):
  """
  Decodes zmq frames from EIGER ZMQ stream into dxtbx format object

  """

  def __init__(self, directory, image_template):
    '''
    Initialize the processor

    '''
    self.header = None
    self.directory = directory
    self.image_template = image_template

  def decode(self, frames):
    """
    Decode and process EIGER ZMQ stream frames

    """
    import json
    header = json.loads(frames[0].bytes)
    if header["htype"].startswith("dheader-"):
      return self.decode_header(frames)
    elif header["htype"].startswith("dimage-"):
      return self.decode_image(frames)
    elif header["htype"].startswith("dseries_end"):
      return self.decode_endofseries(frames)
    else:
      raise RuntimeError("No EIGER ZMQ message received")

  def decode_header(self, frames):
    '''
    Decode a header message

    '''
    self.header = Header(frames, self.directory, self.image_template)
    return self.header

  def decode_image(self, frames):
    '''
    Decode an image object

    '''
    if self.header is None:
      raise RuntimeError('Need header information before reading images')
    return Image(frames, self.header)

  def decode_endofseries(self, frames):
    '''
    Decode an endofseries object

    '''
    return EndOfSeries()
