

class ImageFrame(wx.Frame):
  ''' Low-level image frame class.

  Display a single image Provide methods to manipulate the image frame and
  toggle various bits of information on and off.

  '''

  @property
  def image(self):
    return self._image

  @image.setter
  def image(self, value):
    return self._image = value


# Example
frame = ImageFrame(...)
frame.image = flex.double(flex.grid(width, height))
frame.xrange = (bbox[0], bbox[1])
frame.yrange = (bbox[2], bbox[3])
frame.points = [(xcentre, ycentre)]
frame.contour = [(xi, yi), ...]
frame.show_pixel_values = True
frame.show_points = True
frame.show_contour = True


class ReflectionFrame(wx.Frame):
  ''' Higher-level reflection frame class.

  Display information for a single reflection. Display each slice of the
  reflection in it's own ImageFrame instance and toggle higher-level items on
  and off, such as reflection centroid. Give frame a slider to allow moving
  through slices of reflection shoebox.
  '''
  pass

# Example
frame = ReflectionFrame(...)
frame.data = reflection_table[i]
frame.show_pixel_values = True # Show the pixel values in pixels
frame.show_centroid = True     # Display the centroid on each frame
frame.overlay_mask = True      # Show mask contours on data
frame.display_mode = 'data'    # ('data' | 'mask' | 'background')
frame.slice_mode = 'xy'        # ('xy' | 'xz' | 'yz')


def view_reflection(reflection, **kwargs):
  ''' A simple wrapper function that can be called from my code that takes a
  reflection (as a row of data) and displays it using the viewer. The kwargs
  can be used to configure the properties above. This would be really useful for
  debugging. '''
  pass

# Example
view_reflection(reflection_table[i], show_centroid=False, overlay_mask=True)
