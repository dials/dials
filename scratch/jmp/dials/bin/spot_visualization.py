
class SpotVisualization (object):

    def __init__(self):
        self.vmin = 0
        self.vmax = 1000
        self.background_colour = (1.0, 1.0, 1.0)
        self.roi_colour = (0.19, 0.55, 0.91)
        self.window_size = (500, 500)

    def numpy_array_as_vtk_volume(self, image):
        """Create a VTK volume from the 3D numpy array."""
        import vtk
        import numpy

        # Normalize the array and convert to uint8
        max_image = numpy.max(image)
        min_image = numpy.min(image)
        b = 255.0 / (max_image - min_image)
        a = -min_image * b
        image = (a + b * image).astype(numpy.uint8)

        # Get the shape of the array
        size = image.shape

        # Convery the array to string to import into vtk
        image_string = image.tostring()

        # Create the data importer and import the numpy array
        dataImporter = vtk.vtkImageImport()
        dataImporter.CopyImportVoidPointer(image_string,len(image_string))
        dataImporter.SetDataScalarTypeToUnsignedChar()
        dataImporter.SetNumberOfScalarComponents(1)
        dataImporter.SetDataExtent(0, size[2]-1, 0, size[1]-1, 0, size[0]-1)
        dataImporter.SetWholeExtent(0, size[2]-1, 0, size[1]-1, 0, size[0]-1)

        # Set the volume alpha function
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(0, 0.0)
        alphaChannelFunc.AddPoint(127, 0.2)
        alphaChannelFunc.AddPoint(255, 1.0)

        # Set the volume colour function
        colorFunc = vtk.vtkColorTransferFunction()
        colorFunc.AddRGBPoint(0, 1.0, 1.0, 1.0)
        colorFunc.AddRGBPoint(127, 1.0, 1.0, 1.0)
        colorFunc.AddRGBPoint(255, 1.0, 0.0, 0.0)

        # Set the volume properties
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(colorFunc)
        volumeProperty.SetScalarOpacity(alphaChannelFunc)
        compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()

        # Create the volume mapper
        volumeMapper = vtk.vtkVolumeRayCastMapper()
        volumeMapper.SetVolumeRayCastFunction(compositeFunction)
        volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

        # Create the vtk volume
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)

        # Return the volume
        return volume

    def create_cube_actor_from_roi(self, roi):
        """Create a cube actor from the reflection roi"""
        import vtk

        # Get the bounds from the
        bounds = (roi[0], roi[1]-1, roi[2], roi[3]-1, roi[4], roi[5]-1)

        # Create the cube
        cube = vtk.vtkCubeSource()
        cube.SetBounds(bounds)

        # Create the mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cube.GetOutputPort())

        # Create the actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetRepresentationToWireframe()
        actor.GetProperty().SetColor(self.roi_colour)

        # Return the actor
        return actor

    def render_reflections(self, volume, roi):
        """Render the reflections"""
        import vtk

        # Create the renderer
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(self.background_colour)

        # Add the volume
        if volume != None:
            renderer.AddVolume(volume)

        # Add the regions of interest
        if roi != None:
            if isinstance(roi, list):
                for actor in roi:
                    renderer.AddActor(actor)
            else:
                renderer.AddActor(roi)

        # Create the render window
        render_win = vtk.vtkRenderWindow()
        render_win.AddRenderer(renderer)
        render_win.SetSize(self.window_size)
        render_win.SetWindowName("Reflections (q to quit)")

        # Create the render window interactor
        render_interactor = vtk.vtkRenderWindowInteractor()
        render_interactor.SetRenderWindow(render_win)
        render_interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        render_interactor.Initialize()
        render_win.Render()
        render_interactor.Start()

    def visualize_reflections(self, image, roi):
        """Visualise the reflections"""
        import numpy

        if image != None:
            # Cap image between vmin and vmax
            min_ind = numpy.where(image < self.vmin)
            max_ind = numpy.where(image > self.vmax)
            image[min_ind] = self.vmin
            image[max_ind] = self.vmax

            # Create the VTK volume image
            volume = self.numpy_array_as_vtk_volume(image)
        else:
            volume = None

        # Create the roi actors
        if roi != None:
            if isinstance(roi, list):
                roi_actor = []
                for r in roi:
                    roi_actor.append(self.create_cube_actor_from_roi(r))
            else:
                roi_actor = self.create_cube_actor_from_roi(roi)
        else:
            roi_actor = None

        # Render the reflections
        self.render_reflections(volume, roi_actor)


def create_spot(dim, A = 1.0, sig = (1.0, 1.0, 1.0)):
    import numpy
    from math import exp
    sig_x = sig[2]
    sig_y = sig[1]
    sig_z = sig[0]
    x0 = dim[2] / 2
    y0 = dim[1] / 2
    z0 = dim[0] / 2
    g = numpy.zeros(shape=dim, dtype=numpy.float64)
    for z in range(dim[0]):
        for y in range(dim[1]):
            for x in range(dim[2]):
                g[z, y, x] = A * exp(-(((x - x0)**2) / (2.0 * sig_x**2) +
                                       ((y - y0)**2) / (2.0 * sig_y**2) +
                                       ((z - z0)**2) / (2.0 * sig_z**2)))

    return g

def test_visualization():
    """Test the visualisation"""
    import numpy

    # Create the reflection image
    image = create_spot((10, 10, 10), 255, (5.0, 5.0, 5.0))
    image = image.astype(numpy.uint8)
    roi = (0, 10, 0, 10, 0, 10)

    v = SpotVisualization()
    v.visualize_reflections(image, roi)

if __name__ == '__main__':
    test_visualization()
