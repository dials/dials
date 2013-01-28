def spot(dim, A = 1.0, sig = (1.0, 1.0, 1.0)):
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

def numpy_array_as_vtk_volume(image):
    """Create a VTK volume from the 3D numpy array."""
    import vtk
    import numpy

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
    colorFunc.AddRGBPoint(127, 1.0, 0.0, 0.0)
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

def create_cube_actor_from_roi(roi):
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
    actor.GetProperty().SetColor((0.19, 0.55, 0.91))
    
    # Return the actor
    return actor

def render_reflections(volume, roi):
    """Render the reflections"""
    import vtk

    # Create the renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.AddVolume(volume)
    
    # Add the regions of interest
    if roi:
        if isinstance(roi, list):
            for actor in roi:
                renderer.AddActor(actor)
        else:
            renderer.AddActor(roi)
    
    # Create the render window
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    renderWin.SetSize(500, 500)
    renderWin.SetWindowName("Reflections (q to quit)")

    # Create the render window interactor
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)
    renderInteractor.Initialize()
    renderWin.Render()
    renderInteractor.Start()

def visualize_reflections():
    """Visualise the reflections"""
    import numpy

    # Create the reflection image
    image = spot((10, 10, 10), 255, (5.0, 5.0, 5.0))
    image = image.astype(numpy.uint8)
    
    roi = (0, 10, 0, 10, 0, 10)
    
    # Create the VTK volume image
    volume = numpy_array_as_vtk_volume(image)

    # Create the roi actor
    roi_actor = create_cube_actor_from_roi(roi)

    # Render the reflections
    render_reflections(volume, roi_actor)
    
if __name__ == '__main__':
    visualize_reflections()


