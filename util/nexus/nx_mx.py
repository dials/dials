from __future__ import absolute_import, division, print_function

import math
from copy import deepcopy

from scitbx import matrix

# Extensions to NXMX
#
#  "detector/underload" - trusted_range[0]
#  "detector/timestamp" - epochs
#  "entry/index" - experiment index and shared model indices
#  "entry/template" - save location to external files (e.g. cbf)

schema_url = (
    "https://github.com/nexusformat/definitions/blob/master/applications/NXmx.nxdl.xml"
)


def convert_to_nexus_beam_direction(experiments):
    EPS = 1e-7

    zaxis = matrix.col((0, 0, -1))

    class Dummy:
        pass

    experiments2 = []

    # Save the sharedness as the id
    for exp in experiments:
        exp.beam_id = id(exp.beam)
        exp.detector_id = id(exp.detector)
        exp.goniometer_id = id(exp.goniometer)
        exp.scan_id = id(exp.scan)
        exp.crystal_id = id(exp.crystal)
        exp2 = Dummy()
        exp2.beam_id = id(exp.beam)
        exp2.detector_id = id(exp.detector)
        exp2.goniometer_id = id(exp.goniometer)
        exp2.scan_id = id(exp.scan)
        exp2.crystal_id = id(exp.crystal)
        experiments2.append(exp2)

    # Copy all models
    for exp, exp2 in zip(experiments, experiments2):
        exp2.beam = deepcopy(exp.beam)
        exp2.detector = deepcopy(exp.detector)
        exp2.goniometer = deepcopy(exp.goniometer)
        exp2.scan = deepcopy(exp.scan)
        exp2.crystal = deepcopy(exp.crystal)
        exp2.imageset = exp.imageset
    experiments = experiments2

    # Rotate the beams
    rotations = []
    for exp in experiments:
        d = matrix.col(exp.beam.get_direction()).normalize()
        angle = d.angle(zaxis, deg=False)
        if abs(angle - math.pi) < EPS:
            axis = (1, 0, 0)
        elif abs(angle) < EPS:
            rotations.append(((0, 0, 0), 0))
            continue
        else:
            axis = d.cross(zaxis).normalize()
        exp.beam.rotate_around_origin(axis, angle, deg=False)
        exp.detector.rotate_around_origin(axis, angle, deg=False)
        if exp.goniometer:
            exp.goniometer.rotate_around_origin(axis, angle, deg=False)
        exp.crystal.rotate_around_origin(axis, angle, deg=False)
        d = matrix.col(exp.beam.get_direction())
        assert abs(d.angle(zaxis)) < EPS
        rotations.append((axis, angle))

    # Return converted experiments
    return experiments, rotations


def convert_from_nexus_beam_direction(experiments, rotations):
    EPS = 1e-7

    zaxis = matrix.col((0, 0, -1))
    for exp, (axis, angle) in zip(experiments, rotations):
        d = matrix.col(exp.beam.get_direction()).normalize()
        assert abs(d.angle(zaxis)) < EPS
        exp.beam.rotate_around_origin(axis, angle, deg=False)
        exp.detector.rotate_around_origin(axis, angle, deg=False)
        if exp.goniometer:
            exp.goniometer.rotate_around_origin(axis, angle, deg=False)
        exp.crystal.rotate_around_origin(axis, angle, deg=False)
    return experiments


def polarization_normal_to_stokes(n, p):
    EPS = 1e-7
    ax = matrix.col((1, 0, 0))
    ay = matrix.col((0, 1, 0))
    az = matrix.col((0, 0, 1))
    n = matrix.col(n).normalize()
    assert abs(n.dot(az)) < EPS
    I = 1.0
    X = 0.0
    W = math.atan2(n.dot(ay), n.dot(ax)) - math.pi / 2
    S0 = I
    S1 = I * p * math.cos(2 * W) * math.cos(2 * X)
    S2 = I * p * math.sin(2 * W) * math.cos(2 * X)
    S3 = I * p * math.sin(2 * X)
    return S0, S1, S2, S3


def polarization_stokes_to_normal(S0, S1, S2, S3):
    EPS = 1e-7
    p = math.sqrt(S1 * S1 + S2 * S2 + S3 * S3) / S0
    W = math.atan2(S2, S1) * 0.5
    X = math.atan2(S3, math.sqrt(S1 * S1 + S2 * S2)) * 0.5
    assert abs(X) < EPS
    n = matrix.col((-math.sin(W), math.cos(W), 0.0)).normalize()
    return n, p


def get_nx_class(handle, klass, path):
    if path in handle:
        group = handle[path]
        assert group.attrs["NX_class"] == klass
    else:
        group = handle.create_group(path)
        group.attrs["NX_class"] = klass
    return group


def get_nx_sample(handle, path):
    return get_nx_class(handle, "NXsample", path)


def get_nx_beam(handle, path):
    return get_nx_class(handle, "NXbeam", path)


def get_nx_detector(handle, path):
    return get_nx_class(handle, "NXdetector", path)


def get_nx_detector_module(handle, path):
    return get_nx_class(handle, "NXdetector_module", path)


def get_nx_instrument(handle, path):
    return get_nx_class(handle, "NXinstrument", path)


def get_nx_transformations(handle, path):
    return get_nx_class(handle, "NXtransformations", path)


def get_nx_process(handle, path):
    return get_nx_class(handle, "NXprocess", path)


def get_nx_note(handle, path):
    return get_nx_class(handle, "NXnote", path)


def get_nx_data(handle, path, data):
    handle[path] = data
    return handle[path]


def get_nx_dials(handle, path):
    return get_nx_class(handle, "NXdials", path)


def dump_beam(entry, beam):
    """ Export the beam model. """
    EPS = 1e-7

    # Get the nx_beam
    nx_sample = get_nx_sample(entry, "sample")
    nx_beam = get_nx_beam(nx_sample, "beam")

    # Generate the stokes polarization parameters
    d = matrix.col(beam.get_direction()).normalize()
    n = matrix.col(beam.get_polarization_normal()).normalize()
    p = beam.get_polarization_fraction()
    assert abs(n.dot(d)) < EPS
    assert abs(d.dot(matrix.col((0, 0, -1))) - 1) < EPS

    # Get the polarization in stokes parameters
    S0, S1, S2, S3 = polarization_normal_to_stokes(n, p)

    # Set the beam parameters
    nx_beam["incident_wavelength"] = beam.get_wavelength()
    nx_beam["incident_polarization_stokes"] = (S0, S1, S2, S3)


def dump_detector(entry, detector, beam, imageset, scan):
    EPS = 1e-7

    # Get the detector
    nx_instrument = get_nx_instrument(entry, "instrument")
    nx_detector = get_nx_detector(nx_instrument, "detector")

    # Get some bulk properties
    thickness = detector[0].get_thickness()
    material = detector[0].get_material()
    dtype = detector[0].get_type()
    trusted_range = detector[0].get_trusted_range()

    # Check all panels obey these bulk properties
    assert [abs((p.get_thickness() - thickness) < EPS) for p in detector].count(
        False
    ) == 0
    assert [p.get_material() == material for p in detector].count(False) == 0
    assert [p.get_type() == dtype for p in detector].count(False) == 0
    assert [p.get_trusted_range() == trusted_range for p in detector].count(False) == 0

    # Take the distance and beam centre from the first panel
    # distance = detector[0].get_directed_distance()
    # beam_centre_x, beam_centre_y = detector[0].get_beam_centre(beam.get_s0())

    # Set the detector properties
    nx_detector["sensor_thickness"] = thickness
    nx_detector["sensor_material"] = material
    nx_detector["type"] = dtype
    # nx_detector['distance'] = distance
    # nx_detector['beam_centre_x'] = beam_centre_x
    # nx_detector['beam_centre_y'] = beam_centre_y
    nx_detector["saturation_value"] = int(trusted_range[1])
    nx_detector["underload"] = int(trusted_range[0])  # FIXME non-standard
    nx_detector["description"] = dtype

    # Make up some fake stuff
    if scan is not None:
        nx_detector[
            "timestamp"
        ] = scan.get_epochs().as_numpy_array()  # FIXME non-standard
        nx_detector["dead_time"] = [0.0] * len(scan)
        nx_detector["count_time"] = [0.0] * len(scan)
        nx_detector["frame_time"] = scan.get_exposure_times().as_numpy_array()
        nx_detector["detector_readout_time"] = [0.0] * len(scan)
    else:
        nx_detector["timestamp"] = 0.0
        nx_detector["dead_time"] = 0.0
        nx_detector["count_time"] = 0.0
        nx_detector["frame_time"] = 0.0
        nx_detector["detector_readout_time"] = 0.0
    nx_detector["bit_depth_readout"] = 32

    # Create the detector depends on
    nx_detector["depends_on"] = "."

    # Creat some data for example file
    # data = [im.as_numpy_array() for im in imageset]

    # Create the nx data
    # nx_data = get_nx_data(nx_detector, "data", [[[0]],[[0]],[[0]]])

    # Loop through all the panels
    for i, panel in enumerate(detector):

        # Get some panel attributes
        pixel_size = panel.get_pixel_size()
        image_size = panel.get_image_size()
        origin = matrix.col(panel.get_origin())

        # Get the detector module object
        nx_module = get_nx_detector_module(nx_detector, "module%d" % i)

        # Set the data size
        nx_module["data_size"] = image_size
        nx_module["data_origin"] = (-1, -1)  # FIXME INVALID

        # Set the module offset
        nx_module["module_offset"] = origin.length()
        nx_module["module_offset"].attrs["depends_on"] = "."
        nx_module["module_offset"].attrs["transformation_type"] = "translation"
        nx_module["module_offset"].attrs["offset"] = (0.0, 0.0, 0.0)
        nx_module["module_offset"].attrs["vector"] = origin.normalize()

        # The path for items below
        module_offset_path = str(nx_module["module_offset"].name)

        # Write the fast pixel direction
        nx_module["fast_pixel_direction"] = pixel_size[0]
        nx_module["fast_pixel_direction"].attrs["depends_on"] = module_offset_path
        nx_module["fast_pixel_direction"].attrs["transformation_type"] = "translation"
        nx_module["fast_pixel_direction"].attrs["offset"] = (0.0, 0.0, 0.0)
        nx_module["fast_pixel_direction"].attrs["vector"] = panel.get_fast_axis()

        # Write the slow pixel direction
        nx_module["slow_pixel_direction"] = pixel_size[1]
        nx_module["slow_pixel_direction"].attrs["depends_on"] = module_offset_path
        nx_module["slow_pixel_direction"].attrs["transformation_type"] = "translation"
        nx_module["slow_pixel_direction"].attrs["offset"] = (0.0, 0.0, 0.0)
        nx_module["slow_pixel_direction"].attrs["vector"] = panel.get_slow_axis()


def dump_goniometer(entry, goniometer, scan):
    """ Export the goniometer model. """
    if scan is None or goniometer is None:
        return

    # The angles for each image
    phi0, dphi = scan.get_oscillation(deg=True)
    phi = [phi0 + dphi * i for i in range(len(scan))]

    nx_sample = get_nx_sample(entry, "sample")
    nx_transformations = get_nx_transformations(nx_sample, "transformations")

    # Write out the rotation axis and oscillation
    r = matrix.sqr(goniometer.get_setting_rotation())
    q = r.r3_rotation_matrix_as_unit_quaternion()
    angle, axis = q.unit_quaternion_as_axis_and_angle()
    nx_transformations["setting_rotation"] = angle
    nx_transformations["setting_rotation"].attrs["depends_on"] = "."
    nx_transformations["setting_rotation"].attrs["transformation_type"] = "rotation"
    nx_transformations["setting_rotation"].attrs["offset_units"] = "mm"
    nx_transformations["setting_rotation"].attrs["offset"] = (0.0, 0.0, 0.0)
    nx_transformations["setting_rotation"].attrs["vector"] = axis

    r = matrix.sqr(goniometer.get_fixed_rotation())
    q = r.r3_rotation_matrix_as_unit_quaternion()
    angle, axis = q.unit_quaternion_as_axis_and_angle()
    nx_transformations["fixed_rotation"] = angle
    nx_transformations["fixed_rotation"].attrs["depends_on"] = str(
        nx_transformations["setting_rotation"].name
    )
    nx_transformations["fixed_rotation"].attrs["transformation_type"] = "rotation"
    nx_transformations["fixed_rotation"].attrs["offset_units"] = "mm"
    nx_transformations["fixed_rotation"].attrs["offset"] = (0.0, 0.0, 0.0)
    nx_transformations["fixed_rotation"].attrs["vector"] = axis

    nx_transformations["phi"] = phi
    nx_transformations["phi"].attrs["depends_on"] = str(
        nx_transformations["fixed_rotation"].name
    )
    nx_transformations["phi"].attrs["transformation_type"] = "rotation"
    nx_transformations["phi"].attrs["offset_units"] = "mm"
    nx_transformations["phi"].attrs["offset"] = (0.0, 0.0, 0.0)
    nx_transformations["phi"].attrs["vector"] = goniometer.get_rotation_axis_datum()


def dump_crystal(entry, crystal, scan):
    """ Export the crystal model. """
    from scitbx.array_family import flex

    # Get the sample
    nx_sample = get_nx_sample(entry, "sample")

    # Set the space group
    nx_sample["unit_cell_group"] = crystal.get_space_group().type().hall_symbol()

    # Get the unit cell and orientation matrix in the case of scan varying and
    # scan static models
    if crystal.num_scan_points:
        num = crystal.num_scan_points
        unit_cell = flex.double(flex.grid(num, 6))
        orientation_matrix = flex.double(flex.grid(num, 9))
        for i in range(num):
            __cell = crystal.get_unit_cell_at_scan_point(i).parameters()
            for j in range(6):
                unit_cell[i, j] = __cell[j]
            __matrix = crystal.get_U_at_scan_point(i)
            for j in range(9):
                orientation_matrix[i, j] = __matrix[j]
        orientation_matrix = [
            [
                tuple(orientation_matrix[i : i + 1, 0:3]),
                tuple(orientation_matrix[i : i + 1, 3:6]),
                tuple(orientation_matrix[i : i + 1, 6:9]),
            ]
            for i in range(num)
        ]
        unit_cell = [tuple(unit_cell[i : i + 1, :]) for i in range(num)]
        average_unit_cell = crystal.get_unit_cell().parameters()
        average_orientation_matrix = crystal.get_U()
        average_orientation_matrix = [
            average_orientation_matrix[0:3],
            average_orientation_matrix[3:6],
            average_orientation_matrix[6:9],
        ]

    else:
        unit_cell = [crystal.get_unit_cell().parameters()]
        orientation_matrix = [crystal.get_U()]
        orientation_matrix = [
            [
                tuple(orientation_matrix[0][0:3]),
                tuple(orientation_matrix[0][3:6]),
                tuple(orientation_matrix[0][6:9]),
            ]
        ]
        average_unit_cell = unit_cell[0]
        average_orientation_matrix = orientation_matrix[0]

    # Save the unit cell data
    nx_sample["name"] = "FROM_DIALS"
    nx_sample["unit_cell"] = unit_cell
    nx_sample["unit_cell"].attrs["angles_units"] = "deg"
    nx_sample["unit_cell"].attrs["length_units"] = "angstrom"

    # Save the orientation matrix
    nx_sample["orientation_matrix"] = orientation_matrix

    # Set an average unit cell etc for scan static stuff
    nx_sample["average_unit_cell"] = average_unit_cell
    nx_sample["average_unit_cell"].attrs["angles_units"] = "deg"
    nx_sample["average_unit_cell"].attrs["length_units"] = "angstrom"
    nx_sample["average_orientation_matrix"] = average_orientation_matrix

    # Set depends on
    if scan is not None:
        nx_sample["depends_on"] = str(nx_sample["transformations/phi"].name)
    else:
        nx_sample["depends_on"] = "."


def dump_details(entry):
    from time import strftime

    # Program info
    entry["program_name"] = "dials.export_nxmx"
    entry["program_name"].attrs["version"] = 1
    entry["program_name"].attrs["configuration"] = ""

    # Set some processing information (each program should add itself)
    nx_process = get_nx_process(entry, "process")
    nx_process["program"] = "dials"
    nx_process["version"] = 1
    nx_process["date"] = strftime("%Y-%m-%dT%H:%M:%S")

    nx_note = get_nx_note(nx_process, "spot_finding")
    nx_note["author"] = "dials.find_spots"
    nx_note["date"] = strftime("%Y-%m-%dT%H:%M:%S")
    nx_note["type"] = "text/plain"
    nx_note["description"] = "Spot finding parameters"
    nx_note["data"] = "dials.find_spots imported.expt"
    nx_note["sequence_index"] = 0

    nx_note = get_nx_note(nx_process, "indexing")
    nx_note["author"] = "dials.index"
    nx_note["date"] = strftime("%Y-%m-%dT%H:%M:%S")
    nx_note["type"] = "text/plain"
    nx_note["description"] = "Indexing parameters"
    nx_note["data"] = "dials.index imported.expt strong.refl"
    nx_note["sequence_index"] = 1

    nx_note = get_nx_note(nx_process, "refinement")
    nx_note["author"] = "dials.refine"
    nx_note["date"] = strftime("%Y-%m-%dT%H:%M:%S")
    nx_note["type"] = "text/plain"
    nx_note["description"] = "Refinement parameters"
    nx_note["data"] = "dials.refine indexed.expt indexed.refl"
    nx_note["sequence_index"] = 2

    nx_note = get_nx_note(nx_process, "integration")
    nx_note["author"] = "dials.integrate"
    nx_note["date"] = strftime("%Y-%m-%dT%H:%M:%S")
    nx_note["type"] = "text/plain"
    nx_note["description"] = "Integration parameters"
    nx_note["data"] = "dials.integrate refined.expt refined.refl"
    nx_note["sequence_index"] = 3


def load_beam(entry):
    from dxtbx.model import Beam

    EPS = 1e-7

    # Get the nx_beam
    nx_sample = get_nx_sample(entry, "sample")
    nx_beam = get_nx_beam(nx_sample, "beam")
    wavelength = nx_beam["incident_wavelength"][()]
    S0, S1, S2, S3 = tuple(nx_beam["incident_polarization_stokes"])
    n, p = polarization_stokes_to_normal(S0, S1, S2, S3)
    assert n.dot(matrix.col((0, 0, -1))) < EPS

    # Return the beam model
    return Beam((0, 0, -1), wavelength, 0, 0, n, p, 0, 1)


def load_detector(entry):
    from dxtbx.model import Detector

    # Get the detector module object
    nx_instrument = get_nx_instrument(entry, "instrument")
    nx_detector = get_nx_detector(nx_instrument, "detector")
    assert nx_detector["depends_on"][()] == "."
    material = nx_detector["sensor_material"][()]
    det_type = nx_detector["type"][()]
    thickness = nx_detector["sensor_thickness"][()]
    trusted_range = (nx_detector["underload"][()], nx_detector["saturation_value"][()])

    # The detector model
    detector = Detector()

    i = 0
    while True:
        try:
            module = get_nx_detector_module(nx_detector, "module%d" % i)
        except Exception:
            break
        # Set the data size
        image_size = module["data_size"]

        # Set the module offset
        offset_length = module["module_offset"][()]
        assert module["module_offset"].attrs["depends_on"] == "."
        assert module["module_offset"].attrs["transformation_type"] == "translation"
        assert tuple(module["module_offset"].attrs["offset"]) == (0, 0, 0)
        offset_vector = matrix.col(module["module_offset"].attrs["vector"])
        origin = offset_vector * offset_length

        # Write the fast pixel direction
        module_offset_path = str(module["module_offset"].name)
        pixel_size_x = module["fast_pixel_direction"][()]
        assert module["fast_pixel_direction"].attrs["depends_on"] == module_offset_path
        assert (
            module["fast_pixel_direction"].attrs["transformation_type"] == "translation"
        )
        assert tuple(module["fast_pixel_direction"].attrs["offset"]) == (0, 0, 0)
        fast_axis = tuple(module["fast_pixel_direction"].attrs["vector"])

        # Write the slow pixel direction
        pixel_size_y = module["slow_pixel_direction"][()]
        assert module["slow_pixel_direction"].attrs["depends_on"] == module_offset_path
        assert (
            module["slow_pixel_direction"].attrs["transformation_type"] == "translation"
        )
        assert tuple(module["slow_pixel_direction"].attrs["offset"]) == (0, 0, 0)
        slow_axis = tuple(module["slow_pixel_direction"].attrs["vector"])

        # Get the pixel size and axis vectors
        pixel_size = (pixel_size_x, pixel_size_y)

        # Create the panel
        panel = detector.add_panel()
        panel.set_frame(fast_axis, slow_axis, origin)
        panel.set_pixel_size(pixel_size)
        panel.set_image_size(image_size)
        panel.set_type(det_type)
        panel.set_thickness(thickness)
        panel.set_material(material)
        panel.set_trusted_range(trusted_range)
        i += 1

    # Return the detector and panel
    return detector


def load_goniometer(entry):
    from dxtbx.model import Goniometer

    # Write out the rotation axis and oscillation
    nx_sample = get_nx_sample(entry, "sample")
    try:
        transformations = get_nx_transformations(nx_sample, "transformations")
    except Exception:
        return None
    assert transformations["phi"].attrs["depends_on"] == str(
        transformations["fixed_rotation"].name
    )
    assert transformations["phi"].attrs["transformation_type"] == "rotation"
    assert transformations["phi"].attrs["offset_units"] == "mm"
    assert tuple(transformations["phi"].attrs["offset"]) == (0, 0, 0)
    rotation_axis = tuple(transformations["phi"].attrs["vector"])

    assert transformations["fixed_rotation"].attrs["depends_on"] == str(
        transformations["setting_rotation"].name
    )
    assert transformations["fixed_rotation"].attrs["transformation_type"] == "rotation"
    assert transformations["fixed_rotation"].attrs["offset_units"] == "mm"
    assert tuple(transformations["phi"].attrs["offset"]) == (0, 0, 0)
    axis = matrix.col(transformations["fixed_rotation"].attrs["vector"])
    angle = transformations["fixed_rotation"][()]
    fixed_rotation = axis.axis_and_angle_as_r3_rotation_matrix(angle)

    assert transformations["setting_rotation"].attrs["depends_on"] == "."
    assert (
        transformations["setting_rotation"].attrs["transformation_type"] == "rotation"
    )
    assert transformations["setting_rotation"].attrs["offset_units"] == "mm"
    assert tuple(transformations["phi"].attrs["offset"]) == (0, 0, 0)
    axis = matrix.col(transformations["setting_rotation"].attrs["vector"])
    angle = transformations["setting_rotation"][()]
    setting_rotation = axis.axis_and_angle_as_r3_rotation_matrix(angle)

    # Return the goniometer model
    return Goniometer(rotation_axis, fixed_rotation, setting_rotation)


def load_scan(entry):
    from dxtbx.model import Scan

    # Write out the rotation axis and oscillation
    nx_sample = get_nx_sample(entry, "sample")
    try:
        transformations = get_nx_transformations(nx_sample, "transformations")
    except Exception:
        return None
    phi = transformations["phi"]
    assert transformations["phi"].attrs["transformation_type"] == "rotation"
    assert transformations["phi"].attrs["offset_units"] == "mm"
    assert tuple(transformations["phi"].attrs["offset"]) == (0, 0, 0)
    image_range = (1, len(phi))
    oscillation = (phi[0], phi[1] - phi[0])
    nx_instrument = get_nx_instrument(entry, "instrument")
    nx_detector = get_nx_detector(nx_instrument, "detector")
    exposure_time = nx_detector["frame_time"]
    epochs = nx_detector["timestamp"]
    return Scan(image_range, oscillation, exposure_time, epochs, deg=True)


def load_crystal(entry):
    from dxtbx.model import Crystal
    from scitbx.array_family import flex
    from cctbx import uctbx
    import numpy

    # Get the sample
    nx_sample = get_nx_sample(entry, "sample")

    # Set the space group
    space_group_symbol = nx_sample["unit_cell_group"][()]

    # Get depends on
    if nx_sample["depends_on"][()] != ".":
        assert nx_sample["depends_on"][()] == str(nx_sample["transformations/phi"].name)

    # Read the average unit cell data
    average_unit_cell = flex.double(numpy.array(nx_sample["average_unit_cell"]))
    assert nx_sample["average_unit_cell"].attrs["angles_units"] == "deg"
    assert nx_sample["average_unit_cell"].attrs["length_units"] == "angstrom"
    assert len(average_unit_cell.all()) == 1
    assert len(average_unit_cell) == 6
    average_orientation_matrix = flex.double(
        numpy.array(nx_sample["average_orientation_matrix"])
    )
    assert len(average_orientation_matrix.all()) == 2
    assert average_orientation_matrix.all()[0] == 3
    assert average_orientation_matrix.all()[1] == 3

    # Get the real space vectors
    uc = uctbx.unit_cell(tuple(average_unit_cell))
    U = matrix.sqr(average_orientation_matrix)
    B = matrix.sqr(uc.fractionalization_matrix()).transpose()
    A = U * B
    A = A.inverse()
    real_space_a = A[0:3]
    real_space_b = A[3:6]
    real_space_c = A[6:9]

    # Read the unit cell data
    unit_cell = flex.double(numpy.array(nx_sample["unit_cell"]))
    assert nx_sample["unit_cell"].attrs["angles_units"] == "deg"
    assert nx_sample["unit_cell"].attrs["length_units"] == "angstrom"

    # Read the orientation matrix
    orientation_matrix = flex.double(numpy.array(nx_sample["orientation_matrix"]))
    assert len(unit_cell.all()) == 2
    assert len(orientation_matrix.all()) == 3
    assert unit_cell.all()[0] == orientation_matrix.all()[0]
    assert unit_cell.all()[1] == 6
    assert orientation_matrix.all()[1] == 3
    assert orientation_matrix.all()[2] == 3

    # Construct the crystal model
    crystal = Crystal(real_space_a, real_space_b, real_space_c, space_group_symbol)

    # Sort out scan points
    if unit_cell.all()[0] > 1:
        A_list = []
        for i in range(unit_cell.all()[0]):
            uc = uctbx.unit_cell(tuple(unit_cell[i : i + 1, :]))
            U = matrix.sqr(tuple(orientation_matrix[i : i + 1, :, :]))
            B = matrix.sqr(uc.fractionalization_matrix()).transpose()
            A_list.append(U * B)
        crystal.set_A_at_scan_points(A_list)
    else:
        assert unit_cell.all_eq(average_unit_cell)
        assert orientation_matrix.all_eq(average_orientation_matrix)

    # Return the crystal
    return crystal


def dump(entry, experiments):
    from dxtbx.imageset import ImageSweep

    print("Dumping NXmx")

    # Rotate the experiments such that beam direction is along (0, 0, -1)
    experiments, rotations = convert_to_nexus_beam_direction(experiments)

    # Add the feature
    if "features" in entry:
        features = entry["features"]
        assert features.dtype == "uint64"
        features.resize((len(features) + 1,))
        features[len(features) - 1] = 6
    else:
        import numpy as np

        features = entry.create_dataset(
            "features", (1,), maxshape=(None,), dtype=np.uint64
        )
        features[0] = 6

    exp_names = []

    # Get the experiment
    for index, experiment in enumerate(experiments):

        # Create the entry
        assert ("experiment_%d" % index) not in entry
        nxmx = entry.create_group("experiment_%d" % index)
        nxmx.attrs["NX_class"] = "NXsubentry"
        exp_names.append(str(nxmx.name))

        # Get the dials specific stuff
        nx_dials = get_nx_dials(nxmx, "dials")
        nx_dials["index"] = index
        nx_dials["index"].attrs["source"] = experiment.beam_id
        nx_dials["index"].attrs["detector"] = experiment.detector_id
        if experiment.goniometer is not None:
            nx_dials["index"].attrs["goniometer"] = experiment.goniometer_id
        if experiment.scan is not None:
            nx_dials["index"].attrs["scan"] = experiment.scan_id
        nx_dials["index"].attrs["sample"] = experiment.crystal_id

        # Write out the original orientation (dials specific)
        transformations = get_nx_transformations(nx_dials, "transformations")
        transformations["angle"] = -rotations[index][1]
        transformations["angle"].attrs["transformation_type"] = "rotation"
        transformations["angle"].attrs["vector"] = rotations[index][0]
        transformations["angle"].attrs["offset"] = (0, 0, 0)
        transformations["angle"].attrs["offset_units"] = "mm"
        transformations["angle"].attrs["depends_on"] = "."

        # Create the imageset template
        if experiment.imageset is None:
            nx_dials["template"] = ""
            if experiment.scan is not None:
                nx_dials["template"].attrs["range"] = experiment.scan.get_image_range()
        else:
            from os.path import abspath

            if isinstance(experiment.imageset, ImageSweep):
                template = abspath(experiment.imageset.get_template())
                nx_dials["template"] = template
                nx_dials["template"].attrs["range"] = experiment.scan.get_image_range()
            else:
                template = [
                    abspath(experiment.imageset.get_path(i))
                    for i in range(len(experiment.imageset))
                ]
                nx_dials["template"] = template

        # Create the definition
        definition = nxmx.create_dataset("definition", data="NXmx")
        definition.attrs["version"] = 1
        definition.attrs["URL"] = schema_url

        nxmx["title"] = "FROM_DIALS"

        # Dump the models
        dump_beam(nxmx, experiment.beam)
        dump_detector(
            nxmx,
            experiment.detector,
            experiment.beam,
            experiment.imageset,
            experiment.scan,
        )
        dump_goniometer(nxmx, experiment.goniometer, experiment.scan)
        dump_crystal(nxmx, experiment.crystal, experiment.scan)

    # Dump some details
    dump_details(entry)

    # Link the data
    # nxmx['data'] = nxmx['instrument/detector/data']

    return exp_names


def find_nx_mx_entries(nx_file, entry):
    """
    Find NXmx entries

    """
    hits = []

    def visitor(name, obj):
        if "NX_class" in obj.attrs:
            if obj.attrs["NX_class"] in ["NXentry", "NXsubentry"]:
                if "definition" in obj:
                    if obj["definition"][()] == "NXmx":
                        hits.append(obj)

    nx_file[entry].visititems(visitor)
    return hits


def load(entry, exp_index):
    from dxtbx.model.experiment_list import ExperimentList
    from dxtbx.model.experiment_list import Experiment

    print("Loading NXmx")

    # Check file contains the feature
    assert "features" in entry
    assert 6 in entry["features"][()]

    experiment_list = ExperimentList()

    # Find all the experiments
    entries = find_nx_mx_entries(entry, ".")
    if len(entries) > 1:
        entries = sorted(entries, key=lambda x: x["dials/index"][()])

    assert len(entries) == len(exp_index)
    for nxmx, name in zip(entries, exp_index):
        assert nxmx.name == name

    index = []
    rotations = []
    for name in exp_index:

        # Get the entry
        nxmx = entry.file[name]

        # Get the definition
        definition = nxmx["definition"]
        assert definition[()] == "NXmx"
        assert definition.attrs["version"] == 1

        # Get dials specific stuff
        nx_dials = get_nx_dials(nxmx, "dials")

        # Set index
        b = nx_dials["index"].attrs["source"]
        d = nx_dials["index"].attrs["detector"]
        if "goniometer" in nx_dials["index"].attrs:
            g = nx_dials["index"].attrs["goniometer"]
        else:
            g = None
        if "scan" in nx_dials["index"].attrs:
            s = nx_dials["index"].attrs["scan"]
        else:
            s = None
        c = nx_dials["index"].attrs["sample"]
        index.append((b, d, g, s, c))

        # Get the original orientation (dials specific)
        transformations = get_nx_transformations(nx_dials, "transformations")
        angle = transformations["angle"][()]
        assert transformations["angle"].attrs["transformation_type"] == "rotation"
        axis = transformations["angle"].attrs["vector"]
        assert tuple(transformations["angle"].attrs["offset"]) == (0, 0, 0)
        assert transformations["angle"].attrs["offset_units"] == "mm"
        assert transformations["angle"].attrs["depends_on"] == "."
        rotations.append((axis, angle))

        # Get the template and imageset
        try:
            template = list(nx_dials["template"])
            image_range = None
        except Exception:
            template = nx_dials["template"][()]
            if template == "":
                template = None
            if "range" in nx_dials["template"].attrs:
                image_range = nx_dials["template"].attrs["range"]
            else:
                image_range = None

        # Create the experiment
        experiment = Experiment()

        # Read the models
        experiment.beam = load_beam(nxmx)
        experiment.detector = load_detector(nxmx)
        experiment.goniometer = load_goniometer(nxmx)
        experiment.scan = load_scan(nxmx)
        experiment.crystal = load_crystal(nxmx)

        # Set the image range
        if image_range is not None and experiment.scan is not None:
            num = image_range[1] - image_range[0] + 1
            assert num == len(experiment.scan)
            experiment.scan.set_image_range(image_range)

        # Return the experiment list
        experiment_list.append(experiment)

    # Convert from nexus beam direction
    experiment_list = convert_from_nexus_beam_direction(experiment_list, rotations)

    from collections import defaultdict

    beam = defaultdict(list)
    detector = defaultdict(list)
    goniometer = defaultdict(list)
    scan = defaultdict(list)
    crystal = defaultdict(list)
    for i, ind in enumerate(index):
        beam[ind[0]].append(i)
        detector[ind[1]].append(i)
        goniometer[ind[2]].append(i)
        scan[ind[3]].append(i)
        crystal[ind[4]].append(i)

    # Set all the shared beams
    for value in beam.values():
        b1 = experiment_list[value[0]].beam
        assert all(experiment_list[v].beam == b1 for v in value[1:])
        for v in value[1:]:
            experiment_list[v].beam = b1
    # Set all the shared detectors
    for value in detector.values():
        d1 = experiment_list[value[0]].detector
        assert all(experiment_list[v].detector == d1 for v in value[1:])
        for v in value[1:]:
            experiment_list[v].detector = d1
    # Set all the shared goniometer
    for value in goniometer.values():
        g1 = experiment_list[value[0]].goniometer
        assert all(experiment_list[v].goniometer == g1 for v in value[1:])
        for v in value[1:]:
            experiment_list[v].goniometer = g1
    # Set all the shared scans
    for value in scan.values():
        s1 = experiment_list[value[0]].scan
        assert all(experiment_list[v].scan == s1 for v in value[1:])
        for v in value[1:]:
            experiment_list[v].scan = s1
    # Set all the shared crystals
    for value in crystal.values():
        c1 = experiment_list[value[0]].crystal
        assert all(experiment_list[v].crystal == c1 for v in value[1:])
        for v in value[1:]:
            experiment_list[v].crystal = c1

    return experiment_list
