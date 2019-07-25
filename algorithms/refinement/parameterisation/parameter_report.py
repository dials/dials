from __future__ import absolute_import, division, print_function


class ParameterReporter(object):
    """
    Keeps a record of all the ModelParameterisations and
    ScanVaryingModelParameterisations present and provides access to their
    Parameters and ScanVaryingParameterSets for reporting purposes.

    It is assumed that the provided model parameterisations will be one of five
    types:

    * Detector parameterisation
    * Beam parameterisation
    * Crystal orientation parameterisation
    * Crystal unit cell parameterisation
    * Goniometer setting parameterisation
    """

    def __init__(
        self,
        detector_parameterisations=None,
        beam_parameterisations=None,
        xl_orientation_parameterisations=None,
        xl_unit_cell_parameterisations=None,
        goniometer_parameterisations=None,
    ):

        if detector_parameterisations is None:
            detector_parameterisations = []
        if beam_parameterisations is None:
            beam_parameterisations = []
        if xl_orientation_parameterisations is None:
            xl_orientation_parameterisations = []
        if xl_unit_cell_parameterisations is None:
            xl_unit_cell_parameterisations = []
        if goniometer_parameterisations is None:
            goniometer_parameterisations = []

        # Keep references to all parameterised models
        self._detector_parameterisations = detector_parameterisations
        self._beam_parameterisations = beam_parameterisations
        self._xl_orientation_parameterisations = xl_orientation_parameterisations
        self._xl_unit_cell_parameterisations = xl_unit_cell_parameterisations
        self._goniometer_parameterisations = goniometer_parameterisations

        self._length = self._len()

    def _len(self):

        length = 0
        for model in self._detector_parameterisations:
            length += model.num_free()
        for model in self._beam_parameterisations:
            length += model.num_free()
        for model in self._xl_orientation_parameterisations:
            length += model.num_free()
        for model in self._xl_unit_cell_parameterisations:
            length += model.num_free()
        for model in self._goniometer_parameterisations:
            length += model.num_free()

        return length

    def __len__(self):
        return self._length

    def _indent(self, string):
        return "\n".join("    " + e for e in str(string).split("\n"))

    def __str__(self):

        s = "Parameter Report:\n"
        if self._detector_parameterisations:
            s += "Detector parameters:\n"
            det_plists = [x.get_params() for x in self._detector_parameterisations]
            params = [x for l in det_plists for x in l]
            for p in params:
                tmp = self._indent(p)
                s += tmp + "\n"

        if self._beam_parameterisations:
            s += "Beam parameters:\n"
            beam_plists = [x.get_params() for x in self._beam_parameterisations]
            params = [x for l in beam_plists for x in l]
            for p in params:
                tmp = self._indent(p)
                s += tmp + "\n"

        if self._xl_orientation_parameterisations:
            s += "Crystal orientation parameters:\n"
            xlo_plists = [
                x.get_params() for x in self._xl_orientation_parameterisations
            ]
            params = [x for l in xlo_plists for x in l]
            for p in params:
                tmp = self._indent(p)
                s += tmp + "\n"

        if self._xl_unit_cell_parameterisations:
            s += "Crystal unit cell parameters:\n"
            xluc_plists = [x.get_params() for x in self._xl_unit_cell_parameterisations]
            params = [x for l in xluc_plists for x in l]
            for p in params:
                tmp = self._indent(p)
                s += tmp + "\n"

        if self._goniometer_parameterisations:
            s += "Goniometer parameters:\n"
            gon_plists = [x.get_params() for x in self._goniometer_parameterisations]
            params = [x for l in gon_plists for x in l]
            for p in params:
                tmp = self._indent(p)
                s += tmp + "\n"

        return s

    def varying_params_vs_image_number(self, image_range):
        """Returns a string which is a table of scan-varying parameter values vs
        image number, if scan-varying parameters are present. Otherwise returns
        None"""

        image_numbers = list(range(image_range[0], image_range[1] + 1))
        columns = [TableColumn("Image", image_numbers)]

        for parameterisation in (
            self._detector_parameterisations
            + self._beam_parameterisations
            + self._xl_orientation_parameterisations
            + self._xl_unit_cell_parameterisations
            + self._goniometer_parameterisations
        ):
            for p in parameterisation.get_params():
                try:
                    vals = [
                        parameterisation.get_smoothed_parameter_value(i, p)
                        for i in image_numbers
                    ]
                    columns.append(TableColumn(p.name_stem, vals))
                except AttributeError:
                    continue

        if len(columns) > 1:
            header = "\t".join([e.title for e in columns])
            text = header + "\n"
            for i in range(len(columns[0])):
                vals = "\t".join(["%.6f" % e.values[i] for e in columns])
                text += vals + "\n"
            return text

        else:
            return None

    def get_params(self, only_free=True):
        """return a concatenated list of parameters from each of the components
        in the global model"""

        global_p_list = []
        for parameterisation in (
            self._detector_parameterisations
            + self._beam_parameterisations
            + self._xl_orientation_parameterisations
            + self._xl_unit_cell_parameterisations
            + self._goniometer_parameterisations
        ):
            global_p_list.extend(parameterisation.get_params(only_free))

        return global_p_list


class TableColumn(object):
    """Bucket to store data to be used for constructing tables to print."""

    def __init__(self, title, values):

        self._title = title
        self._values = values

    def __len__(self):
        return len(self._values)

    @property
    def title(self):
        return self._title

    @property
    def values(self):
        return self._values
