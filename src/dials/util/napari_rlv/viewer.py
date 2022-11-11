from __future__ import annotations

from collections import namedtuple
from enum import Enum
from math import ceil, floor, pi

import napari
import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, Label
from napari.experimental import link_layers

import libtbx.phil
from dxtbx import flumpy
from libtbx import Auto
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import minimum_covering_sphere

from dials.util.napari_rlv import Render3d

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.reciprocal_lattice.phil_scope

show_rotation_axis = False
  .type = bool
show_beam_vector = False
  .type = bool
show_reciprocal_cell = True
  .type = bool
label_nearest_point = False
  .type = bool
marker_size = Auto
  .type = int(value_min=1)
autospin = False
  .type = bool
model_view_matrix = None
  .type = floats(size=16)
""",
    process_includes=True,
)


class OutlierStatus(Enum):
    all = 0
    inliers = 1
    outliers = 2


class ReflectionStatus(Enum):
    all = 0
    unindexed = 1
    indexed = 2
    integrated = 3


@magicgui(auto_call=True, d_min={"label": "high resolution", "step": 0.05})
def rlv_display(
    viewer: napari.Viewer,
    marker_size: int,
    d_min: float,
    z_min: int,
    z_max: int,
    outlier_status=OutlierStatus.all,
    reflection_status=ReflectionStatus.all,
):
    for layer in viewer.layers:
        if layer.name.startswith("relps"):
            layer.size[:] = marker_size
            shown = layer.properties["res"] >= d_min
            shown = shown & (layer.properties["z"] >= z_min)
            shown = shown & (layer.properties["z"] <= z_max)
            if outlier_status is OutlierStatus.inliers:
                shown = shown & ~layer.properties["outlier_status"]
            elif outlier_status is OutlierStatus.outliers:
                shown = shown & layer.properties["outlier_status"]
            if reflection_status is ReflectionStatus.unindexed:
                shown = shown & layer.properties["unindexed_status"]
            if reflection_status is ReflectionStatus.indexed:
                shown = shown & layer.properties["indexed_status"]
            if reflection_status is ReflectionStatus.integrated:
                shown = shown & layer.properties["integrated_status"]
            layer.shown = shown
            layer.refresh()


class ReciprocalLatticeViewer(Render3d):
    def __init__(self, parent, id, title, size, settings, napari_viewer, *args, **kwds):
        Render3d.__init__(self, settings=settings)

        self.rlv_window = RLVWindow(
            settings=self.settings,
        )

        self.napari_viewer = napari_viewer
        self._rlv_layers = {}

    @magicgui(auto_call=True)
    def rlv_geometry(self, invert_rotation_axis: bool, crystal_frame: bool):

        # Set values
        self.settings.reverse_phi = invert_rotation_axis
        self.settings.crystal_frame = crystal_frame

        self.update_settings()

        self.update_layers()

    def add_rlv_widgets(self):
        # Add the rlv_display widget and set values and limits
        self.napari_viewer.window.add_dock_widget(rlv_display, name="rlv display")
        rlv_display.marker_size.value = self.settings.marker_size
        rlv_display.d_min.value = self.settings.d_min
        rlv_display.d_min.min = self.settings.d_min
        rlv_display.z_min.min = floor(self.settings.z_min)
        rlv_display.z_min.max = ceil(self.settings.z_max)
        rlv_display.z_min.value = self.settings.z_min
        rlv_display.z_max.min = floor(self.settings.z_min)
        rlv_display.z_max.max = ceil(self.settings.z_max)
        rlv_display.z_max.value = self.settings.z_max
        if self.settings.outlier_display == "inliers":
            rlv_display.outlier_status.value = OutlierStatus.inliers
        if self.settings.outlier_display == "outliers":
            rlv_display.outlier_status.value = OutlierStatus.outliers

        # Add the rlv_geometry widget and set values
        self.napari_viewer.window.add_dock_widget(
            self.rlv_geometry, name="rlv geometry"
        )
        self.rlv_geometry.invert_rotation_axis.value = self.settings.reverse_phi
        self.rlv_geometry.crystal_frame.value = self.settings.crystal_frame

        # not ideal workaround https://forum.image.sc/t/magicgui-widget-spacing/66954/3?u=dagewa
        rlv_display.max_height = 200
        self.rlv_geometry.max_height = 200

        # Add the relp status display
        self.relp_status = Label(value="")
        container = Container(
            widgets=[
                self.relp_status,
            ]
        )
        self.napari_viewer.window.add_dock_widget(container, name="relp status")

    def update_layers(self):
        """Add or update the layers data"""

        # Get the reciprocal cells for drawing
        cells = self.rlv_window.draw_cells()

        # Add relps and cell (if present) as points and shapes layers for each id
        for exp_id in sorted(set(self.rlv_window.points_data["id"]), reverse=True):
            sel = self.rlv_window.points_data["id"] == exp_id

            # Convert points and colors to numpy arrays
            points = flumpy.to_numpy(self.rlv_window.points.select(sel))
            colors = flumpy.to_numpy(self.rlv_window.colors.select(sel))

            # Napari uses a left-handed coordinate system, so invert Z
            # https://forum.image.sc/t/3d-view-coordinate-system-is-left-handed/66995
            points[:, 2] *= -1

            id_data = self.rlv_window.points_data["id"].select(sel)
            x, y, z = self.rlv_window.points_data["xyz"].select(sel).parts()
            panel_data = self.rlv_window.points_data["panel"].select(sel)
            d_spacing_data = self.rlv_window.points_data["d_spacing"].select(sel)
            outlier_status = self.rlv_window.points_data["outlier_status"].select(sel)
            unindexed_status = self.rlv_window.points_data["unindexed_status"].select(
                sel
            )
            indexed_status = self.rlv_window.points_data["indexed_status"].select(sel)
            integrated_status = self.rlv_window.points_data["integrated_status"].select(
                sel
            )
            point_properties = {
                "id": flumpy.to_numpy(id_data),
                "x": flumpy.to_numpy(x).round(1),
                "y": flumpy.to_numpy(y).round(1),
                "z": flumpy.to_numpy(z).round(1),
                "panel": flumpy.to_numpy(panel_data),
                "res": flumpy.to_numpy(d_spacing_data).round(3),
                "outlier_status": flumpy.to_numpy(outlier_status),
                "unindexed_status": flumpy.to_numpy(unindexed_status),
                "indexed_status": flumpy.to_numpy(indexed_status),
                "integrated_status": flumpy.to_numpy(integrated_status),
            }

            if "miller_index" in self.rlv_window.points_data and exp_id != -1:
                h, k, l = (
                    self.rlv_window.points_data["miller_index"]
                    .select(sel)
                    .as_vec3_double()
                    .parts()
                )
                h = h.iround()
                k = k.iround()
                l = l.iround()
                point_properties["h"] = flumpy.to_numpy(h)
                point_properties["k"] = flumpy.to_numpy(k)
                point_properties["l"] = flumpy.to_numpy(l)

            layer_name = f"relps id: {exp_id}"
            if layer_name in self._rlv_layers:
                # Update existing layer
                relps_layer = self._rlv_layers[layer_name]
                relps_layer.data = points
                relps_layer.properties = point_properties
                relps_layer.refresh()
            else:
                # Add new layer
                relps_layer = self.napari_viewer.add_points(
                    points,
                    properties=point_properties,
                    face_color=colors,
                    size=self.settings.marker_size,
                    name=layer_name,
                )
                relps_layer.blending = "opaque"
                try:
                    relps_layer.canvas_size_limits = (2, 100)
                except AttributeError:
                    pass
                # https://forum.image.sc/t/adjusting-depth-fading/68181/4
                try:
                    # Remove except block when this fix is released
                    # https://github.com/napari/napari/issues/4683
                    relps_layer.antialiasing = 0.3
                except AttributeError:
                    relps_layer._antialias = 0.3

                self._rlv_layers[layer_name] = relps_layer

                @relps_layer.mouse_drag_callbacks.append
                def show_point_info(layer, event):
                    point_index = layer.get_value(
                        position=event.position,
                        view_direction=event.view_direction,
                        dims_displayed=event.dims_displayed,
                        world=True,
                    )
                    msg = [f"Active layer: {layer.name}"]
                    if point_index is not None:
                        x = layer.properties["x"][point_index]
                        y = layer.properties["y"][point_index]
                        z = layer.properties["z"][point_index]
                        h = layer.properties["h"][point_index]
                        k = layer.properties["k"][point_index]
                        l = layer.properties["l"][point_index]
                        outlier = layer.properties["outlier_status"][point_index]
                        panel = layer.properties["panel"][point_index]
                        d_min = layer.properties["res"][point_index]

                        msg.append(f"panel: {panel}")
                        msg.append(f"xyz: {x:.1f} {y:.1f} {z:.1f}")
                        msg.append(f"res: {d_min:.2f} Å")
                        msg.append(f"outlier: {outlier}")
                        if layer.properties["indexed_status"][point_index]:
                            msg.append(f"hkl: ({h} {k} {l})")
                        msg = "\n".join(msg)
                        self.relp_status.value = msg

            # Add the cell as a shapes layer, if it exists
            cell = cells.get(exp_id)
            if cell:
                layer_name = f"cell id: {exp_id}"
                if layer_name in self._rlv_layers:
                    # Update existing layer
                    cell_layer = self._rlv_layers[layer_name]
                    cell_layer.data = cell.lines
                    cell_layer.refresh()
                else:
                    # Create new layer
                    labels = ["a*", "b*", "c*"] + [""] * 9
                    properties = {"label": labels}
                    text_parameters = {
                        "text": "label",
                        "size": 14,
                        "color": np.array(cell.colors),
                        "anchor": "center",
                        "translation": -5,
                    }
                    # Set "translation": [0, -5, 0], when https://github.com/napari/napari/pull/5321 is released
                    cell_layer = self.napari_viewer.add_shapes(
                        cell.lines,
                        shape_type="line",
                        edge_width=2,
                        edge_color=np.array(cell.colors),
                        name=layer_name,
                        properties=properties,
                        text=text_parameters,
                    )
                    self._rlv_layers[layer_name] = cell_layer

                # Link the visibility of the relps and cell layer
                link_layers([relps_layer, cell_layer], ("visible",))

        # Determine suitable scale factor for lines. Code extracted from draw_axis
        if self.rlv_window.minimum_covering_sphere is None:
            self.rlv_window.update_minimum_covering_sphere()
        s = self.rlv_window.minimum_covering_sphere
        scale = max(max(s.box_max()), abs(min(s.box_min())))

        # Calculate rotation axis line, taking into account Napari's left-handed coordinate system
        axis = self.rlv_window.rotation_axis
        axis_line = np.array(
            [[0, 0, 0], [axis[0] * scale, axis[1] * scale, -1 * axis[2] * scale]]
        )

        # Draw the rotation axis
        if "axis" in self._rlv_layers:
            # Update existing layer
            axis_layer = self._rlv_layers["axis"]
            for i, (_, b) in enumerate(zip(axis_layer.data, axis_line)):
                axis_layer.data[i] = b
            axis_layer.refresh()
        else:
            # Create new layer
            labels = [
                "phi",
            ]
            properties = {"label": labels}
            text_parameters = {
                "text": "label",
                "size": 14,
                "color": "white",
                "anchor": "center",
                "translation": -20,
            }
            # Set "translation": [0, -20, 0], when https://github.com/napari/napari/pull/5321 is released
            axis_layer = self.napari_viewer.add_shapes(
                [
                    axis_line,
                ],
                shape_type="line",
                edge_width=2,
                edge_color="white",
                name="axis",
                properties=properties,
                text=text_parameters,
            )
            self._rlv_layers["axis"] = axis_layer

        # Calculate the beam vector line, taking into account Napari's left-handed coordinate system
        beam_vector = self.rlv_window.beam_vector
        beam_vector_line = np.array(
            [
                [0, 0, 0],
                [
                    beam_vector[0] * scale,
                    beam_vector[1] * scale,
                    -1 * beam_vector[2] * scale,
                ],
            ]
        )

        # Draw the beam vector
        if "beam_vector" in self._rlv_layers:
            # Update existing layer
            beam_vector_layer = self._rlv_layers["beam_vector"]
            for i, (_, b) in enumerate(zip(beam_vector_layer.data, beam_vector_line)):
                beam_vector_layer.data[i] = b
            beam_vector_layer.refresh()
        else:
            # Create new layer
            labels = [
                "beam",
            ]
            properties = {"label": labels}
            text_parameters = {
                "text": "label",
                "size": 14,
                "color": "white",
                "anchor": "center",
                "translation": -20,
            }
            # Set "translation": [0, -20, 0], when https://github.com/napari/napari/pull/5321 is released
            beam_vector_layer = self.napari_viewer.add_shapes(
                [
                    beam_vector_line,
                ],
                shape_type="line",
                edge_width=2,
                edge_color="white",
                name="beam_vector",
                properties=properties,
                text=text_parameters,
            )
            self._rlv_layers["beam_vector"] = beam_vector_layer

        return

    def load_models(self, experiments, reflections):
        Render3d.load_models(self, experiments, reflections)
        if self.settings.beam_centre is not None:

            pass
        if self.settings.marker_size is Auto:
            max_radius = max(self.reflections["rlp"].norms())
            volume = 4 / 3 * pi * max_radius**3
            density = len(self.reflections) / volume
            # Set marker size depending on relp density, where
            # 1000 < density < 20000 ==> max_size < marker_size < min_size
            # XXX this does not take into account narrow wedges!
            min_size, max_size = 1, 10
            grad = (max_size - min_size) / (20000 - 1000)
            intercept = max_size - 1000 * grad
            marker_size = grad * density + intercept
            marker_size = max(marker_size, min_size)
            marker_size = min(marker_size, max_size)
            self.settings.marker_size = marker_size

    def set_points(self):
        Render3d.set_points(self)

    def update_settings(self, *args, **kwds):
        self.set_beam_centre(self.settings.beam_centre_panel, self.settings.beam_centre)
        self.map_points_to_reciprocal_space()
        self.set_points()


class RLVWindow:
    def __init__(self, settings, *args, **kwds):
        self.settings = settings
        self.points = flex.vec3_double()
        self.colors = None
        self.palette = None
        self.rotation_axis = None
        self.beam_vector = None
        self.recip_latt_vectors = None
        self.recip_crystal_vectors = None
        self.flag_show_minimum_covering_sphere = False
        self.minimum_covering_sphere = None
        self.field_of_view_y = 0.001

    def set_points(self, points):
        self.points = points
        self.points_display_list = None
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()

    def set_points_data(self, reflections):
        dstar = reflections["rlp"].norms()
        dstar.set_selected(dstar == 0, 1e-8)
        self.points_data = {
            "panel": reflections["panel"],
            "id": reflections["id"],
            "xyz": reflections["xyzobs.px.value"],
            "d_spacing": 1 / dstar,
            "outlier_status": reflections.get_flags(reflections.flags.centroid_outlier),
            "unindexed_status": reflections.get_flags(reflections.flags.strong)
            & ~reflections.get_flags(reflections.flags.indexed),
            "indexed_status": reflections.get_flags(reflections.flags.indexed),
            "integrated_status": reflections.get_flags(reflections.flags.integrated),
        }
        if "miller_index" in reflections:
            self.points_data["miller_index"] = reflections["miller_index"]

    def set_colors(self, colors):
        assert len(colors) == len(self.points)
        self.colors = colors

    def set_palette(self, palette):
        self.palette = palette

    def set_rotation_axis(self, axis):
        self.rotation_axis = axis

    def set_beam_vector(self, beam):
        self.beam_vector = beam

    def set_reciprocal_lattice_vectors(self, vectors_per_crystal):
        self.recip_latt_vectors = vectors_per_crystal

    def set_reciprocal_crystal_vectors(self, vectors_per_crystal):
        self.recip_crystal_vectors = vectors_per_crystal

    def update_minimum_covering_sphere(self):
        n_points = min(1000, self.points.size())
        isel = flex.random_permutation(self.points.size())[:n_points]
        self.minimum_covering_sphere = minimum_covering_sphere(self.points.select(isel))

    def draw_cells(self):
        """Create a dictionary mapping experiment id to lists of lines and
        colours for reciprocal unit cells"""

        result = {}
        CellDrawing = namedtuple("CellDrawing", ["lines", "colors"])
        if self.settings.show_reciprocal_cell:
            # if we don't have one sort of vector we don't have the other either
            vectors = self.recip_latt_vectors
            if self.settings.crystal_frame:
                vectors = self.recip_crystal_vectors

            if vectors:
                for i, axes in enumerate(vectors):
                    if self.settings.experiment_ids:
                        if i not in self.settings.experiment_ids:
                            continue
                    j = (i + 1) % self.palette.size()
                    result[i] = CellDrawing(
                        lines=self.cell_edges(axes), colors=[self.palette[j]] * 12
                    )

        return result

    def cell_edges(self, axes):
        astar, bstar, cstar = axes[0], axes[1], axes[2]

        # Invert Z to account for napari's LH coordinate system
        M = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, -1))
        astar = M * astar
        bstar = M * bstar
        cstar = M * cstar

        farpoint = astar + bstar + cstar

        lines = [
            np.array([[0, 0, 0], [*astar.elems]]),
            np.array([[0, 0, 0], [*bstar.elems]]),
            np.array([[0, 0, 0], [*cstar.elems]]),
            np.array([[*astar.elems], [*(astar + bstar).elems]]),
            np.array([[*astar.elems], [*(astar + cstar).elems]]),
            np.array([[*bstar.elems], [*(bstar + astar).elems]]),
            np.array([[*bstar.elems], [*(bstar + cstar).elems]]),
            np.array([[*cstar.elems], [*(cstar + astar).elems]]),
            np.array([[*cstar.elems], [*(cstar + bstar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - astar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - bstar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - cstar).elems]]),
        ]
        return lines
