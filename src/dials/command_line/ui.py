# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_SET_DISPATCHER_NAME dials.ui

"""
ui.py — A step-by-step graphical front end for the DIALS
macromolecular crystallography data processing suite.

This GUI wraps the command-line DIALS programs used in the
"Processing in Detail" workflow tutorials (CCP4/DLS & CCP4/APS 2024
workshops, see https://github.com/graeme-winter/dials_tutorials):

    dials.import
    dials.find_spots
    dials.search_beam_position   (optional)
    dials.index
    dials.refine_bravais_settings (optional)
    dials.refine
    dials.integrate
    dials.symmetry
    dials.scale
    dials.merge / dials.export

as well as the interactive viewing tools:

    dials.show
    dials.image_viewer
    dials.reciprocal_lattice_viewer
    dials.report

The GUI does not re-implement any DIALS functionality itself: it only
constructs the correct command line for each stage, runs it as a
subprocess in a chosen working directory, streams the live output to
the screen, and then reads back the DIALS-generated `dials.<program>.log`
file to build a short, readable digest of the result (RMSDs, % indexed,
resolution, space group, merging statistics, etc). Every processing
stage and viewer described in the tutorial WORKFLOW.md is reachable
from this interface, and the intermediate .expt / .refl files can be
freely substituted at every stage, so any of the tutorials on that page
(the basic single-sweep workflow, the "if you're impatient" TLDR path,
optimising the beam centre, forcing a Bravais lattice/space group, and
scaling anomalous vs native data) can be worked through interactively.

Requirements
------------
Python 3.8+ with wxPython (install with `pip install wxPython`, or
`libtbx.pip install wxPython` inside a DIALS/cctbx environment — DIALS
already ships wxPython for its own viewers, so it is usually present in a
sourced DIALS environment). DIALS itself must already be installed and
set up in the environment the GUI is launched from (i.e. `dials.import`
etc. must be on $PATH) - this GUI is a front end, not a replacement, for
that installation.

Run with:

    python3 dials_gui.py

"""

from __future__ import annotations

import base64
import glob
import io
import json
import math
import os
import queue
import re
import shutil
import subprocess
import threading
import webbrowser
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import wx

# matplotlib is an optional dependency: the live-plotting "Plots" tab is
# only offered if it (and its wxAgg backend) import successfully. The rest
# of the GUI — the whole pipeline, log summaries, dials.report button —
# works without it, so a missing matplotlib degrades gracefully to "no
# Plots tab" rather than failing to start.
try:
    import matplotlib

    matplotlib.use("WXAgg")
    from matplotlib.backends.backend_wxagg import (
        FigureCanvasWxAgg as FigureCanvas,
    )
    from matplotlib.backends.backend_wxagg import (
        NavigationToolbar2WxAgg as NavigationToolbar,
    )
    from matplotlib.figure import Figure

    HAVE_MPL = True
except Exception:  # pragma: no cover - depends on environment
    HAVE_MPL = False


# --------------------------------------------------------------------------
# Application / taskbar icon
#
# Embedded as base64-encoded PNG data so the GUI is a single self-contained
# file with no external icon asset to lose track of. Decoded lazily (once)
# into a wx.Icon the first time a frame asks for it.
# --------------------------------------------------------------------------

APP_ICON_PNG_BASE64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAGwAAABrCAYAAACSY2d1AAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQA"
    "APoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAhGVYSWZNTQAqAAAACAAFARIAAwAAAAEAAQAAARoABQAAAAEA"
    "AABKARsABQAAAAEAAABSASgAAwAAAAEAAgAAh2kABAAAAAEAAABaAAAAAAAAAJYAAAABAAAAlgAAAAEAA6ABAAMA"
    "AAABAAEAAKACAAQAAAABAAAAbKADAAQAAAABAAAAawAAAAA5rXDkAAAACXBIWXMAABcSAAAXEgFnn9JSAAABWWlU"
    "WHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0"
    "az0iWE1QIENvcmUgNi4wLjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkv"
    "MDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAg"
    "ICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOk9y"
    "aWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgIDwvcmRmOkRlc2NyaXB0aW9uPgogICA8L3JkZjpS"
    "REY+CjwveDp4bXBtZXRhPgoZXuEHAAAiSElEQVR4AdVdaYxk11U+r+rV2lW9Ts/usWe8jMdLbBNjRyGYiACyQAoS"
    "RPwjEj/4g4AfSPxA8AMk+MMfFJGIPyCCUIQCIggJCYEVoQSSOBln4vEae+LZPXtPT3d1Vdf++L7z3qt+VfWW+6qX"
    "qbnq6np1313Puefcc88991zrqS/+ueOII37o49fGxqbcvHFH+lvRYiHeD3zyP5lMRo4cXpZCIe+/Nv52vDIzliWC"
    "j//buID7mLDT7cmtW3dls7GpbTdtip3Nyv79i1IuFyOzEA53V9dlFZ/RYPsRigAk7AFLtVpjDFk+gpiez8GQBtDE"
    "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"
)


_app_icon_cache: Optional["wx.Icon"] = None


def get_app_icon() -> "wx.Icon":
    """Return the application icon, decoding the embedded PNG on first use."""
    global _app_icon_cache
    if _app_icon_cache is None:
        png_bytes = base64.b64decode(APP_ICON_PNG_BASE64)
        image = wx.Image(io.BytesIO(png_bytes), wx.BITMAP_TYPE_PNG)
        bitmap = wx.Bitmap(image)
        icon = wx.Icon()
        icon.CopyFromBitmap(bitmap)
        _app_icon_cache = icon
    return _app_icon_cache


# --------------------------------------------------------------------------
# Step / field definitions
# --------------------------------------------------------------------------


@dataclass
class ExtraField:
    """A single extra command-line parameter exposed in the GUI."""

    key: str
    label: str
    kind: str = "entry"  # "entry", "check", "combo"
    default: str = ""
    choices: Optional[List[str]] = None
    help: str = ""
    # For "check" fields: the value emitted when ticked. Defaults to
    # "True" (so the arg is `key=True`); set to e.g. "false" for a toggle
    # like joint=false that should emit `joint=false` when ticked.
    check_value: str = "True"

    def build_arg(self, value) -> Optional[str]:
        if self.kind == "check":
            if value:
                return f"{self.key}={self.check_value}"
            return None
        value = (value or "").strip()
        if not value:
            return None
        return f"{self.key}={value}"


@dataclass
class InputSpec:
    label: str
    default: str


@dataclass
class StepDef:
    id: str
    title: str
    program: str
    help: str
    inputs: List[InputSpec]
    extra_fields: List[ExtraField] = field(default_factory=list)
    outputs: List[str] = field(default_factory=list)
    log_file: str = ""
    optional: bool = False
    is_import: bool = False
    # Some steps (merge/export) choose their program dynamically.
    dynamic: bool = False
    # Which live-plot view (if any) this step supports. One of
    # "find_spots", "refine", "integrate", "scale" — or "" for steps with
    # no Plots tab. Drives whether select_step() adds a "Plots" tab and
    # which _refresh_plots_* method it uses.
    plot_kind: str = ""


STEPS: List[StepDef] = [
    StepDef(
        id="import",
        title="1. Import",
        program="dials.import",
        help=(
            "Read image headers and write imported.expt describing the "
            "experiment geometry (detector, beam, goniometer, scan). "
            "Nothing else happens at this stage - if the beam centre, "
            "distance or wavelength look wrong here, fix it before going "
            "further."
        ),
        inputs=[],
        is_import=True,
        extra_fields=[
            ExtraField(
                "image_range",
                "Image range (start,end)",
                "entry",
                help="e.g. 1,1200 - leave blank to use all images",
            ),
        ],
        outputs=["imported.expt"],
        log_file="dials.import.log",
    ),
    StepDef(
        id="find_spots",
        title="2. Find Spots",
        program="dials.find_spots",
        help=(
            "Scan every image for strong reflections and write strong.refl. "
            "This is one of the two most time-consuming steps as every "
            "image is read and processed."
        ),
        inputs=[InputSpec("Experiment file", "imported.expt")],
        extra_fields=[
            ExtraField("nproc", "nproc (blank = all cores)", "entry"),
        ],
        outputs=["strong.refl"],
        log_file="dials.find_spots.log",
        plot_kind="find_spots",
    ),
    StepDef(
        id="search_beam",
        title="3. Search Beam Position (optional)",
        program="dials.search_beam_position",
        help=(
            "Refine the beam centre against the found spots, writing "
            "optimised.expt. Useful if the reciprocal lattice viewer shows "
            "spots that don't line up. The shift reported should usually "
            "be small for a well-calibrated beamline."
        ),
        inputs=[
            InputSpec("Experiment file", "imported.expt"),
            InputSpec("Reflection file", "strong.refl"),
        ],
        outputs=["optimised.expt"],
        log_file="dials.search_beam_position.log",
        optional=True,
    ),
    StepDef(
        id="index",
        title="4. Index",
        program="dials.index",
        help=(
            "Identify the crystal lattice from the strong spot positions "
            "and assign Miller indices, writing indexed.expt / "
            "indexed.refl. Use the experiment file from Import, or from "
            "Search Beam Position if you ran it. Set space_group / "
            "unit_cell here if known, or after inspecting Bravais lattice "
            "options below. For MULTIPLE crystals (many imported sweeps "
            "that do not share an orientation matrix), tick 'multi-crystal "
            "(joint=false)' so each sweep is indexed independently in one "
            "run - see the Cows/Pigs/People workflow."
        ),
        inputs=[
            InputSpec("Experiment file", "imported.expt"),
            InputSpec("Reflection file", "strong.refl"),
        ],
        extra_fields=[
            ExtraField(
                "joint",
                "multi-crystal (joint=false)",
                "check",
                check_value="false",
                help="index many crystals independently in one run",
            ),
            ExtraField("space_group", "space_group", "entry"),
            ExtraField(
                "unit_cell", "unit_cell", "entry", help="e.g. 78,78,78,90,90,90"
            ),
            ExtraField(
                "max_lattices",
                "max_lattices",
                "entry",
                help="set to 2 (say) if % indexed is low and a "
                "second lattice is suspected",
            ),
        ],
        outputs=["indexed.expt", "indexed.refl"],
        log_file="dials.index.log",
        plot_kind="index",
    ),
    StepDef(
        id="bravais",
        title="5. Bravais Lattice Determination (optional)",
        program="dials.refine_bravais_settings",
        help=(
            "List every Bravais lattice approximately consistent with the "
            "triclinic cell from indexing, with the RMS deviation each "
            "would introduce. If a solution's rmsd is not noticeably worse "
            "than the triclinic one, its lattice / space group is a good "
            "candidate. The simplest way to use a solution is to re-run "
            "Index above with space_group set accordingly (rather than "
            "using bravais_setting_N.expt directly, which needs the "
            "reflections re-indexed to match)."
        ),
        inputs=[
            InputSpec("Experiment file", "indexed.expt"),
            InputSpec("Reflection file", "indexed.refl"),
        ],
        outputs=[],
        log_file="dials.refine_bravais_settings.log",
        optional=True,
    ),
    StepDef(
        id="refine",
        title="6. Refine",
        program="dials.refine",
        help=(
            "Re-refine the crystal/detector/beam models, including "
            "scan-varying refinement of the crystal, writing refined.expt "
            "/ refined.refl. RMSDs should improve slightly relative to "
            "the end of indexing."
        ),
        inputs=[
            InputSpec("Experiment file", "indexed.expt"),
            InputSpec("Reflection file", "indexed.refl"),
        ],
        outputs=["refined.expt", "refined.refl"],
        log_file="dials.refine.log",
        plot_kind="refine",
    ),
    StepDef(
        id="integrate",
        title="7. Integrate",
        program="dials.integrate",
        help=(
            "Build a reflection profile model and integrate the "
            "background-subtracted intensity of every predicted "
            "reflection, writing integrated.expt / integrated.refl. This "
            "is the most computationally expensive step."
        ),
        inputs=[
            InputSpec("Experiment file", "refined.expt"),
            InputSpec("Reflection file", "refined.refl"),
        ],
        extra_fields=[
            ExtraField(
                "prediction.d_min",
                "prediction.d_min",
                "entry",
                help="optional resolution limit, e.g. 1.8",
            ),
            ExtraField("nproc", "nproc (blank = all cores)", "entry"),
        ],
        outputs=["integrated.expt", "integrated.refl"],
        log_file="dials.integrate.log",
        plot_kind="integrate",
    ),
    StepDef(
        id="symmetry",
        title="8. Symmetry (single crystal)",
        program="dials.symmetry",
        help=(
            "Assess spot positions and intensities to identify symmetry "
            "operations present in the data, compose these into a "
            "candidate Laue group / space group, and write "
            "symmetrized.expt / symmetrized.refl. Not needed if the "
            "correct space group was already set at Index. For MULTIPLE "
            "crystals use Cosym (below) instead - it determines symmetry "
            "and resolves indexing ambiguity across all data sets at once."
        ),
        inputs=[
            InputSpec("Experiment file", "integrated.expt"),
            InputSpec("Reflection file", "integrated.refl"),
        ],
        outputs=["symmetrized.expt", "symmetrized.refl"],
        log_file="dials.symmetry.log",
    ),
    StepDef(
        id="cosym",
        title="8b. Cosym (multi-crystal)",
        program="dials.cosym",
        help=(
            "For MULTIPLE crystals: determine the Patterson symmetry AND "
            "resolve indexing ambiguity across all data sets simultaneously "
            "(replaces dials.symmetry). Aligns the lattices in reciprocal "
            "space, estimates the crystal symmetry, and writes "
            "symmetrized.expt / symmetrized.refl plus dials.cosym.html. Run "
            "this instead of Symmetry when you indexed with joint=false."
        ),
        inputs=[
            InputSpec("Experiment file", "integrated.expt"),
            InputSpec("Reflection file", "integrated.refl"),
        ],
        outputs=["symmetrized.expt", "symmetrized.refl"],
        log_file="dials.cosym.log",
        optional=True,
    ),
    StepDef(
        id="correlation_matrix",
        title="8c. Correlation Matrix (multi-crystal)",
        program="dials.correlation_matrix",
        help=(
            "For MULTIPLE crystals: measure the pairwise similarity of the "
            "data sets and cluster the isomorphous ones using the OPTICS "
            "algorithm, writing dials.correlation_matrix.html (with the "
            "correlation / cos-angle matrices, dendrograms and cluster "
            "assignments). Normally run after Cosym on symmetrized data; you "
            "can ALSO run it after Scale by ticking 'use scaled data' (which "
            "switches the inputs to scaled.expt/.refl and writes to "
            "dials.correlation_matrix.scaled.html so it doesn't overwrite the "
            "earlier run). 'output clusters' (on by default) writes "
            "cluster_0.expt/.refl, cluster_1.expt/.refl, ... which can then "
            "be scaled independently (see Scale). The Plots tab visualises "
            "the matrices, reachability and cluster coordinates from the HTML."
        ),
        inputs=[
            InputSpec("Experiment file", "symmetrized.expt"),
            InputSpec("Reflection file", "symmetrized.refl"),
        ],
        extra_fields=[
            ExtraField(
                "use_scaled",
                "use scaled data (run after scaling)",
                "check",
                help="switch inputs to scaled.expt/.refl and write to "
                "dials.correlation_matrix.scaled.html",
            ),
            ExtraField(
                "significant_clusters.output",
                "output clusters (write cluster_N.expt/.refl)",
                "check",
                default="True",
                help="on by default; needed to scale clusters separately",
            ),
        ],
        outputs=[],
        log_file="dials.correlation_matrix.log",
        optional=True,
        plot_kind="correlation_matrix",
    ),
    StepDef(
        id="scale",
        title="9. Scale",
        program="dials.scale",
        help=(
            "Correct for radiation damage, beam intensity changes and "
            "sample absorption, writing scaled.expt / scaled.refl and "
            "dials.scale.html. Merging statistics and the error model are "
            "printed at the end - this is where you find out about the "
            "final quality of the data. Tick 'anomalous' for anomalous "
            "data (e.g. SAD/MAD). For MULTIPLE clusters from the "
            "Correlation Matrix step, use the cluster selector below to "
            "scale each cluster_N independently - each run writes its own "
            "scaled_cluster_N.* / dials.scale.cluster_N.* so nothing is "
            "overwritten."
        ),
        inputs=[
            InputSpec("Experiment file", "symmetrized.expt"),
            InputSpec("Reflection file", "symmetrized.refl"),
        ],
        extra_fields=[
            ExtraField("anomalous", "anomalous", "check"),
            ExtraField(
                "absorption_level",
                "absorption_level",
                "combo",
                choices=["", "low", "medium", "high"],
                help="low (~1%, default), medium (~5%), high (~25%)",
            ),
            ExtraField(
                "d_min",
                "d_min",
                "entry",
                help="optional resolution cutoff from CC-half fit",
            ),
        ],
        outputs=["scaled.expt", "scaled.refl", "dials.scale.html"],
        log_file="dials.scale.log",
        plot_kind="scale",
    ),
    StepDef(
        id="merge_export",
        title="10. Merge / Export",
        program="dials.merge",
        help=(
            "Produce a final MTZ file for downstream use. 'merge' writes "
            "a scaled and merged MTZ (dials.merge) - use this for most "
            "structure solution / molecular replacement pipelines. "
            "'export' writes the scaled but unmerged MTZ (dials.export)."
        ),
        inputs=[
            InputSpec("Experiment file", "scaled.expt"),
            InputSpec("Reflection file", "scaled.refl"),
        ],
        extra_fields=[
            ExtraField(
                "mode", "mode", "combo", default="merge", choices=["merge", "export"]
            ),
            ExtraField(
                "d_min",
                "d_min",
                "entry",
                help="optional resolution cutoff suggested by scaling",
            ),
        ],
        outputs=[],
        log_file="",
        dynamic=True,
    ),
]


TOOLS = [
    ("dials.show", "dials.show", ["Experiment / reflection file(s)"]),
    (
        "dials.image_viewer",
        "dials.image_viewer",
        ["Experiment file", "Reflection file (optional)"],
    ),
    (
        "dials.reciprocal_lattice_viewer",
        "dials.reciprocal_lattice_viewer",
        ["Experiment file", "Reflection file"],
    ),
    ("dials.report", "dials.report", ["Experiment file", "Reflection file"]),
]


# --------------------------------------------------------------------------
# Log summarisation helpers
# --------------------------------------------------------------------------

_TABLE_LINE = re.compile(r"^\s*[+|].*[+|]\s*$")
_KEYWORD_LINES = re.compile(
    r"(Best solution|Space group|Unit cell|Laue group|% indexed|"
    r"num images|num stills|sequences|Writing experiments|Writing "
    r"reflections|reflections indexed|resolution|RMSD|Suggested|"
    r"High resolution limit|Low resolution limit|Completeness|"
    r"Multiplicity|I/sigma|CC half|Rmerge|Rmeas|Rpim|Error model|"
    r"estimated I/sigma|shift|Reindex operator)",
    re.IGNORECASE,
)


def summarise_log(text: str, max_blocks: int = 8) -> str:
    """Produce a short, readable digest of a DIALS log file."""

    if not text.strip():
        return "(no log output captured yet)"

    lines = text.splitlines()
    out: List[str] = []
    seen_keyword_lines = set()
    blocks_found = 0
    i = 0
    n = len(lines)
    while i < n and blocks_found < max_blocks:
        line = lines[i]
        if _TABLE_LINE.match(line):
            # capture the whole table block
            block = []
            # walk backwards to include a header line above a leading '+--' row
            j = i
            while (
                j > 0
                and lines[j - 1].strip()
                and not _TABLE_LINE.match(lines[j - 1])
                and len(block) < 1
            ):
                j -= 1
            start = max(i - 2, 0)
            k = i
            while k < n and (
                _TABLE_LINE.match(lines[k]) or "|" in lines[k] or not lines[k].strip()
            ):
                block.append(lines[k])
                k += 1
                if k - i > 40:
                    break
            out.extend(lines[start:i])
            out.extend(block)
            out.append("")
            blocks_found += 1
            i = k
            continue
        if _KEYWORD_LINES.search(line) and line.strip() not in seen_keyword_lines:
            out.append(line.rstrip())
            seen_keyword_lines.add(line.strip())
        i += 1

    if not out:
        # fall back to the tail of the log
        out = ["(no recognised summary patterns found - showing tail of log)", ""]
        out.extend(lines[-40:])

    return "\n".join(out)


# --------------------------------------------------------------------------
# Live-plot data extraction
# --------------------------------------------------------------------------
#
# The functions below turn a step's cumulative stdout / log text into
# structured numeric series for the "Plots" tab. They are all pure
# functions (no Tk / matplotlib dependency) and, crucially, tolerant of
# *partial* input: they can be called repeatedly while a step is still
# running so the plots update live, parsing whatever complete data exists
# so far and ignoring the rest.


def _split_table_row(line: str) -> Optional[List[str]]:
    """Split a '| a | b | c |' pipe-table row into its cell strings, or
    return None if this line isn't a data-bearing pipe row (e.g. it's a
    +----+ / |----| border, or not a table row at all)."""
    s = line.strip()
    if not s.startswith("|"):
        return None
    if set(s) <= set("|+-= "):
        return None
    return [c.strip() for c in s.strip("|").split("|")]


_FIND_SPOTS_RE = re.compile(
    r"Found\s+(\d+)\s+strong\s+pixels\s+on\s+image\s+(\d+)", re.IGNORECASE
)


def parse_find_spots(text: str) -> Dict[str, List[float]]:
    """find_spots: 'Found N strong pixels on image M' -> {image, pixels}
    sorted by image number."""
    by_image: Dict[int, int] = {}
    for m in _FIND_SPOTS_RE.finditer(text):
        by_image[int(m.group(2))] = int(m.group(1))
    images = sorted(by_image)
    return {
        "image": [float(i) for i in images],
        "pixels": [float(by_image[i]) for i in images],
    }


# find_spots processes one imageset at a time; each is introduced by a
# banner block like:
#   --------------------------------------------------------------------
#   Finding strong spots in imageset 0
#   --------------------------------------------------------------------
# and then emits its own "Found N strong pixels on image M" lines (with M
# a per-imageset frame number, restarting at 1 for each imageset). We split
# on these banners so each imageset is a separate series on the plot,
# captioned by its imageset number.
#
# The real DIALS wording is "Finding strong spots IN imageset N" (see
# dials.algorithms.spot_finding.finder). We accept "in" or "on", and allow
# the phrase to appear with or without the word "strong", to be robust to
# small wording changes across DIALS versions.
_FIND_SPOTS_IMAGESET_RE = re.compile(
    r"Finding\s+(?:strong\s+)?spots\s+(?:in|on)\s+imageset\s+(\d+)",
    re.IGNORECASE,
)


def parse_find_spots_by_imageset(text: str) -> List[Dict[str, object]]:
    """Split find_spots output into per-imageset series.

    Returns a list of {'imageset': int, 'image': [...], 'pixels': [...]}
    in the order the imagesets appear. Any 'Found ...' lines that occur
    before the first banner (or if there are no banners at all) are
    collected under imageset None as a single fallback series, so a
    single-sweep run - which has no banner - still plots.

    Tolerant of streaming: the currently-processing (last) imageset simply
    has fewer points until it finishes."""
    lines = text.splitlines()
    series: List[Dict[str, object]] = []
    current: Optional[Dict[str, object]] = None

    def _new(imageset: Optional[int]) -> Dict[str, object]:
        d: Dict[str, object] = {"imageset": imageset, "image": [], "pixels": []}
        series.append(d)
        return d

    for line in lines:
        banner = _FIND_SPOTS_IMAGESET_RE.search(line)
        if banner:
            current = _new(int(banner.group(1)))
            continue
        m = _FIND_SPOTS_RE.search(line)
        if m:
            if current is None:
                current = _new(None)
            current["image"].append(float(m.group(2)))  # type: ignore[union-attr]
            current["pixels"].append(float(m.group(1)))  # type: ignore[union-attr]

    # Drop any empty banner-only series (e.g. a banner seen but no spot
    # lines yet is fine to keep; but a trailing empty one adds nothing).
    return [s for s in series if s["image"]] or series


_REFINE_HEADER_RE = re.compile(r"Refinement steps", re.IGNORECASE)


def parse_refine_steps(text: str) -> Dict[str, List[float]]:
    """refine: the 'Refinement steps' table -> {step, rmsd_x, rmsd_y,
    rmsd_phi}. Reads rows after the *last* 'Refinement steps' header seen
    (there can be more than one macrocycle)."""
    lines = text.splitlines()
    steps: List[float] = []
    rmsd_x: List[float] = []
    rmsd_y: List[float] = []
    rmsd_phi: List[float] = []

    header_idx = None
    for i, line in enumerate(lines):
        if _REFINE_HEADER_RE.search(line):
            header_idx = i
    if header_idx is None:
        return {"step": [], "rmsd_x": [], "rmsd_y": [], "rmsd_phi": []}

    for line in lines[header_idx + 1 :]:
        cells = _split_table_row(line)
        if cells is None:
            if steps and not line.strip():
                break
            continue
        if len(cells) < 5:
            continue
        try:
            step = float(int(cells[0]))
            x = float(cells[2])
            y = float(cells[3])
            phi = float(cells[4])
        except ValueError:
            continue
        steps.append(step)
        rmsd_x.append(x)
        rmsd_y.append(y)
        rmsd_phi.append(phi)

    return {"step": steps, "rmsd_x": rmsd_x, "rmsd_y": rmsd_y, "rmsd_phi": rmsd_phi}


_REFINE_GROUP_RE = re.compile(
    r"Selected group of experiments to refine with original ids:\s*([0-9,\s]+)",
    re.IGNORECASE,
)


def _parse_one_refine_table(
    lines: List[str], start: int
) -> Tuple[Dict[str, List[float]], int]:
    """Parse a single 'Refinement steps' table whose header is at/after
    `start`. Returns (table_dict, index_after_table). table_dict is empty
    if no data rows were found."""
    tbl = {"step": [], "rmsd_x": [], "rmsd_y": [], "rmsd_phi": []}
    j = start
    n = len(lines)
    started = False
    while j < n:
        cells = _split_table_row(lines[j])
        if cells is None:
            if started and not lines[j].strip():
                break
            j += 1
            continue
        if len(cells) < 5:
            j += 1
            continue
        try:
            step = float(int(cells[0]))
            x = float(cells[2])
            y = float(cells[3])
            phi = float(cells[4])
        except ValueError:
            j += 1
            continue
        started = True
        tbl["step"].append(step)
        tbl["rmsd_x"].append(x)
        tbl["rmsd_y"].append(y)
        tbl["rmsd_phi"].append(phi)
        j += 1
    return tbl, j


def parse_all_refine_steps(text: str) -> List[Dict[str, List[float]]]:
    """refine (multi-crystal joint=false): one convergence table per
    refinement RUN, in order.

    A run is delimited by a
    'Selected group of experiments to refine with original ids: <ids>'
    line. Within each run, refinement may print several 'Refinement steps'
    tables (e.g. static then scan-varying macrocycles); we keep only the
    LAST one in the run (the final convergence). This avoids the earlier
    bug where every macrocycle table counted as a separate 'run', inflating
    the run count well beyond the number of experiments actually refined.

    Each returned dict is {'step','rmsd_x','rmsd_y','rmsd_phi','ids'} where
    'ids' is the original experiment id string from the group marker (or ''
    if the run wasn't introduced by a marker - e.g. single-crystal refine,
    which has no such marker and yields a single run from its lone table).

    RMSD columns may be mm (single sweep) or px (multi-crystal); we take
    columns 2/3/4 after the integer step in column 0."""
    lines = text.splitlines()
    n = len(lines)

    # Find the group-marker positions.
    markers = [
        (i, m.group(1).strip())
        for i, line in enumerate(lines)
        for m in [_REFINE_GROUP_RE.search(line)]
        if m
    ]

    tables: List[Dict[str, List[float]]] = []

    if not markers:
        # No group markers (single-crystal refine, or an older/661 DIALS
        # format): fall back to keeping the LAST 'Refinement steps' table in
        # the whole output, so we don't over-count macrocycles. (Previously
        # this returned every table; the final one is the meaningful
        # convergence for a single run.)
        last = None
        i = 0
        while i < n:
            if _REFINE_HEADER_RE.search(lines[i]):
                tbl, j = _parse_one_refine_table(lines, i + 1)
                if tbl["step"]:
                    last = tbl
                i = max(j, i + 1)
            else:
                i += 1
        if last is not None:
            last["ids"] = ""
            tables.append(last)
        return tables

    # With markers: for each run (marker i .. next marker), keep the LAST
    # 'Refinement steps' table found within that span.
    for idx, (mstart, ids) in enumerate(markers):
        mend = markers[idx + 1][0] if idx + 1 < len(markers) else n
        last = None
        i = mstart + 1
        while i < mend:
            if _REFINE_HEADER_RE.search(lines[i]):
                tbl, j = _parse_one_refine_table(lines, i + 1)
                if tbl["step"]:
                    last = tbl
                i = max(j, i + 1)
            else:
                i += 1
        if last is not None:
            last["ids"] = ids
            tables.append(last)
    return tables


_FRAMES_RE = re.compile(r"Frames:\s*(\d+)\s*->\s*(\d+)")


def parse_integrate_blocks(text: str) -> List[Tuple[int, int]]:
    """integrate: the block table near the start -> list of
    (frame_from, frame_to). [] until the table appears. Identified by a
    header row containing 'Frame From' / 'Frame To'."""
    lines = text.splitlines()
    blocks: List[Tuple[int, int]] = []
    in_table = False
    header_cols: Optional[List[str]] = None

    for line in lines:
        cells = _split_table_row(line)
        if cells is None:
            if in_table and blocks and not line.strip():
                break
            continue
        lowered = [c.lower() for c in cells]
        if "frame from" in lowered and "frame to" in lowered:
            header_cols = lowered
            in_table = True
            blocks = []
            continue
        if not in_table or header_cols is None:
            continue
        if len(cells) != len(header_cols):
            continue
        try:
            ff_idx = header_cols.index("frame from")
            ft_idx = header_cols.index("frame to")
            blocks.append((int(cells[ff_idx]), int(cells[ft_idx])))
        except (ValueError, IndexError):
            continue

    return blocks


def parse_integrate_progress(text: str) -> Dict[str, int]:
    """integrate: count 'Frames: A -> B' block completions so far and
    report the most recent range. The block loop runs twice (profile
    modelling, then integration), so 'completed' can exceed the number of
    blocks; the plotting code accounts for this."""
    matches = _FRAMES_RE.findall(text)
    if not matches:
        return {"completed": 0, "last_from": 0, "last_to": 0}
    last_from, last_to = matches[-1]
    return {
        "completed": len(matches),
        "last_from": int(last_from),
        "last_to": int(last_to),
    }


# Multi-crystal indexing (joint=false) processes one imageset at a time and
# prints a line like "Indexing imageset id 19 (20/36)". The (k/N) is a
# reliable progress measure: k of N imagesets started. We take the last
# such line seen.
_INDEX_PROGRESS_RE = re.compile(
    r"Indexing\s+imageset\s+id\s+(\d+)\s*\(\s*(\d+)\s*/\s*(\d+)\s*\)",
    re.IGNORECASE,
)


def parse_index_progress(text: str) -> Optional[Dict[str, int]]:
    """index (multi): parse 'Indexing imageset id <id> (k/N)' progress
    lines. Returns {'imageset_id': id, 'done': k, 'total': N} for the most
    recent line, or None if no such line has appeared (single-crystal
    indexing, or not started yet)."""
    matches = _INDEX_PROGRESS_RE.findall(text)
    if not matches:
        return None
    iset, done, total = matches[-1]
    return {"imageset_id": int(iset), "done": int(done), "total": int(total)}


_SUMMARY_VS_IMAGE_RE = re.compile(r"Summary vs image number", re.IGNORECASE)
_INTEGRATE_SUMMARY_KEYS = [
    "image",
    "n_full",
    "n_part",
    "n_over",
    "n_ice",
    "n_sum",
    "n_prf",
    "ibg",
    "isigi_sum",
    "isigi_prf",
    "cc_prf",
    "rmsd_xy",
]


def parse_integrate_summary(text: str) -> Dict[str, List[float]]:
    """integrate (end): the 'Summary vs image number' table -> per-column
    series keyed by image number (one row per image)."""
    lines = text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if _SUMMARY_VS_IMAGE_RE.search(line):
            header_idx = i
    result: Dict[str, List[float]] = {k: [] for k in _INTEGRATE_SUMMARY_KEYS}
    if header_idx is None:
        return result

    started = False
    for line in lines[header_idx + 1 :]:
        cells = _split_table_row(line)
        if cells is None:
            if started and not line.strip():
                break
            continue
        if len(cells) < 13:
            continue
        try:
            vals = [
                float(int(cells[1])),  # image
                float(int(cells[2])),  # n_full
                float(int(cells[3])),  # n_part
                float(int(cells[4])),  # n_over
                float(int(cells[5])),  # n_ice
                float(int(cells[6])),  # n_sum
                float(int(cells[7])),  # n_prf
                float(cells[8]),  # ibg
                float(cells[9]),  # isigi_sum
                float(cells[10]),  # isigi_prf
                float(cells[11]),  # cc_prf
                float(cells[12]),  # rmsd_xy
            ]
        except (ValueError, IndexError):
            continue
        started = True
        for k, v in zip(_INTEGRATE_SUMMARY_KEYS, vals):
            result[k].append(v)

    return result


_MERGING_HEADER_RE = re.compile(r"Merging statistics by resolution bin", re.IGNORECASE)
_MERGING_KEYS = [
    "d_max",
    "d_min",
    "inv_d2",
    "mult",
    "completeness",
    "i_mean",
    "i_over_sigma",
    "r_merge",
    "r_meas",
    "r_pim",
    "r_anom",
    "cc_half",
    "cc_anom",
]


def _merging_float(token: str) -> Optional[float]:
    """Parse a merging-stats numeric token, tolerating a trailing
    significance marker '*' (e.g. '1.000*' or '0.114*')."""
    token = token.strip().rstrip("*")
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def parse_scale_merging(text: str) -> Dict[str, object]:
    """scale: the 'Merging statistics by resolution bin' block -> per-bin
    series. Because resolution bins are non-linear, an 'inv_d2' series is
    provided (1 / d^2, d = geometric mean of the bin's d_min and d_max)
    for use as the X axis. The trailing overall-summary row (spanning the
    whole resolution range) is separated out under key 'overall' so it
    doesn't distort the per-bin plots.

    This is a whitespace-columnar table, NOT pipe-delimited."""
    lines = text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if _MERGING_HEADER_RE.search(line):
            header_idx = i
    result: Dict[str, object] = {k: [] for k in _MERGING_KEYS}
    if header_idx is None:
        return result

    rows: List[List[float]] = []
    started = False
    for line in lines[header_idx + 1 :]:
        stripped = line.strip()
        if not stripped:
            if started:
                break
            continue
        if "d_max" in stripped or "d_min" in stripped:
            continue
        tokens = stripped.split()
        if len(tokens) != 14:
            if started:
                break
            continue
        vals = [_merging_float(t) for t in tokens]
        if any(v is None for v in vals):
            if started:
                break
            continue
        started = True
        rows.append(vals)  # type: ignore[arg-type]

    if not rows:
        return result

    overall_row = None
    per_bin_rows = rows
    if len(rows) >= 2:
        last = rows[-1]
        others = rows[:-1]
        min_dmin = min(r[1] for r in others)
        widest_bin_span = max(r[0] - r[1] for r in others)
        span_last = last[0] - last[1]
        if last[1] <= min_dmin + 1e-6 and span_last > widest_bin_span + 1e-6:
            overall_row = last
            per_bin_rows = others

    for r in per_bin_rows:
        (
            d_max,
            d_min,
            _nobs,
            _nuniq,
            mult,
            comp,
            i_mean,
            i_sig,
            r_mrg,
            r_meas,
            r_pim,
            r_anom,
            cc12,
            ccano,
        ) = r
        d_geo = math.sqrt(d_min * d_max)
        result["d_max"].append(d_max)  # type: ignore[union-attr]
        result["d_min"].append(d_min)  # type: ignore[union-attr]
        result["inv_d2"].append(1.0 / (d_geo * d_geo))  # type: ignore[union-attr]
        result["mult"].append(mult)  # type: ignore[union-attr]
        result["completeness"].append(comp)  # type: ignore[union-attr]
        result["i_mean"].append(i_mean)  # type: ignore[union-attr]
        result["i_over_sigma"].append(i_sig)  # type: ignore[union-attr]
        result["r_merge"].append(r_mrg)  # type: ignore[union-attr]
        result["r_meas"].append(r_meas)  # type: ignore[union-attr]
        result["r_pim"].append(r_pim)  # type: ignore[union-attr]
        result["r_anom"].append(r_anom)  # type: ignore[union-attr]
        result["cc_half"].append(cc12)  # type: ignore[union-attr]
        result["cc_anom"].append(ccano)  # type: ignore[union-attr]

    if overall_row is not None:
        result["overall"] = {
            "d_max": overall_row[0],
            "d_min": overall_row[1],
            "completeness": overall_row[5],
            "i_over_sigma": overall_row[7],
            "cc_half": overall_row[12],
        }

    return result


# --------------------------------------------------------------------------
# Multi-crystal (COWS_PIGS_PEOPLE) live-plot data extraction
# --------------------------------------------------------------------------
#
# When many data sets are imported and indexed with joint=False, find_spots,
# refine and integrate each produce output covering N data sets. These
# helpers split that output per data set so the Plots tab can show one page
# per data set. They reuse the single-data-set parsers above where the
# per-image / per-step data isn't itself tagged by data set.


_IMAGESET_RE = re.compile(r"imageset\s+(\d+)", re.IGNORECASE)
_SWEEP_COUNT_RE = re.compile(r"sweep:\s*(\d+)", re.IGNORECASE)


def parse_find_spots_histograms(text: str) -> Dict[int, int]:
    """find_spots (multi): the per-imageset histogram headers
    'NNNN spots found on 100 images' tagged 'for imageset K' -> {imageset:
    total spots}. Used to detect how many data sets there are and to show
    a per-data-set spot-total summary. Returns {} if no such headers."""
    result: Dict[int, int] = {}
    lines = text.splitlines()
    current = None
    for line in lines:
        m = _IMAGESET_RE.search(line)
        if m and "histogram" in line.lower():
            current = int(m.group(1))
            continue
        if current is not None:
            m2 = re.match(r"\s*(\d+)\s+spots found", line)
            if m2:
                result[current] = int(m2.group(1))
                current = None
    return result


_RMSD_BY_EXP_RE = re.compile(r"RMSDs?\s+by\s+experiment", re.IGNORECASE)


def parse_refine_by_experiment(text: str) -> Dict[str, List[float]]:
    """refine (multi): the 'RMSDs by experiment' table ->
    {exp, nref, rmsd_x, rmsd_y, rmsd_z}. One row per experiment/data set.
    Columns are (px)/(px)/(images) for multi-crystal refine. Reads rows
    after the last such header."""
    lines = text.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if _RMSD_BY_EXP_RE.search(line):
            header_idx = i
    keys = ["exp", "nref", "rmsd_x", "rmsd_y", "rmsd_z"]
    result: Dict[str, List[float]] = {k: [] for k in keys}
    if header_idx is None:
        return result

    started = False
    for line in lines[header_idx + 1 :]:
        cells = _split_table_row(line)
        if cells is None:
            if started and not line.strip():
                break
            continue
        if len(cells) < 5:
            continue
        try:
            vals = [
                float(int(cells[0])),  # exp id
                float(int(cells[1])),  # nref
                float(cells[2]),  # rmsd_x
                float(cells[3]),  # rmsd_y
                float(cells[4]),  # rmsd_z
            ]
        except ValueError:
            continue
        started = True
        for k, v in zip(keys, vals):
            result[k].append(v)

    return result


def integrate_summary_by_dataset(text: str) -> Dict[int, Dict[str, List[float]]]:
    """integrate (multi): split the 'Summary vs image number' data by its
    first column (the imageset / data set ID, 0..N) into
    {dataset_id: {image, n_full, ...}}.

    DIALS may print the summary either as one big table or as one table per
    imageset (each with its own 'Summary vs image number' header). Either
    way we want *every* such block, keyed by the ID column - so unlike a
    single-table parser we do NOT anchor on the last header only, and we do
    NOT stop at the first blank line (which would end after the first
    block). Instead we scan the whole text and treat any 13+ column pipe
    row whose first two cells are integers as a data row, ignoring header /
    border / unit rows. This is robust to multiple blocks and to streaming
    (partial last block).

    Rows are grouped by ID; within each ID they're kept in the order seen
    (image number order as DIALS emits them)."""
    out: Dict[int, Dict[str, List[float]]] = {}
    if _SUMMARY_VS_IMAGE_RE.search(text) is None:
        return out

    for line in text.splitlines():
        cells = _split_table_row(line)
        if cells is None or len(cells) < 13:
            continue
        # A data row: first cell is the dataset ID (int), second is the
        # image number (int). Header rows ("ID","Image",...) and the
        # "(sum)"/"(prf)" unit-continuation row fail these int parses and
        # are skipped.
        try:
            ds = int(cells[0])
            row = [
                float(int(cells[1])),  # image
                float(int(cells[2])),  # n_full
                float(int(cells[3])),  # n_part
                float(int(cells[4])),  # n_over
                float(int(cells[5])),  # n_ice
                float(int(cells[6])),  # n_sum
                float(int(cells[7])),  # n_prf
                float(cells[8]),  # ibg
                float(cells[9]),  # isigi_sum
                float(cells[10]),  # isigi_prf
                float(cells[11]),  # cc_prf
                float(cells[12]),  # rmsd_xy
            ]
        except (ValueError, IndexError):
            continue
        d = out.setdefault(ds, {k: [] for k in _INTEGRATE_SUMMARY_KEYS})
        for k, v in zip(_INTEGRATE_SUMMARY_KEYS, row):
            d[k].append(v)

    return out


_CLUSTER_HEAD_RE = re.compile(r"^\s*Cluster\s+(\d+)\s*$")
_CLUSTER_DATASETS_RE = re.compile(r"Datasets:\s*([0-9,\s]+)")
_CLUSTER_COMPLETENESS_RE = re.compile(r"Completeness:\s*([0-9.]+)")
_CLUSTER_MULTIPLICITY_RE = re.compile(r"Multiplicity:\s*([0-9.]+)")


def parse_cluster_list(text: str) -> List[Dict[str, object]]:
    """dials.correlation_matrix (stdout): parse the
    'Cluster N / Number of datasets / Completeness / Multiplicity /
    Datasets:...' blocks into a list of
    {'id': N, 'datasets': [...], 'completeness': f, 'multiplicity': f}.
    Tolerant of streaming: returns whatever complete-enough blocks exist."""
    clusters: List[Dict[str, object]] = []
    lines = text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        m = _CLUSTER_HEAD_RE.match(lines[i])
        if not m:
            i += 1
            continue
        cid = int(m.group(1))
        block = {"id": cid, "datasets": [], "completeness": None, "multiplicity": None}
        j = i + 1
        while j < n and not _CLUSTER_HEAD_RE.match(lines[j]):
            comp = _CLUSTER_COMPLETENESS_RE.search(lines[j])
            if comp:
                block["completeness"] = float(comp.group(1))
            mult = _CLUSTER_MULTIPLICITY_RE.search(lines[j])
            if mult:
                block["multiplicity"] = float(mult.group(1))
            ds = _CLUSTER_DATASETS_RE.search(lines[j])
            if ds:
                nums = [int(x) for x in re.findall(r"\d+", ds.group(1))]
                block["datasets"] = nums
            j += 1
        clusters.append(block)
        i = j
    return clusters


# --------------------------------------------------------------------------
# dials.correlation_matrix.html embedded-Plotly-JSON extraction
# --------------------------------------------------------------------------
#
# The HTML written by dials.correlation_matrix embeds several
# `var graphs_<name> = { ...JSON... };` assignments that feed Plotly. We
# pull those JSON objects out by name (balanced-brace scan, so nested
# objects are handled) and parse them with the stdlib json module - no JS
# engine needed. The graphs we know how to render with matplotlib are
# listed in CORRMAT_KNOWN_GRAPHS.


def _extract_json_object(text: str, start_brace: int) -> Optional[str]:
    """Given the index of an opening '{' in text, return the substring up
    to and including its matching '}', respecting braces inside strings."""
    depth = 0
    in_str = False
    escape = False
    for i in range(start_brace, len(text)):
        c = text[i]
        if in_str:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                in_str = False
            continue
        if c == '"':
            in_str = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[start_brace : i + 1]
    return None


def extract_corrmat_graphs(html: str) -> Dict[str, dict]:
    """Return {graph_name: parsed_json_dict} for every
    `var graphs_X = {...}` assignment in a dials.correlation_matrix.html.
    Skips any blob that fails to parse."""
    out: Dict[str, dict] = {}
    for m in re.finditer(r"var\s+(graphs_[A-Za-z0-9_]+)\s*=\s*", html):
        name = m.group(1)
        brace = html.find("{", m.end())
        if brace == -1:
            continue
        blob = _extract_json_object(html, brace)
        if blob is None:
            continue
        try:
            out[name] = json.loads(blob)
        except json.JSONDecodeError:
            continue
    return out


def corrmat_matrix(blob: dict) -> Optional[Dict[str, object]]:
    """cc_cluster / cos_angle_cluster blob -> {'z','order','title'}."""
    data = blob.get("data", [])
    heatmap = next((t for t in data if t.get("type") == "heatmap"), None)
    if heatmap is None or "z" not in heatmap:
        return None
    layout = blob.get("layout", {})
    order = (layout.get("xaxis", {}) or {}).get("ticktext")
    title = (heatmap.get("colorbar", {}) or {}).get("title", "correlation")
    return {"z": heatmap["z"], "order": order, "title": title}


def corrmat_cluster_series(blob: dict) -> List[Dict[str, object]]:
    """reachability / cosym-coordinates blob -> list of per-cluster series
    [{'name','x','y','color'}]. Infinity y-values become None."""
    out = []
    for t in blob.get("data", []):
        ys = []
        for v in t.get("y", []):
            if isinstance(v, float) and (v == float("inf") or v != v):
                ys.append(None)
            else:
                ys.append(v)
        out.append(
            {
                "name": t.get("name", ""),
                "x": t.get("x", []),
                "y": ys,
                "color": (t.get("marker", {}) or {}).get("color"),
            }
        )
    return out


def corrmat_xy(blob: dict) -> Optional[Dict[str, object]]:
    """Single-trace line/bar blob (dimensions, rij histogram) ->
    {'x','y','type','title','xtitle','ytitle'}."""
    data = blob.get("data", [])
    if not data:
        return None
    tr = data[0]
    layout = blob.get("layout", {})
    return {
        "x": tr.get("x", []),
        "y": tr.get("y", []),
        "type": tr.get("type", "line"),
        "title": layout.get("title", ""),
        "xtitle": (layout.get("xaxis", {}) or {}).get("title", ""),
        "ytitle": (layout.get("yaxis", {}) or {}).get("title", ""),
    }


# --------------------------------------------------------------------------
# Small value-holder helpers
# --------------------------------------------------------------------------
#
# The pure command-building / plot-source helpers below read GUI field
# values through a uniform `.get()` interface (a hold-over from the Tk
# StringVar/BooleanVar API). Rather than thread wx widget references through
# all of that logic, each editable field is wrapped in one of these tiny
# adapters exposing `.get()` (and `.set()` where the code writes back). The
# adapter is backed by a live wx control, so `.get()` always reflects the
# current on-screen value and `.set()` updates the control.


class _WidgetVar:
    """Adapter exposing `.get()` / `.set()` over a wx control that has
    GetValue/SetValue (TextCtrl, ComboBox, CheckBox)."""

    def __init__(self, ctrl, cast=lambda v: v):
        self.ctrl = ctrl
        self._cast = cast

    def get(self):
        return self._cast(self.ctrl.GetValue())

    def set(self, value):
        self.ctrl.SetValue(value)


class _FalseVar:
    """Stand-in for a variable that always reports False / empty - used as a
    safe default when a named field may not exist yet."""

    def get(self):
        return False


class ProcessRunner:
    """Runs a command in a background thread, streaming stdout lines to a
    queue so the wx main loop can poll it without blocking."""

    def __init__(self, cmd: List[str], cwd: str):
        self.cmd = cmd
        self.cwd = cwd
        self.q: "queue.Queue[tuple]" = queue.Queue()
        self.proc: Optional[subprocess.Popen] = None
        self._thread: Optional[threading.Thread] = None

    def start(self):
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self):
        try:
            self.proc = subprocess.Popen(
                self.cmd,
                cwd=self.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except FileNotFoundError as exc:
            self.q.put(("line", f"ERROR: could not run {self.cmd[0]}: {exc}\n"))
            self.q.put(("done", -1))
            return
        except Exception as exc:  # pragma: no cover - defensive
            self.q.put(("line", f"ERROR launching process: {exc}\n"))
            self.q.put(("done", -1))
            return

        assert self.proc.stdout is not None
        for line in self.proc.stdout:
            self.q.put(("line", line))
        rc = self.proc.wait()
        self.q.put(("done", rc))

    def terminate(self):
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
            except Exception:
                pass


# --------------------------------------------------------------------------
# Main application
# --------------------------------------------------------------------------


STATUS_ICONS = {
    "pending": "\u2b1c",  # white square
    "running": "\u23f3",  # hourglass
    "done": "\u2705",  # check mark
    "failed": "\u274c",  # cross mark
}


class DialsFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="DIALS Workflow GUI", size=(1180, 760))
        self.SetIcon(get_app_icon())

        self.workdir_value = os.getcwd()
        self.status = {s.id: "pending" for s in STEPS}
        self.selected_step: Optional[StepDef] = None
        self.field_vars: dict = {}  # step id -> {field key: _WidgetVar}
        self.input_vars: dict = {}  # step id -> [_WidgetVar per input]
        self.image_files: List[str] = []

        self.runner: Optional[ProcessRunner] = None
        self.running_step_id: Optional[str] = None

        # Live-plot state. `live_output` accumulates the raw stdout of the
        # currently running step so the plot parsers (which want the whole
        # text so far) can be re-run on each poll. The plot widgets are
        # rebuilt per select_step(); `plot_canvas` is None when the current
        # step has no Plots tab or matplotlib is unavailable.
        self.live_output: str = ""
        self.plot_canvas = None  # FigureCanvasWxAgg or None
        self.plot_figure = None  # matplotlib Figure or None
        self.plot_status_label = None  # wx.StaticText or None
        self.progress_bar = None  # wx.Gauge or None
        self.progress_label = None  # wx.StaticText or None
        self.plot_page_combo = None  # wx.ComboBox or None
        self.scale_cluster_combo = None  # wx.ComboBox or None
        self.log_cluster_combo = None  # wx.ComboBox or None
        # Redrawing the figure on every streamed line is wasteful; only
        # redraw every Nth poll or on completion.
        self._poll_tick = 0

        # Widgets rebuilt per step and referenced elsewhere.
        self.output_text = None
        self.summary_text = None
        self.log_text = None
        self.command_preview = None
        self.report_status_label = None
        self.run_button = None
        self.stop_button = None
        self.report_button = None
        self.extra_params_ctrl = None
        self.import_listbox = None
        self.notebook = None

        self._build_layout()
        self._check_dials_available()

        # Poll timer for the running subprocess. wx has no direct analogue of
        # Tk's after()-scheduled polling loop; a repeating timer polls the
        # ProcessRunner queue while a step runs and is otherwise idle.
        self._timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self._on_timer, self._timer)

        self.select_step(STEPS[0])

        # On startup, silently pick up any existing pipeline progress in the
        # initial working directory (cwd) so the GUI reflects work already
        # done there, as if it had been run through the GUI. Silent so it
        # doesn't nag when starting in an empty directory; the user can also
        # re-run this any time via the "Load state from working dir" button.
        try:
            self._load_state_from_workdir(announce=False)
        except Exception:
            pass

    # ---------------------------------------------------------- workdir --
    @property
    def workdir(self):
        # Kept as a property so the many `self.workdir.get()` call sites in
        # the ported pure-logic helpers continue to work unchanged.
        parent = self

        class _WD:
            def get(_self):
                return parent.workdir_value

            def set(_self, v):
                parent.workdir_value = v
                if parent.workdir_ctrl is not None:
                    parent.workdir_ctrl.SetValue(v)

        return _WD()

    # ---------------------------------------------------------- top bar --
    def _build_layout(self):
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)

        # --- top bar ---
        top = wx.BoxSizer(wx.HORIZONTAL)
        top.Add(
            wx.StaticText(panel, label="Working directory:"),
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.ALL,
            4,
        )
        self.workdir_ctrl = wx.TextCtrl(panel, value=self.workdir_value, size=(480, -1))
        top.Add(self.workdir_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
        self.workdir_ctrl.Bind(
            wx.EVT_TEXT,
            lambda _e: setattr(self, "workdir_value", self.workdir_ctrl.GetValue()),
        )
        browse_btn = wx.Button(panel, label="Browse...")
        browse_btn.Bind(wx.EVT_BUTTON, lambda _e: self._choose_workdir())
        top.Add(browse_btn, 0, wx.ALL, 4)
        load_btn = wx.Button(panel, label="Load state from working dir")
        load_btn.Bind(wx.EVT_BUTTON, lambda _e: self._load_state_from_workdir())
        top.Add(load_btn, 0, wx.ALL, 4)
        self.dials_status_label = wx.StaticText(panel, label="")
        self.dials_status_label.SetForegroundColour(wx.RED)
        top.Add(self.dials_status_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 8)
        outer.Add(top, 0, wx.EXPAND)

        # --- body: sidebar + main ---
        body = wx.BoxSizer(wx.HORIZONTAL)

        side = wx.BoxSizer(wx.VERTICAL)
        hdr = wx.StaticText(panel, label="Pipeline steps")
        hdr.SetFont(hdr.GetFont().Bold())
        side.Add(hdr, 0, wx.ALL, 4)

        self.step_buttons: dict = {}
        for s in STEPS:
            row = wx.BoxSizer(wx.HORIZONTAL)
            icon = wx.StaticText(panel, label=STATUS_ICONS["pending"], size=(20, -1))
            row.Add(icon, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
            btn = wx.Button(panel, label=s.title, size=(260, -1))
            btn.Bind(wx.EVT_BUTTON, lambda _e, st=s: self.select_step(st))
            row.Add(btn, 1, wx.EXPAND)
            side.Add(row, 0, wx.EXPAND | wx.BOTTOM, 1)
            self.step_buttons[s.id] = (btn, icon)

        side.Add(wx.StaticLine(panel), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 8)
        vh = wx.StaticText(panel, label="Viewing tools")
        vh.SetFont(vh.GetFont().Bold())
        side.Add(vh, 0, wx.ALL, 4)
        for label, program, arg_labels in TOOLS:
            b = wx.Button(panel, label=label)
            b.Bind(
                wx.EVT_BUTTON,
                lambda _e, p=program, a=arg_labels: self.launch_tool(p, a),
            )
            side.Add(b, 0, wx.EXPAND | wx.BOTTOM, 1)

        side.Add(wx.StaticLine(panel), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 8)
        reset_btn = wx.Button(panel, label="Reset all step statuses")
        reset_btn.Bind(wx.EVT_BUTTON, lambda _e: self._reset_statuses())
        side.Add(reset_btn, 0, wx.EXPAND)

        body.Add(side, 0, wx.EXPAND | wx.ALL, 4)

        # main panel holds the per-step notebook, rebuilt by select_step().
        self.main_panel = wx.Panel(panel)
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_panel.SetSizer(self.main_sizer)
        body.Add(self.main_panel, 1, wx.EXPAND | wx.ALL, 6)

        outer.Add(body, 1, wx.EXPAND)
        panel.SetSizer(outer)
        self._root_panel = panel

    def _check_dials_available(self):
        if shutil.which("dials.import") is None:
            self.dials_status_label.SetLabel(
                "Warning: dials.import not found on $PATH - "
                "make sure your DIALS environment is set up."
            )
            self.dials_status_label.SetForegroundColour(wx.RED)
        else:
            self.dials_status_label.SetLabel("DIALS found on $PATH")
            self.dials_status_label.SetForegroundColour(wx.Colour(0, 128, 0))
        self.dials_status_label.GetParent().Layout()

    def _choose_workdir(self):
        dlg = wx.DirDialog(
            self, "Choose working directory", defaultPath=self.workdir_value
        )
        if dlg.ShowModal() == wx.ID_OK:
            d = dlg.GetPath()
            self.workdir_value = d
            self.workdir_ctrl.SetValue(d)
            # Reflect any existing progress in the newly-chosen directory.
            self.image_files = []
            try:
                self._load_state_from_workdir(announce=False)
            except Exception:
                pass
        dlg.Destroy()

    def _reset_statuses(self):
        for s in STEPS:
            self.status[s.id] = "pending"
            _, icon = self.step_buttons[s.id]
            icon.SetLabel(STATUS_ICONS["pending"])

    def _step_outputs_present(self, step: StepDef) -> bool:
        """True if this step looks 'done' judging by files in the working
        directory: its declared output files exist, or (for steps that
        declare none) its log file exists."""
        workdir = self.workdir.get()

        def here(name: str) -> bool:
            return os.path.exists(os.path.join(workdir, name))

        # merge/export writes an MTZ or a log depending on the mode; treat
        # either produced log as done.
        if step.id == "merge_export":
            return (
                here("dials.merge.log")
                or here("dials.export.log")
                or here("merged.mtz")
                or here("scaled.mtz")
            )
        # scale may have run per-cluster (no plain scaled.expt) - accept any
        # scale log/result as evidence it ran.
        if step.id == "scale":
            if here("scaled.expt") or here("dials.scale.log"):
                return True
            try:
                for nm in os.listdir(workdir):
                    if re.match(r"dials\.scale\.cluster_\d+\.log$", nm) or re.match(
                        r"scaled_cluster_\d+\.expt$", nm
                    ):
                        return True
            except OSError:
                pass
            return False
        # correlation_matrix declares no outputs; use its HTML/log.
        if step.id == "correlation_matrix":
            return (
                here("dials.correlation_matrix.html")
                or here("dials.correlation_matrix.log")
                or here("dials.correlation_matrix.scaled.html")
            )

        if step.outputs:
            return all(here(o) for o in step.outputs if o.endswith((".expt", ".refl")))
        if step.log_file:
            return here(step.log_file)
        return False

    def _load_state_from_workdir(self, announce: bool = True):
        """Infer pipeline progress from files already in the working
        directory and update the step status icons accordingly, as if the
        steps had been run through the GUI. Also repopulates the Import file
        list from imported.expt if present. Non-destructive: it only marks
        steps done where evidence exists; others are left pending."""
        workdir = self.workdir.get()
        if not os.path.isdir(workdir):
            if announce:
                wx.MessageBox(
                    f"Working directory does not exist:\n{workdir}",
                    "Load state",
                    wx.OK | wx.ICON_WARNING,
                )
            return 0

        done = 0
        for s in STEPS:
            if self._step_outputs_present(s):
                self.status[s.id] = "done"
                self.step_buttons[s.id][1].SetLabel(STATUS_ICONS["done"])
                done += 1
            else:
                # don't clobber a 'running' state; otherwise reset to pending
                if self.status.get(s.id) != "running":
                    self.status[s.id] = "pending"
                    self.step_buttons[s.id][1].SetLabel(STATUS_ICONS["pending"])

        # If Import ran, reflect imported.expt as the import 'file' so the
        # Import command preview and downstream defaults make sense. We only
        # set this if the user hasn't already queued specific images.
        if (
            os.path.exists(os.path.join(workdir, "imported.expt"))
            and not self.image_files
        ):
            self.image_files = ["imported.expt"]

        # Refresh the currently-displayed step so its Log/Plots/inputs pick
        # up whatever is now on disk.
        if self.selected_step is not None:
            self.select_step(self.selected_step)

        if announce:
            wx.MessageBox(
                f"Marked {done} step(s) as done based on files in\n{workdir}",
                "Load state",
                wx.OK | wx.ICON_INFORMATION,
            )
        return done

    # ---------------------------------------------------- step display --
    def select_step(self, step: StepDef):
        self.selected_step = step

        # Rebuild the main panel's notebook from scratch for this step.
        self.main_sizer.Clear(delete_windows=True)
        self.notebook = wx.Notebook(self.main_panel)

        setup_tab = wx.Panel(self.notebook)
        output_tab = wx.Panel(self.notebook)
        summary_tab = wx.Panel(self.notebook)
        log_tab = wx.Panel(self.notebook)
        self.notebook.AddPage(setup_tab, "Setup & Run")
        self.notebook.AddPage(output_tab, "Live Output")
        self.notebook.AddPage(summary_tab, "Summary")
        self.notebook.AddPage(log_tab, "Full Log")

        self._build_setup_tab(setup_tab, step)
        self.output_text = self._make_readonly_text(output_tab)
        self.summary_text = self._make_readonly_text(summary_tab)
        self.log_text = self._make_readonly_text(log_tab)

        # Reset per-step plot state first so a step without plots doesn't
        # inherit a stale canvas.
        self.plot_canvas = None
        self.plot_figure = None
        self.plot_status_label = None
        self.progress_bar = None
        self.progress_label = None
        self.plot_page_combo = None
        if step.plot_kind and HAVE_MPL:
            plots_tab = wx.Panel(self.notebook)
            self.notebook.AddPage(plots_tab, "Plots")
            self._build_plots_tab(plots_tab, step)
        elif step.plot_kind and not HAVE_MPL:
            plots_tab = wx.Panel(self.notebook)
            self.notebook.AddPage(plots_tab, "Plots")
            s = wx.BoxSizer(wx.VERTICAL)
            lbl = wx.StaticText(
                plots_tab,
                label=(
                    "Live plots for this step need matplotlib, which isn't "
                    "installed in this environment.\n\nInstall it into your "
                    "DIALS/Python environment (e.g. `pip install matplotlib` "
                    "or `libtbx.pip install matplotlib`) and restart the GUI "
                    "to enable the Plots tab. Everything else — running the "
                    "step, the log Summary, and the 'Run and show report in "
                    "web browser' button — works without it."
                ),
            )
            lbl.Wrap(760)
            s.Add(lbl, 0, wx.ALL, 8)
            plots_tab.SetSizer(s)

        # Full Log tab controls. For the scale step, add a selector to choose
        # WHICH log to show: the plain dials.scale.log (default) or any
        # completed cluster's dials.scale.cluster_N.log. This is independent
        # of the Plots-tab cluster view.
        self.log_cluster_combo = None
        log_ctrl = wx.BoxSizer(wx.HORIZONTAL)
        refresh_log_btn = wx.Button(log_tab, label="Refresh from log file")
        refresh_log_btn.Bind(
            wx.EVT_BUTTON, lambda _e, st=step: self._refresh_log_tab(st)
        )
        log_ctrl.Add(refresh_log_btn, 0, wx.ALL, 2)
        if step.id == "scale":
            result_clusters = self._scale_result_clusters()
            choices = ["dials.scale.log (default)"] + [
                f"cluster_{c}" for c in result_clusters
            ]
            log_ctrl.Add(
                wx.StaticText(log_tab, label="   Log:"),
                0,
                wx.ALIGN_CENTER_VERTICAL | wx.LEFT,
                6,
            )
            self.log_cluster_combo = wx.ComboBox(
                log_tab,
                choices=choices,
                value=choices[0],
                style=wx.CB_READONLY,
                size=(200, -1),
            )
            log_ctrl.Add(self.log_cluster_combo, 0, wx.ALL, 2)
            self.log_cluster_combo.Bind(
                wx.EVT_COMBOBOX,
                lambda _e, st=step: self._refresh_log_tab(st),
            )
        # Insert the log controls above the read-only text in the log tab.
        log_sizer = log_tab.GetSizer()
        log_sizer.Insert(0, log_ctrl, 0, wx.EXPAND)
        log_tab.Layout()

        self._refresh_log_tab(step)
        # If a log already exists for this step (e.g. re-selecting a step that
        # ran earlier), populate the plots from it immediately.
        if step.plot_kind and HAVE_MPL:
            self._refresh_plots_from_text(step, self._plot_source_text(step))

        self.main_sizer.Add(self.notebook, 1, wx.EXPAND)
        self.main_panel.Layout()

    def _make_readonly_text(self, parent) -> wx.TextCtrl:
        sizer = parent.GetSizer()
        if sizer is None:
            sizer = wx.BoxSizer(wx.VERTICAL)
            parent.SetSizer(sizer)
        txt = wx.TextCtrl(
            parent,
            style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP | wx.HSCROLL,
        )
        txt.SetFont(wx.Font(wx.FontInfo(10).Family(wx.FONTFAMILY_TELETYPE)))
        sizer.Add(txt, 1, wx.EXPAND | wx.ALL, 4)
        return txt

    def _set_text(self, widget: wx.TextCtrl, content: str):
        widget.ChangeValue(content)

    def _append_text(self, widget: wx.TextCtrl, content: str):
        widget.AppendText(content)

    def _build_setup_tab(self, parent, step: StepDef):
        sizer = wx.BoxSizer(wx.VERTICAL)
        parent.SetSizer(sizer)

        help_lbl = wx.StaticText(parent, label=step.help)
        help_lbl.Wrap(760)
        sizer.Add(help_lbl, 0, wx.ALL, 6)

        self.input_vars[step.id] = []
        if step.is_import:
            self._build_import_inputs(parent, sizer, step)
        else:
            for spec in step.inputs:
                row = wx.BoxSizer(wx.HORIZONTAL)
                row.Add(
                    wx.StaticText(parent, label=spec.label, size=(170, -1)),
                    0,
                    wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                    2,
                )
                ctrl = wx.TextCtrl(parent, value=spec.default, size=(360, -1))
                row.Add(ctrl, 1, wx.ALL, 2)
                sizer.Add(row, 0, wx.EXPAND)
                var = _WidgetVar(ctrl)
                self.input_vars[step.id].append(var)
                ctrl.Bind(wx.EVT_TEXT, lambda _e: self._update_command_preview())

        self.field_vars[step.id] = {}

        # Scale step: cluster selector for multi-crystal cluster scaling.
        # If cluster_N.expt/.refl files exist (written by the Correlation
        # Matrix step with 'output clusters' ticked), let the user pick one to
        # scale independently; picking a cluster rewrites the input files and
        # adds distinct output.* names so runs don't overwrite.
        self.scale_cluster_combo = None
        if step.id == "scale":
            clusters = self._available_clusters()
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(
                wx.StaticText(parent, label="Cluster to scale", size=(170, -1)),
                0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                2,
            )
            choices = ["(none - use inputs above)"] + [f"cluster_{c}" for c in clusters]
            self.scale_cluster_combo = wx.ComboBox(
                parent,
                choices=choices,
                value=choices[0],
                style=wx.CB_READONLY,
                size=(220, -1),
            )
            row.Add(self.scale_cluster_combo, 0, wx.ALL, 2)
            if clusters:
                note = (
                    f"{len(clusters)} cluster(s) found: "
                    f"{', '.join(str(c) for c in clusters)}"
                )
            else:
                note = (
                    "(no cluster_N files yet - run Correlation Matrix "
                    "with 'output clusters')"
                )
            note_lbl = wx.StaticText(parent, label=note)
            note_lbl.SetForegroundColour(wx.Colour(128, 128, 128))
            row.Add(note_lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 8)
            sizer.Add(row, 0, wx.EXPAND | wx.TOP, 4)

            def _on_cluster_change(_e):
                c = self._selected_cluster()
                exp_var, refl_var = self.input_vars["scale"][:2]
                if c is not None:
                    exp_var.set(f"cluster_{c}.expt")
                    refl_var.set(f"cluster_{c}.refl")
                else:
                    exp_var.set("symmetrized.expt")
                    refl_var.set("symmetrized.refl")
                self._update_command_preview()

            self.scale_cluster_combo.Bind(wx.EVT_COMBOBOX, _on_cluster_change)

        if step.extra_fields:
            ph = wx.StaticText(parent, label="Parameters:")
            ph.SetFont(ph.GetFont().Bold())
            sizer.Add(ph, 0, wx.ALL, 4)
        for f in step.extra_fields:
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(
                wx.StaticText(parent, label=f.label, size=(170, -1)),
                0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                2,
            )
            if f.kind == "check":
                checked = str(f.default).strip().lower() in ("true", "1", "yes")
                cb = wx.CheckBox(parent)
                cb.SetValue(checked)
                row.Add(cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
                var = _WidgetVar(cb)
                cb.Bind(wx.EVT_CHECKBOX, lambda _e: self._update_command_preview())
            elif f.kind == "combo":
                cb = wx.ComboBox(
                    parent, value=f.default, choices=f.choices or [], size=(160, -1)
                )
                row.Add(cb, 0, wx.ALL, 2)
                var = _WidgetVar(cb)
                cb.Bind(wx.EVT_COMBOBOX, lambda _e: self._update_command_preview())
                cb.Bind(wx.EVT_TEXT, lambda _e: self._update_command_preview())
            else:
                ctrl = wx.TextCtrl(parent, value=f.default, size=(220, -1))
                row.Add(ctrl, 0, wx.ALL, 2)
                var = _WidgetVar(ctrl)
                ctrl.Bind(wx.EVT_TEXT, lambda _e: self._update_command_preview())
            if f.help:
                hl = wx.StaticText(parent, label=f.help)
                hl.SetForegroundColour(wx.Colour(128, 128, 128))
                row.Add(hl, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 8)
            sizer.Add(row, 0, wx.EXPAND)
            self.field_vars[step.id][f.key] = var

        sizer.Add(
            wx.StaticText(parent, label="Additional parameters (free text):"),
            0,
            wx.LEFT | wx.TOP,
            6,
        )
        self.extra_params_ctrl = wx.TextCtrl(parent, value="", size=(560, -1))
        self.extra_params_ctrl.Bind(
            wx.EVT_TEXT, lambda _e: self._update_command_preview()
        )
        sizer.Add(self.extra_params_ctrl, 0, wx.LEFT | wx.BOTTOM, 6)

        cph = wx.StaticText(parent, label="Command preview:")
        cph.SetFont(cph.GetFont().Bold())
        sizer.Add(cph, 0, wx.LEFT | wx.TOP, 6)
        self.command_preview = wx.StaticText(parent, label="")
        self.command_preview.SetForegroundColour(wx.BLUE)
        sizer.Add(self.command_preview, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)

        # For correlation_matrix, toggling 'use scaled data' changes both the
        # input files and which HTML/log files the Plots/Log tabs read.
        if step.id == "correlation_matrix":
            us = self.field_vars[step.id].get("use_scaled")
            if us is not None:

                def _on_use_scaled(_e):
                    exp_var, refl_var = self.input_vars["correlation_matrix"][:2]
                    if bool(us.get()):
                        exp_var.set("scaled.expt")
                        refl_var.set("scaled.refl")
                    else:
                        exp_var.set("symmetrized.expt")
                        refl_var.set("symmetrized.refl")
                    self._update_command_preview()
                    self._refresh_log_tab(step)
                    if step.plot_kind and HAVE_MPL:
                        self._refresh_plots_from_text(
                            step, self._plot_source_text(step)
                        )

                us.ctrl.Bind(wx.EVT_CHECKBOX, _on_use_scaled)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.run_button = wx.Button(
            parent,
            label=f"Run {step.program}" + (" (optional)" if step.optional else ""),
        )
        self.run_button.Bind(wx.EVT_BUTTON, lambda _e: self.run_step(step))
        btn_row.Add(self.run_button, 0, wx.ALL, 2)
        self.stop_button = wx.Button(parent, label="Stop")
        self.stop_button.Enable(False)
        self.stop_button.Bind(wx.EVT_BUTTON, lambda _e: self.stop_running())
        btn_row.Add(self.stop_button, 0, wx.ALL, 2)
        self.report_button = wx.Button(
            parent, label="Run and show report in web browser"
        )
        self.report_button.Bind(
            wx.EVT_BUTTON, lambda _e: self.run_and_show_report(step)
        )
        btn_row.Add(self.report_button, 0, wx.ALL, 2)
        sizer.Add(btn_row, 0, wx.TOP, 8)

        self.report_status_label = wx.StaticText(parent, label="")
        self.report_status_label.SetForegroundColour(wx.Colour(128, 128, 128))
        sizer.Add(self.report_status_label, 0, wx.ALL, 4)

        self._update_command_preview()

    def _build_import_inputs(self, parent, sizer, step: StepDef):
        sizer.Add(
            wx.StaticText(parent, label="Image files / master file(s):"),
            0,
            wx.LEFT | wx.TOP,
            4,
        )
        self.import_listbox = wx.ListBox(parent, size=(-1, 100), style=wx.LB_SINGLE)
        for f in self.image_files:
            self.import_listbox.Append(f)
        sizer.Add(self.import_listbox, 0, wx.EXPAND | wx.ALL, 4)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        browse = wx.Button(parent, label="Browse files...")
        browse.Bind(wx.EVT_BUTTON, lambda _e: self._browse_images())
        btn_row.Add(browse, 0, wx.ALL, 2)
        addglob = wx.Button(parent, label="Add glob pattern...")
        addglob.Bind(wx.EVT_BUTTON, lambda _e: self._add_glob_pattern())
        btn_row.Add(addglob, 0, wx.ALL, 2)
        clear = wx.Button(parent, label="Clear")
        clear.Bind(wx.EVT_BUTTON, lambda _e: self._clear_images())
        btn_row.Add(clear, 0, wx.ALL, 2)
        sizer.Add(btn_row, 0)

    def _browse_images(self):
        dlg = wx.FileDialog(
            self,
            "Select image / master files",
            defaultDir=self.workdir.get(),
            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_FILE_MUST_EXIST,
        )
        if dlg.ShowModal() == wx.ID_OK:
            for f in dlg.GetPaths():
                self.image_files.append(f)
                if self.import_listbox is not None:
                    self.import_listbox.Append(f)
            self._update_command_preview()
        dlg.Destroy()

    def _add_glob_pattern(self):
        dlg = wx.TextEntryDialog(
            self,
            "Enter a glob pattern (e.g. ../data/CIX*gz or ../data/ins10_?.nxs).\n"
            "The pattern is passed to dials.import as-is (not expanded here), "
            "so it's fine for it to match thousands of images:",
            "Glob pattern",
        )
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        pattern = dlg.GetValue().strip()
        dlg.Destroy()
        if not pattern:
            return
        # Pass the pattern through verbatim - dials.import does its own shell-
        # style expansion, and for large sweeps expanding here would put
        # thousands of paths on the command line (and in the listbox). Just do
        # a quick, non-authoritative count as a sanity hint to the user.
        try:
            n = len(glob.glob(pattern))
        except Exception:
            n = None
        self.image_files.append(pattern)
        hint = "" if n is None else f"  [matches {n} file(s) now]"
        if self.import_listbox is not None:
            self.import_listbox.Append(pattern + hint)
        self._update_command_preview()

    def _clear_images(self):
        self.image_files = []
        if self.import_listbox is not None:
            self.import_listbox.Clear()
        self._update_command_preview()

    # ---------------------------------------------------- command build --
    def _build_command(self, step: StepDef) -> List[str]:
        args: List[str] = []

        cluster = self._selected_cluster()
        # correlation_matrix "use scaled data" pseudo-toggle
        cm_use_scaled = step.id == "correlation_matrix" and bool(
            self.field_vars.get(step.id, {}).get("use_scaled", _FalseVar()).get()
        )

        if step.is_import:
            args.extend(self.image_files)
        else:
            # The input fields are the single source of truth. For cluster
            # scaling the "Cluster to scale" selector has filled them with
            # cluster_N.expt/.refl, and for correlation_matrix the "use scaled
            # data" toggle has set them to scaled.expt/.refl - so we just read
            # the fields here.
            for var in self.input_vars[step.id]:
                v = var.get().strip()
                if v:
                    args.append(v)

        for f in step.extra_fields:
            var = self.field_vars[step.id][f.key]
            val = var.get()
            arg = f.build_arg(val)
            if arg:
                args.append(arg)

        if self.extra_params_ctrl is not None:
            extra_text = self.extra_params_ctrl.GetValue().strip()
            if extra_text:
                args.extend(extra_text.split())

        program = step.program
        if step.dynamic and step.id == "merge_export":
            mode = self.field_vars[step.id].get("mode")
            mode_val = mode.get() if mode else "merge"
            program = "dials.export" if mode_val == "export" else "dials.merge"
            # remove the 'mode' pseudo-parameter, it isn't a real dials param
            args = [a for a in args if not a.startswith("mode=")]

        # Cluster scaling: redirect all outputs to cluster-tagged names so
        # repeated Scale runs (one per cluster) don't overwrite each other.
        # This mirrors the tutorial's "mkdir 0 1 2; scale in each" but keeps
        # everything in one working directory.
        if step.id == "scale" and cluster is not None:
            args.extend(
                [
                    f"output.experiments=scaled_cluster_{cluster}.expt",
                    f"output.reflections=scaled_cluster_{cluster}.refl",
                    f"output.html=dials.scale.cluster_{cluster}.html",
                    f"output.log=dials.scale.cluster_{cluster}.log",
                ]
            )

        # correlation_matrix on scaled data: drop the GUI-only 'use_scaled'
        # pseudo-flag and redirect the HTML/log to '.scaled.' names so this
        # run doesn't overwrite the earlier (post-cosym) correlation matrix.
        if step.id == "correlation_matrix":
            args = [a for a in args if not a.startswith("use_scaled=")]
            if cm_use_scaled:
                args.extend(
                    [
                        "output.html=dials.correlation_matrix.scaled.html",
                        "output.log=dials.correlation_matrix.scaled.log",
                    ]
                )

        return [program] + args

    def _update_command_preview(self):
        step = self.selected_step
        if step is None or self.command_preview is None:
            return
        try:
            cmd = self._build_command(step)
        except Exception:
            return
        self.command_preview.SetLabel(" ".join(cmd))
        self.command_preview.Wrap(880)
        parent = self.command_preview.GetParent()
        if parent is not None:
            parent.Layout()

    # -------------------------------------------------------- run step --
    def run_step(self, step: StepDef):
        if self.runner is not None and self.running_step_id is not None:
            wx.MessageBox(
                "Another step is currently running - please wait or stop it.",
                "Busy",
                wx.OK | wx.ICON_WARNING,
            )
            return

        workdir = self.workdir.get()
        if not os.path.isdir(workdir):
            wx.MessageBox(
                f"{workdir} is not a directory",
                "Invalid directory",
                wx.OK | wx.ICON_ERROR,
            )
            return

        cmd = self._build_command(step)
        if step.is_import and not self.image_files:
            wx.MessageBox(
                "Please add at least one image / master file.",
                "No input files",
                wx.OK | wx.ICON_ERROR,
            )
            return

        self._set_text(self.output_text, "")
        self._append_text(self.output_text, f"$ {' '.join(cmd)}\n\n")
        self._set_text(self.summary_text, "(running...)")

        # Reset live-plot accumulation for this run and clear any stale figure
        # so plots build up fresh as output streams in.
        self.live_output = ""
        self._poll_tick = 0
        if step.plot_kind and HAVE_MPL:
            self._refresh_plots_from_text(step, "")

        self.status[step.id] = "running"
        self.step_buttons[step.id][1].SetLabel(STATUS_ICONS["running"])
        if self.run_button is not None:
            self.run_button.Enable(False)
        if self.stop_button is not None:
            self.stop_button.Enable(True)

        self.runner = ProcessRunner(cmd, workdir)
        self.running_step_id = step.id
        self._running_step = step
        self.runner.start()
        # Poll the runner's queue via a repeating timer (~10 Hz).
        self._timer.Start(100)

    def stop_running(self):
        if self.runner:
            self.runner.terminate()

    def _on_timer(self, _evt):
        step = getattr(self, "_running_step", None)
        if self.runner is None or step is None:
            self._timer.Stop()
            return
        got_line = False
        try:
            while True:
                kind, payload = self.runner.q.get_nowait()
                if kind == "line":
                    self._append_text(self.output_text, payload)
                    self.live_output += payload
                    got_line = True
                elif kind == "done":
                    # Flush a final plot update from everything streamed, then
                    # fall through to _finish_step (which also re-reads the
                    # on-disk log for the definitive version).
                    self._timer.Stop()
                    if step.plot_kind and HAVE_MPL:
                        self._refresh_plots_from_text(step, self.live_output)
                    self._finish_step(step, payload)
                    return
        except queue.Empty:
            pass

        # Live plot update, throttled: redraw roughly every ~0.5s worth of
        # polls when new output has arrived, rather than on every single line
        # (integration alone emits thousands of lines).
        if got_line and step.plot_kind and HAVE_MPL:
            self._poll_tick += 1
            if self._poll_tick % 5 == 0:
                self._refresh_plots_from_text(step, self.live_output)

    def _finish_step(self, step: StepDef, returncode: int):
        ok = returncode == 0
        self.status[step.id] = "done" if ok else "failed"
        self.step_buttons[step.id][1].SetLabel(STATUS_ICONS["done" if ok else "failed"])
        if self.run_button is not None:
            self.run_button.Enable(True)
        if self.stop_button is not None:
            self.stop_button.Enable(False)
        self.runner = None
        self.running_step_id = None
        self._running_step = None
        # If a scale cluster run just finished, point the Plots "view cluster"
        # page at it so the results are shown immediately (and the newly-
        # written log is now on disk to build the page list from). Keep the
        # Full Log tab on its current selection (default plain), but refresh
        # that selector's choices so the new cluster is now pickable there too.
        if step.id == "scale" and ok:
            run_cluster = self._selected_cluster()
            if run_cluster is not None and self.plot_page_combo is not None:
                self._select_plot_page(f"cluster_{run_cluster}")
            combo = self.log_cluster_combo
            if combo is not None:
                choices = ["dials.scale.log (default)"] + [
                    f"cluster_{c}" for c in self._scale_result_clusters()
                ]
                cur = combo.GetValue()
                combo.Set(choices)
                combo.SetValue(cur if cur in choices else choices[0])
        self._refresh_log_tab(step)
        # Definitive plot update from the on-disk log (the streamed stdout and
        # the log file should agree, but the log is canonical; the "Summary vs
        # image number" and merging-stats tables in particular are written at
        # the very end).
        if step.plot_kind and HAVE_MPL:
            src = self._plot_source_text(step)
            if step.plot_kind == "correlation_matrix":
                self._refresh_plots_from_text(step, src)
            else:
                self._refresh_plots_from_text(
                    step, src if src.strip() else self.live_output
                )
        if not ok:
            self._append_text(
                self.output_text, f"\n[process exited with code {returncode}]\n"
            )

    # ------------------------------------------------------- clusters --
    def _available_clusters(self) -> List[int]:
        """Scan the working directory for cluster_N.expt files (written by
        dials.correlation_matrix significant_clusters.output=True) and return
        the sorted list of cluster indices N."""
        workdir = self.workdir.get()
        out = []
        try:
            for name in os.listdir(workdir):
                m = re.match(r"cluster_(\d+)\.expt$", name)
                if m:
                    out.append(int(m.group(1)))
        except OSError:
            return []
        return sorted(out)

    def _selected_cluster(self) -> Optional[int]:
        """Return the cluster index currently chosen in the scale step's
        'Cluster to scale' selector (what to RUN next), or None if 'none' /
        not applicable."""
        combo = self.scale_cluster_combo
        if combo is None:
            return None
        val = combo.GetValue()
        m = re.match(r"cluster_(\d+)$", val or "")
        return int(m.group(1)) if m else None

    def _scale_result_clusters(self) -> List[int]:
        """Cluster indices that already have a scale result on disk
        (dials.scale.cluster_N.log). These are the clusters whose results can
        be VIEWED in the Plots / Full Log tabs, independent of which cluster
        is queued to run."""
        workdir = self.workdir.get()
        out = []
        try:
            for name in os.listdir(workdir):
                m = re.match(r"dials\.scale\.cluster_(\d+)\.log$", name)
                if m:
                    out.append(int(m.group(1)))
        except OSError:
            return []
        return sorted(out)

    def _scale_view_cluster(self) -> Optional[int]:
        """Which cluster's *results* the scale Plots/Log tabs should show,
        taken from the Plots-tab page selector. 'plain' or 'all' -> None (the
        non-cluster dials.scale.log). This is separate from _selected_cluster
        (the run target) so you can review cluster 0's results while cluster 1
        is queued to run."""
        page = self._current_plot_page()
        m = re.match(r"cluster_(\d+)$", page or "")
        return int(m.group(1)) if m else None

    def _scale_log_name(self) -> str:
        """The scale log filename to READ for the Plots / Full Log tabs.

        Prefers the Plots-tab 'view cluster' selection (so completed clusters
        can be reviewed while another is queued). Falls back to the run-target
        cluster (useful mid-run before the page list is built), then to the
        plain dials.scale.log."""
        view = self._scale_view_cluster()
        if view is not None:
            return f"dials.scale.cluster_{view}.log"
        run = self._selected_cluster()
        page = self._current_plot_page()
        if run is not None and page in (None, "", "all", "plain"):
            if page != "plain":
                return f"dials.scale.cluster_{run}.log"
        return "dials.scale.log"

    def _scale_log_view_name(self) -> str:
        """The scale log filename the FULL LOG tab should show, from its own
        'Log:' selector (default dials.scale.log). Independent of the Plots
        tab's cluster view."""
        combo = self.log_cluster_combo
        if combo is not None:
            m = re.match(r"cluster_(\d+)$", combo.GetValue() or "")
            if m:
                return f"dials.scale.cluster_{m.group(1)}.log"
        return "dials.scale.log"

    def _read_workdir_file(self, name: str) -> str:
        """Read a file from the working directory by name. '' if absent."""
        path = os.path.join(self.workdir.get(), name)
        if not os.path.exists(path):
            return ""
        try:
            with open(path, "r", errors="replace") as fh:
                return fh.read()
        except OSError:
            return ""

    def _corrmat_html_text(self) -> str:
        """Read the correlation-matrix HTML from the working directory (the
        source for the correlation_matrix Plots tab). If 'use scaled data' is
        ticked, read the '.scaled.' variant this GUI writes for post-scaling
        runs; otherwise the default. '' if absent."""
        use_scaled = bool(
            self.field_vars.get("correlation_matrix", {})
            .get("use_scaled", _FalseVar())
            .get()
        )
        name = (
            "dials.correlation_matrix.scaled.html"
            if use_scaled
            else "dials.correlation_matrix.html"
        )
        return self._read_workdir_file(name)

    def _plot_source_text(self, step: StepDef) -> str:
        """The text a step's Plots tab should parse: the correlation_matrix
        step plots from its HTML output (which carries the Plotly JSON blobs),
        every other plot step from its .log file."""
        if step.plot_kind == "correlation_matrix":
            return self._corrmat_html_text()
        return self._current_log_text(step)

    def _corrmat_log_name(self) -> str:
        """The log filename dials.correlation_matrix writes given the current
        'use scaled data' selection."""
        use_scaled = bool(
            self.field_vars.get("correlation_matrix", {})
            .get("use_scaled", _FalseVar())
            .get()
        )
        return (
            "dials.correlation_matrix.scaled.log"
            if use_scaled
            else "dials.correlation_matrix.log"
        )

    def _current_log_text(self, step: StepDef) -> str:
        """Read back the on-disk log for this step (respecting the dynamic
        merge/export log-name choice and cluster scaling). Returns '' if not
        present."""
        log_name = step.log_file
        if step.dynamic and step.id == "merge_export":
            mode = self.field_vars.get(step.id, {}).get("mode")
            mode_val = mode.get() if mode else "merge"
            log_name = "dials.export.log" if mode_val == "export" else "dials.merge.log"
        elif step.id == "scale":
            log_name = self._scale_log_name()
        elif step.id == "correlation_matrix":
            log_name = self._corrmat_log_name()
        if not log_name:
            return ""
        return self._read_workdir_file(log_name)

    def _refresh_log_tab(self, step: StepDef):
        log_name = step.log_file
        if step.dynamic and step.id == "merge_export":
            mode = self.field_vars.get(step.id, {}).get("mode")
            mode_val = mode.get() if mode else "merge"
            log_name = "dials.export.log" if mode_val == "export" else "dials.merge.log"
        elif step.id == "scale":
            # The Full Log tab has its own 'Log:' selector (default
            # dials.scale.log); it is independent of the Plots-tab cluster view.
            log_name = self._scale_log_view_name()
        elif step.id == "correlation_matrix":
            log_name = self._corrmat_log_name()

        text = ""
        if log_name:
            path = os.path.join(self.workdir.get(), log_name)
            if os.path.exists(path):
                try:
                    with open(path, "r", errors="replace") as fh:
                        text = fh.read()
                except OSError as exc:
                    text = f"(could not read {path}: {exc})"
            else:
                text = f"(log file not found yet: {path})"
        else:
            text = "(no fixed log filename for this step)"

        if self.log_text is not None:
            self._set_text(self.log_text, text)
        if self.summary_text is not None:
            self._set_text(self.summary_text, summarise_log(text))

    # ------------------------------------------------------------ plots --
    def _build_plots_tab(self, parent, step: StepDef):
        """Build the Plots tab for a plot-capable step: an embedded
        matplotlib canvas (with the standard navigation toolbar), a status
        line, an optional page selector (per data set for find_spots/refine/
        integrate in multi-crystal mode, per cluster for scale), and — for
        integration — a live block-processing progress bar."""
        sizer = wx.BoxSizer(wx.VERTICAL)
        parent.SetSizer(sizer)

        top = wx.BoxSizer(wx.HORIZONTAL)
        self.plot_status_label = wx.StaticText(
            parent,
            label=(
                "(no data yet - run this step, or plots will fill in live as it runs)"
            ),
        )
        self.plot_status_label.SetForegroundColour(wx.Colour(128, 128, 128))
        top.Add(self.plot_status_label, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
        refresh_label = (
            "Refresh plots from HTML"
            if step.plot_kind == "correlation_matrix"
            else "Refresh plots from log"
        )
        refresh_btn = wx.Button(parent, label=refresh_label)
        refresh_btn.Bind(
            wx.EVT_BUTTON,
            lambda _e: self._refresh_plots_from_text(
                step, self._plot_source_text(step)
            ),
        )
        top.Add(refresh_btn, 0, wx.ALL, 2)
        sizer.Add(top, 0, wx.EXPAND)

        # Page selector: for steps that can span multiple data sets or
        # clusters, a combobox to flip between one page of plots each. Rebuilt/
        # populated lazily as data arrives (see _update_plot_pages).
        self.plot_page_combo = None
        if step.plot_kind in ("find_spots", "refine", "integrate", "scale"):
            page_row = wx.BoxSizer(wx.HORIZONTAL)
            label = "Cluster" if step.plot_kind == "scale" else "Data set"
            page_row.Add(
                wx.StaticText(parent, label=f"{label}:", size=(70, -1)),
                0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                2,
            )
            self.plot_page_combo = wx.ComboBox(
                parent,
                choices=["all"],
                value="all",
                style=wx.CB_READONLY,
                size=(180, -1),
            )
            page_row.Add(self.plot_page_combo, 0, wx.ALL, 2)

            def _on_page_change(_e, s=step):
                # Use the canonical source for the current page. For scale,
                # _plot_scale reads the selected cluster's log itself (pass "").
                # For others, prefer the live stream only while THIS step is
                # actively running; otherwise read the on-disk log, so
                # reviewing a completed step still parses correctly.
                if s.plot_kind == "scale":
                    self._refresh_plots_from_text(s, "")
                    self._refresh_log_tab(s)
                else:
                    running = self.running_step_id == s.id
                    src = (
                        self.live_output
                        if running and self.live_output
                        else self._plot_source_text(s)
                    )
                    self._refresh_plots_from_text(s, src)

            self.plot_page_combo.Bind(wx.EVT_COMBOBOX, _on_page_change)

            if step.plot_kind == "scale":
                hint = wx.StaticText(
                    parent, label="(pick a completed cluster to view its results)"
                )
                hint.SetForegroundColour(wx.Colour(128, 128, 128))
                page_row.Add(hint, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 8)
            sizer.Add(page_row, 0, wx.EXPAND | wx.TOP, 2)

        # Integration and multi-crystal indexing get a live progress bar.
        if step.plot_kind in ("integrate", "index"):
            prog_row = wx.BoxSizer(wx.HORIZONTAL)
            initial = (
                "Blocks: waiting..."
                if step.plot_kind == "integrate"
                else "Indexing: waiting..."
            )
            self.progress_label = wx.StaticText(parent, label=initial, size=(280, -1))
            prog_row.Add(self.progress_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
            self.progress_bar = wx.Gauge(parent, range=100, size=(400, -1))
            prog_row.Add(self.progress_bar, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 6)
            sizer.Add(prog_row, 0, wx.EXPAND | wx.TOP, 4)

        self.plot_figure = Figure(figsize=(7.5, 5.0), dpi=100)
        self.plot_canvas = FigureCanvas(parent, -1, self.plot_figure)
        sizer.Add(self.plot_canvas, 1, wx.EXPAND | wx.ALL, 2)
        toolbar = NavigationToolbar(self.plot_canvas)
        toolbar.Realize()
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.plot_canvas.draw()

    def _update_plot_pages(self, options: List[str]):
        """Refresh the page-selector combobox's choices, preserving the
        current selection if still valid. `options` is e.g. ['all','0',
        '1',...]. No-op if there's no selector or the options are unchanged."""
        combo = self.plot_page_combo
        if combo is None:
            return
        if list(combo.GetStrings()) == options:
            return
        cur = combo.GetValue()
        combo.Set(options)
        combo.SetValue(cur if cur in options else (options[0] if options else "all"))

    def _select_plot_page(self, value: str):
        combo = self.plot_page_combo
        if combo is None:
            return
        if value in combo.GetStrings():
            combo.SetValue(value)

    def _current_plot_page(self) -> str:
        combo = self.plot_page_combo
        return combo.GetValue() if combo is not None else "all"

    def _refresh_plots_from_text(self, step: StepDef, text: str):
        """Re-parse `text` for this step and redraw the figure. Safe to call
        repeatedly (live) and with partial/empty text. Dispatches on
        step.plot_kind."""
        if not HAVE_MPL or self.plot_figure is None or self.plot_canvas is None:
            return
        if self.selected_step is None or self.selected_step.id != step.id:
            return  # user navigated away; don't draw onto another step's tab

        kind = step.plot_kind
        try:
            if kind == "find_spots":
                self._plot_find_spots(text)
            elif kind == "index":
                self._plot_index(text)
            elif kind == "refine":
                self._plot_refine(text)
            elif kind == "integrate":
                self._plot_integrate(text)
            elif kind == "scale":
                self._plot_scale(text)
            elif kind == "correlation_matrix":
                self._plot_correlation_matrix(text)
        except Exception as exc:  # pragma: no cover - defensive redraw guard
            self._set_plot_status(f"(plot error: {exc})")
            return
        self.plot_canvas.draw_idle()

    def _set_plot_status(self, msg: str):
        if self.plot_status_label is not None:
            self.plot_status_label.SetLabel(msg)
            self.plot_status_label.Wrap(700)
            parent = self.plot_status_label.GetParent()
            if parent is not None:
                parent.Layout()

    def _set_progress(self, fraction: Optional[float], label: str):
        """Update the wx.Gauge (0..100) and its label. fraction None -> 0."""
        if self.progress_bar is not None:
            val = 0 if fraction is None else int(max(0.0, min(1.0, fraction)) * 100)
            self.progress_bar.SetValue(val)
        if self.progress_label is not None:
            self.progress_label.SetLabel(label)

    def _plot_index(self, text: str):
        """Index step: for multi-crystal indexing (joint=false) DIALS prints
        'Indexing imageset id <id> (k/N)' as it works through the imagesets.
        Drive a progress bar from the (k/N) count (reliable). There's no
        per-image line graph for indexing, so the figure just carries a
        short explanatory note; the progress bar is the real content."""
        prog = parse_index_progress(text)

        if self.progress_bar is not None and self.progress_label is not None:
            if prog is not None and prog["total"] > 0:
                self._set_progress(
                    prog["done"] / prog["total"],
                    f"Indexing imageset {prog['imageset_id']} "
                    f"({prog['done']}/{prog['total']})",
                )
            else:
                self._set_progress(
                    0,
                    "Indexing: waiting... (per-imageset progress appears for "
                    "multi-crystal joint=false runs)",
                )

        fig = self.plot_figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.axis("off")
        if prog is not None and prog["total"] > 0:
            frac = prog["done"] / prog["total"]
            msg = (
                f"Indexing multiple crystals\n\n"
                f"{prog['done']} / {prog['total']} imagesets started "
                f"({frac * 100:.0f}%)\n\n"
                f"most recent: imageset id {prog['imageset_id']}"
            )
            self._set_plot_status(f"indexing {prog['done']}/{prog['total']} imagesets")
        else:
            msg = (
                "Indexing.\n\nFor multiple crystals (joint=false), a progress "
                "bar tracks the\n'Indexing imageset id <id> (k/N)' output "
                "above.\n\nFor a single crystal there is no per-imageset "
                "progress;\ncheck the Summary / Full Log tabs for the result."
            )
            self._set_plot_status("(no multi-crystal indexing progress yet)")
        ax.text(
            0.5, 0.5, msg, ha="center", va="center", fontsize=11, transform=ax.transAxes
        )
        fig.tight_layout()

    def _plot_find_spots(self, text: str):
        series = parse_find_spots_by_imageset(text)
        self._update_plot_pages(["all"])
        fig = self.plot_figure
        fig.clear()
        ax = fig.add_subplot(111)

        # Any series with actual points?
        nonempty = [s for s in series if s["image"]]
        if not nonempty:
            self._set_plot_status("(waiting for 'Found N strong pixels' lines...)")
            ax.set_title("Strong pixels per image")
            ax.set_xlabel("Image number")
            ax.set_ylabel("Strong pixels")
            fig.tight_layout()
            return

        # One line per imageset, all on the same axes so earlier imagesets
        # persist as the run proceeds. Image numbers restart at 1 for each
        # imageset, so the shared X axis is the per-imageset image number
        # and each imageset is a separate line (the legend carries the
        # imageset number from the 'Finding strong spots in imageset N'
        # banner). With more imagesets than the default colour cycle (10),
        # colours would repeat, so we draw from a larger colormap and cycle
        # marker shapes too, keeping every line visually distinct.
        try:
            import matplotlib as _mpl

            # matplotlib.colormaps (>=3.5) replaces the deprecated
            # matplotlib.cm.get_cmap; fall back for very old versions.
            try:
                cmap = _mpl.colormaps["tab20"]
            except (AttributeError, KeyError):
                import matplotlib.cm as _cm

                cmap = _cm.get_cmap("tab20")
        except Exception:
            cmap = None
        markers = [".", "o", "s", "^", "v", "D", "x", "+", "*", "<", ">", "p"]

        n = len(nonempty)
        labelled = 0
        for i, s in enumerate(nonempty):
            iset = s["imageset"]
            label = f"imageset {iset}" if iset is not None else "imageset"
            color = cmap(i % 20) if cmap is not None else None
            marker = (
                markers[(i // 20) % len(markers)]
                if n > 20
                else markers[i % len(markers)]
            )
            ax.plot(
                s["image"],
                s["pixels"],
                marker=marker,
                markersize=3,
                linewidth=1,
                color=color,
                label=label,
            )
            labelled += 1

        ax.set_title("Strong pixels found per image (one line per imageset)")
        ax.set_xlabel("Image number (within imageset)")
        ax.set_ylabel("Number of strong pixels")
        ax.grid(True, alpha=0.3)
        # Only show a legend when it stays readable; for many imagesets the
        # legend would swamp the plot, so cap it and note the count instead.
        if labelled <= 12:
            ax.legend(fontsize=7, ncol=2 if labelled > 6 else 1)

        n_sets = sum(1 for s in nonempty if s["imageset"] is not None)
        total_images = sum(len(s["image"]) for s in nonempty)
        if n_sets > 1:
            note = f"{n_sets} imagesets, {total_images} images total"
        elif n_sets == 1:
            note = f"imageset {nonempty[0]['imageset']}, {total_images} images"
        else:
            note = f"{total_images} images processed"
        self._set_plot_status(note)
        fig.tight_layout()

    def _plot_refine(self, text: str):
        fig = self.plot_figure
        fig.clear()

        # Parse every "Refinement steps" table. Single-crystal refine emits
        # one (or a few macrocycle) table(s); multi-crystal joint=false
        # refine emits one per experiment/run. We show the RMSD-vs-step
        # CONVERGENCE for each run (not just the final RMSDs), and the page
        # selector picks which run. This is the fix for "only shows final
        # RMSD per experiment": each run's full convergence is available.
        tables = parse_all_refine_steps(text)

        if not tables:
            self._update_plot_pages(["all"])
            ax = fig.add_subplot(111)
            self._set_plot_status("(waiting for the 'Refinement steps' table...)")
            ax.set_title("Refinement RMSDs vs step")
            ax.set_xlabel("Refinement step")
            fig.tight_layout()
            return

        # One page per refinement run. Label by the original experiment id
        # from the "Selected group..." marker when available ("run 1 (id
        # 0)"), else just "run k". Single-crystal refine has one run.
        def _label(k, tbl):
            ids = tbl.get("ids", "")
            return f"run {k + 1} (id {ids})" if ids != "" else f"run {k + 1}"

        page_labels = [_label(k, t) for k, t in enumerate(tables)]
        self._update_plot_pages(page_labels)
        page = self._current_plot_page()
        # Map the selected page label back to a table index.
        sel = 0
        if page in page_labels:
            sel = page_labels.index(page)
        data = tables[sel]

        steps = data["step"]
        # RMSD_X/Y are a length (mm for a single sweep, px for multi-crystal
        # refine); RMSD_Phi/Z is angular (deg) or images. Two positional on
        # the left axis, the third on a secondary right axis.
        ax = fig.add_subplot(111)
        ax.plot(steps, data["rmsd_x"], marker="o", label="RMSD_X")
        ax.plot(steps, data["rmsd_y"], marker="s", label="RMSD_Y")
        ax.set_xlabel("Refinement step")
        ax.set_ylabel("Positional RMSD (mm or px)")
        ax.grid(True, alpha=0.3)

        ax2 = ax.twinx()
        ax2.plot(
            steps, data["rmsd_phi"], marker="^", color="tab:green", label="RMSD_Phi/Z"
        )
        ax2.set_ylabel("Angular RMSD (deg) / RMSD_Z (images)")

        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8)

        if len(tables) > 1:
            ax.set_title(f"Refinement RMSDs vs step - run {sel + 1} of {len(tables)}")
            self._set_plot_status(
                f"run {sel + 1}/{len(tables)}: {len(steps)} refinement steps "
                f"(use the Data set selector to switch run)"
            )
        else:
            ax.set_title("Refinement RMSDs vs step")
            self._set_plot_status(f"{len(steps)} refinement steps")
        fig.tight_layout()

    def _plot_integrate(self, text: str):
        blocks = parse_integrate_blocks(text)
        progress = parse_integrate_progress(text)

        # --- live progress bar (block processing) ---
        if self.progress_bar is not None and self.progress_label is not None:
            n_blocks = len(blocks)
            if n_blocks:
                # The block loop runs twice (profile modelling, then
                # integration). Show progress within the current pass.
                done = progress["completed"]
                pass_no = 1 if done <= n_blocks else 2
                in_pass = done if done <= n_blocks else done - n_blocks
                in_pass = min(in_pass, n_blocks)
                label = (
                    f"Pass {pass_no}/2 - block {in_pass}/{n_blocks}"
                    if done
                    else f"Blocks: 0/{n_blocks}"
                )
                if progress["last_to"]:
                    label += (
                        f"  (frames {progress['last_from']} -> {progress['last_to']})"
                    )
                self._set_progress(in_pass / n_blocks, label)
            else:
                self._set_progress(0, "Blocks: waiting for block table...")

        # --- end-of-integration line graphs (Summary vs image number) ---
        fig = self.plot_figure
        fig.clear()

        # Split the summary table by data set (ID column). In single-crystal
        # runs there's just one data set (id 0); in multi-crystal runs there
        # are many, and the page selector picks which one to show.
        by_ds = integrate_summary_by_dataset(text)
        if by_ds:
            ds_ids = sorted(by_ds)
            # Page options: one page per data set. (No "all" overlay - with
            # many data sets that would be unreadable; flip between them.)
            self._update_plot_pages([str(d) for d in ds_ids])
            page = self._current_plot_page()
            try:
                sel = int(page)
            except ValueError:
                sel = ds_ids[0]
            if sel not in by_ds:
                sel = ds_ids[0]
            summary = by_ds[sel]
            ds_note = f"  |  data set {sel} of {len(ds_ids)}" if len(ds_ids) > 1 else ""
        else:
            summary = parse_integrate_summary(text)
            self._update_plot_pages(["all"])
            ds_note = ""

        if not summary["image"]:
            ax = fig.add_subplot(111)
            n = len(blocks)
            if n:
                self._set_plot_status(
                    f"Integrating {n} blocks - per-image summary plots will "
                    "appear when integration finishes."
                )
            else:
                self._set_plot_status("(waiting for integration to start...)")
            ax.set_title("Integration summary vs image (pending)")
            ax.set_xlabel("Image number")
            fig.tight_layout()
            return

        img = summary["image"]
        # Four stacked panels of the most useful per-image diagnostics.
        ax1 = fig.add_subplot(221)
        ax1.plot(img, summary["isigi_sum"], label="I/sigI (sum)", linewidth=1)
        ax1.plot(img, summary["isigi_prf"], label="I/sigI (prf)", linewidth=1)
        ax1.set_title("I/sigma vs image", fontsize=9)
        ax1.set_xlabel("Image", fontsize=8)
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)

        ax2 = fig.add_subplot(222)
        ax2.plot(img, summary["n_full"], label="# full", linewidth=1)
        ax2.plot(img, summary["n_part"], label="# part", linewidth=1)
        ax2.set_title("Reflection counts vs image", fontsize=9)
        ax2.set_xlabel("Image", fontsize=8)
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)

        ax3 = fig.add_subplot(223)
        ax3.plot(img, summary["cc_prf"], color="tab:purple", linewidth=1)
        ax3.set_title("CC prf vs image", fontsize=9)
        ax3.set_xlabel("Image", fontsize=8)
        ax3.grid(True, alpha=0.3)

        ax4 = fig.add_subplot(224)
        ax4.plot(img, summary["rmsd_xy"], color="tab:red", linewidth=1)
        ax4.set_title("RMSD XY vs image", fontsize=9)
        ax4.set_xlabel("Image", fontsize=8)
        ax4.grid(True, alpha=0.3)

        self._set_plot_status(f"Integration summary over {len(img)} images{ds_note}")
        fig.tight_layout()

    def _plot_scale(self, text: str):
        # Build the "view cluster" page list from clusters that have a scale
        # result on disk, so completed clusters can be reviewed even while a
        # different cluster is queued to run. Pages: "plain" (the
        # non-cluster dials.scale.log, if present) plus "cluster_N" for each
        # dials.scale.cluster_N.log found. If nothing cluster-specific
        # exists yet, fall back to a single "all" page using the text passed
        # in (covers the live-run and single-crystal cases).
        result_clusters = self._scale_result_clusters()
        plain_exists = os.path.exists(
            os.path.join(self.workdir.get(), "dials.scale.log")
        )
        pages: List[str] = []
        if plain_exists:
            pages.append("plain")
        pages.extend(f"cluster_{c}" for c in result_clusters)

        if pages:
            self._update_plot_pages(pages)
            # Read the log for the currently-selected view page (this is what
            # lets you flip between completed clusters' results).
            view = self._scale_view_cluster()
            if view is not None:
                src = self._read_workdir_file(f"dials.scale.cluster_{view}.log")
                cluster_note = f"  |  cluster {view}"
            else:
                src = self._read_workdir_file("dials.scale.log")
                cluster_note = "  |  (unclustered scale)"
            # Prefer freshly-streamed text if it clearly contains the merging
            # table and the on-disk log doesn't yet (mid-run).
            if (
                "Merging statistics by resolution bin" in (text or "")
                and "Merging statistics by resolution bin" not in src
            ):
                src = text
            text = src
        else:
            self._update_plot_pages(["all"])
            cluster = self._selected_cluster()
            cluster_note = f"  |  cluster {cluster}" if cluster is not None else ""

        data = parse_scale_merging(text)
        fig = self.plot_figure
        fig.clear()
        inv = data["inv_d2"]  # type: ignore[assignment]
        if not inv:
            ax = fig.add_subplot(111)
            self._set_plot_status(
                "(waiting for 'Merging statistics by resolution bin'..."
                + (cluster_note + ")" if cluster_note else ")")
            )
            ax.set_title("Merging statistics vs resolution (pending)")
            fig.tight_layout()
            return

        def _res_ticks(ax):
            """Label the 1/d^2 x-axis with the actual resolution (d, in A)
            at each tick so the non-linear axis stays readable."""
            import numpy as _np  # matplotlib always brings numpy

            xt = _np.linspace(min(inv), max(inv), 6)
            ax.set_xticks(xt)
            ax.set_xticklabels([f"{(1.0 / _np.sqrt(t)):.2f}" for t in xt])
            ax.set_xlabel("Resolution d (A)  [x axis linear in 1/d^2]", fontsize=8)

        # Panel 1: <I/sigI>
        ax1 = fig.add_subplot(221)
        ax1.plot(inv, data["i_over_sigma"], marker=".", linewidth=1)  # type: ignore[index]
        ax1.set_title("<I/sigma> vs resolution", fontsize=9)
        ax1.set_ylabel("<I/sigI>", fontsize=8)
        _res_ticks(ax1)
        ax1.grid(True, alpha=0.3)

        # Panel 2: CC1/2 and CC_anom
        ax2 = fig.add_subplot(222)
        ax2.plot(inv, data["cc_half"], marker=".", label="CC1/2", linewidth=1)  # type: ignore[index]
        ax2.plot(inv, data["cc_anom"], marker=".", label="CC_anom", linewidth=1)  # type: ignore[index]
        ax2.set_title("CC vs resolution", fontsize=9)
        ax2.set_ylabel("CC", fontsize=8)
        ax2.legend(fontsize=7)
        _res_ticks(ax2)
        ax2.grid(True, alpha=0.3)

        # Panel 3: R-factors
        ax3 = fig.add_subplot(223)
        ax3.plot(inv, data["r_merge"], marker=".", label="Rmerge", linewidth=1)  # type: ignore[index]
        ax3.plot(inv, data["r_meas"], marker=".", label="Rmeas", linewidth=1)  # type: ignore[index]
        ax3.plot(inv, data["r_pim"], marker=".", label="Rpim", linewidth=1)  # type: ignore[index]
        ax3.set_title("R-factors vs resolution", fontsize=9)
        ax3.set_ylabel("R", fontsize=8)
        ax3.legend(fontsize=7)
        _res_ticks(ax3)
        ax3.grid(True, alpha=0.3)

        # Panel 4: completeness & multiplicity
        ax4 = fig.add_subplot(224)
        ax4.plot(
            inv,
            data["completeness"],
            marker=".",
            color="tab:green",  # type: ignore[index]
            label="Completeness (%)",
            linewidth=1,
        )
        ax4.set_ylabel("Completeness (%)", fontsize=8)
        ax4b = ax4.twinx()
        ax4b.plot(
            inv,
            data["mult"],
            marker=".",
            color="tab:orange",  # type: ignore[index]
            label="Multiplicity",
            linewidth=1,
        )
        ax4b.set_ylabel("Multiplicity", fontsize=8)
        ax4.set_title("Completeness & multiplicity", fontsize=9)
        _res_ticks(ax4)
        ax4.grid(True, alpha=0.3)

        overall = data.get("overall")
        if isinstance(overall, dict):
            self._set_plot_status(
                f"{len(inv)} resolution bins{cluster_note}  |  overall: "
                f"d_min {overall['d_min']:.2f} A, "
                f"I/sigI {overall['i_over_sigma']:.1f}, "
                f"CC1/2 {overall['cc_half']:.3f}"
            )
        else:
            self._set_plot_status(f"{len(inv)} resolution bins{cluster_note}")
        fig.tight_layout()

    def _plot_correlation_matrix(self, text: str):
        """Plot the diagnostics embedded in dials.correlation_matrix.html.

        Unlike the other plot kinds, the source here is the HTML file
        (passed in as `text`), not a .log - it carries Plotly JSON blobs we
        parse and re-render with matplotlib. Shows: the correlation and
        cos-angle matrices (as heatmaps), the OPTICS reachability plot
        (coloured per cluster), the cosym PCA coordinates (per cluster),
        the dimensions residual curve and the Rij histogram. The page
        selector isn't used here (it's a fixed multi-panel view)."""
        import numpy as _np

        self._update_plot_pages(["all"])
        fig = self.plot_figure
        fig.clear()

        graphs = extract_corrmat_graphs(text) if text else {}
        if not graphs:
            ax = fig.add_subplot(111)
            self._set_plot_status(
                "(run dials.correlation_matrix, or use 'Refresh plots from "
                "log' - reads dials.correlation_matrix.html)"
            )
            ax.set_title("Correlation matrix analysis (pending)")
            fig.tight_layout()
            return

        panels = []  # (draw_fn, present?) collected then laid out on a grid

        cc = graphs.get("graphs_cc_cluster")
        cos = graphs.get("graphs_cos_angle_cluster")
        reach = graphs.get("graphs_reachability")
        coords = graphs.get("graphs_cosym_coordinates_principal_components")
        dims = graphs.get("graphs_dimensions")
        rij = graphs.get("graphs_cosym_rij_histogram_sg")

        def draw_matrix(ax, blob, default_title):
            m = corrmat_matrix(blob)
            if m is None:
                return
            z = _np.array(m["z"])
            im = ax.imshow(z, cmap="YlOrRd", aspect="auto")
            ax.set_title(m.get("title") or default_title, fontsize=9)
            ax.tick_params(labelsize=6)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        def draw_clusters(ax, blob, title, scatter):
            series = corrmat_cluster_series(blob)
            for s in series:
                # plotly ("rgb(0.53,0,0.59)") are 0-1 floats and matplotlib
                # wants 0-1 too, but the format differs, so skip explicit
                # colour and rely on the cycle for robustness.
                if scatter:
                    ax.scatter(s["x"], s["y"], s=10, label=s["name"])
                else:
                    ys = [_np.nan if v is None else v for v in s["y"]]
                    ax.bar(s["x"], ys, label=s["name"])
            ax.set_title(title, fontsize=9)
            ax.tick_params(labelsize=6)
            ax.legend(fontsize=6)

        def draw_xy(ax, blob, logy=False):
            xy = corrmat_xy(blob)
            if xy is None:
                return
            if xy["type"] == "bar":
                ax.bar(
                    xy["x"],
                    xy["y"],
                    width=(xy["x"][1] - xy["x"][0]) * 0.9 if len(xy["x"]) > 1 else 0.02,
                )
            else:
                ax.plot(xy["x"], xy["y"], marker=".", linewidth=1)
            if logy:
                ax.set_yscale("log")
            ax.set_title(xy["title"] or "", fontsize=9)
            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.3)

        if cc is not None:
            panels.append(lambda ax: draw_matrix(ax, cc, "Correlation matrix"))
        if cos is not None:
            panels.append(lambda ax: draw_matrix(ax, cos, "cos(angle) matrix"))
        if reach is not None:
            panels.append(
                lambda ax: draw_clusters(
                    ax, reach, "OPTICS reachability", scatter=False
                )
            )
        if coords is not None:
            panels.append(
                lambda ax: draw_clusters(
                    ax, coords, "Cosym PCA coordinates", scatter=True
                )
            )
        if dims is not None:
            panels.append(lambda ax: draw_xy(ax, dims, logy=True))
        if rij is not None:
            panels.append(lambda ax: draw_xy(ax, rij, logy=False))

        n = len(panels)
        if n == 0:
            ax = fig.add_subplot(111)
            ax.set_title("No recognised correlation-matrix graphs found", fontsize=9)
            self._set_plot_status("(no plottable graphs in the HTML)")
            fig.tight_layout()
            return

        ncols = 2
        nrows = (n + ncols - 1) // ncols
        for i, draw in enumerate(panels):
            ax = fig.add_subplot(nrows, ncols, i + 1)
            try:
                draw(ax)
            except Exception:
                ax.set_title("(failed to draw)", fontsize=8)

        # Cluster summary from the cc blob's cluster dict, if present.
        clusters = (cc or {}).get("clusters", {}) if cc else {}
        self._set_plot_status(
            f"correlation_matrix: {n} graphs from HTML"
            + (f"  |  {len(clusters)} dendrogram nodes" if clusters else "")
        )
        fig.tight_layout()

    # --------------------------------------------------------- dials.report --
    def _report_files_for_step(self, step: StepDef) -> List[str]:
        """Work out which .expt/.refl files best represent the result of this
        stage, to hand to `dials.report`. Prefers the stage's own freshly-
        written outputs; falls back to pairing a single new output with the
        other file type from the current input fields; falls back again to the
        step's current inputs for stages (like Bravais lattice determination or
        Merge/Export) that don't themselves write a fresh .expt/.refl pair."""

        expt_out = next((o for o in step.outputs if o.endswith(".expt")), None)
        refl_out = next((o for o in step.outputs if o.endswith(".refl")), None)

        input_values = [v.get().strip() for v in self.input_vars.get(step.id, [])]

        if expt_out and refl_out:
            return [expt_out, refl_out]
        if expt_out and not refl_out:
            refl_in = next((v for v in input_values if v.endswith(".refl")), "")
            return [f for f in [expt_out, refl_in] if f]
        if refl_out and not expt_out:
            expt_in = next((v for v in input_values if v.endswith(".expt")), "")
            return [f for f in [expt_in, refl_out] if f]

        if step.is_import:
            return ["imported.expt"]

        return [v for v in input_values if v]

    def run_and_show_report(self, step: StepDef):
        workdir = self.workdir.get()
        if not os.path.isdir(workdir):
            wx.MessageBox(
                f"{workdir} is not a directory",
                "Invalid directory",
                wx.OK | wx.ICON_ERROR,
            )
            return
        if shutil.which("dials.report") is None:
            wx.MessageBox(
                "dials.report was not found on $PATH.",
                "Not found",
                wx.OK | wx.ICON_ERROR,
            )
            return

        files = self._report_files_for_step(step)
        if not files:
            wx.MessageBox(
                "No experiment/reflection files are set for this step yet - "
                "fill in the input fields above (or run the step first).",
                "No files",
                wx.OK | wx.ICON_WARNING,
            )
            return

        cmd = ["dials.report"] + files
        self.report_button.Enable(False)
        self.report_status_label.SetLabel(f"Running: {' '.join(cmd)} ...")

        def worker():
            try:
                result = subprocess.run(
                    cmd, cwd=workdir, capture_output=True, text=True
                )
            except Exception as exc:
                wx.CallAfter(
                    self._report_finished,
                    step,
                    False,
                    f"Could not launch dials.report: {exc}",
                )
                return

            html_path = os.path.join(workdir, "dials.report.html")
            if result.returncode != 0 or not os.path.exists(html_path):
                tail = (
                    (result.stdout or "")[-1500:] + "\n" + (result.stderr or "")[-1500:]
                )
                wx.CallAfter(
                    self._report_finished,
                    step,
                    False,
                    f"dials.report exited with code {result.returncode}:\n{tail}",
                )
                return

            wx.CallAfter(self._report_finished, step, True, html_path)

        threading.Thread(target=worker, daemon=True).start()

    def _report_finished(self, step: StepDef, ok: bool, detail: str):
        if self.report_button is not None:
            self.report_button.Enable(True)
        if ok:
            if self.report_status_label is not None:
                self.report_status_label.SetLabel(
                    f"Opened {detail} in your web browser."
                )
            webbrowser.open(f"file://{detail}")
        else:
            if self.report_status_label is not None:
                self.report_status_label.SetLabel(
                    "dials.report failed - see error dialog."
                )
            wx.MessageBox(detail, "dials.report failed", wx.OK | wx.ICON_ERROR)

    # ------------------------------------------------------------ tools --
    def launch_tool(self, program: str, arg_labels: List[str]):
        step = self.selected_step
        defaults = []
        if step is not None:
            if step.is_import:
                defaults = ["imported.expt"]
            else:
                defaults = [v.get() for v in self.input_vars.get(step.id, [])]
        while len(defaults) < len(arg_labels):
            defaults.append("")

        dialog = wx.Dialog(self, title=f"Launch {program}", size=(560, -1))
        dsizer = wx.BoxSizer(wx.VERTICAL)
        ctrls = []
        for label, default in zip(arg_labels, defaults):
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(
                wx.StaticText(dialog, label=label, size=(220, -1)),
                0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                4,
            )
            ctrl = wx.TextCtrl(dialog, value=default, size=(300, -1))
            row.Add(ctrl, 1, wx.ALL, 4)
            dsizer.Add(row, 0, wx.EXPAND)
            ctrls.append(ctrl)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        launch_btn = wx.Button(dialog, wx.ID_OK, label="Launch")
        cancel_btn = wx.Button(dialog, wx.ID_CANCEL, label="Cancel")
        btn_row.Add(launch_btn, 0, wx.ALL, 4)
        btn_row.Add(cancel_btn, 0, wx.ALL, 4)
        dsizer.Add(btn_row, 0, wx.ALIGN_LEFT)
        dialog.SetSizerAndFit(dsizer)

        if dialog.ShowModal() == wx.ID_OK:
            args = [c.GetValue().strip() for c in ctrls if c.GetValue().strip()]
            cmd = [program] + args
            workdir = self.workdir.get()
            if shutil.which(program) is None:
                wx.MessageBox(
                    f"{program} was not found on $PATH.",
                    "Not found",
                    wx.OK | wx.ICON_ERROR,
                )
                dialog.Destroy()
                return
            try:
                subprocess.Popen(cmd, cwd=workdir)
            except Exception as exc:
                wx.MessageBox(str(exc), "Error launching tool", wx.OK | wx.ICON_ERROR)
                dialog.Destroy()
                return
            if program == "dials.report":
                # dials.report writes dials.report.html into the cwd; give it
                # a moment then try to open it in a browser.
                def _open_report():
                    html_path = os.path.join(workdir, "dials.report.html")
                    if os.path.exists(html_path):
                        webbrowser.open(f"file://{html_path}")

                wx.CallLater(4000, _open_report)
        dialog.Destroy()


def main():
    app = wx.App(False)
    frame = DialsFrame()
    frame.Show()
    app.MainLoop()


if __name__ == "__main__":
    main()
