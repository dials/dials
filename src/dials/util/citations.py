from __future__ import annotations

from dataclasses import dataclass

import dials.util.version

dials_version = dials.util.version.dials_version()


@dataclass
class SoftwareAndCitation:
    """
    Defines software citation data for use in mmcif files,
    but could be useful elsewhere.
    """

    journal_abbrev: str
    journal_volume: int
    journal_issue: int | str  # Uses empty string if no issue number.
    page_first: int
    page_last: int
    year: int
    title: str
    software_name: str
    software_version: str
    software_type: str
    software_classification: str
    software_description: str

    def mmcif_software_loop_data(self) -> tuple:
        return (
            self.software_name,
            self.software_version,
            self.software_type,
            self.software_classification,
            self.software_description,
        )

    def mmcif_citation_loop_data(self) -> tuple:
        return (
            self.journal_abbrev,
            self.journal_volume,
            self.journal_issue,
            self.page_first,
            self.page_last,
            self.year,
            self.title,
        )


dials_citation = SoftwareAndCitation(
    "Acta Cryst. D",
    74,
    2,
    85,
    97,
    2018,
    "DIALS: implementation and evaluation of a new integration package",
    "DIALS",
    dials_version,
    "package",
    "data processing",
    "Data processing and integration within the DIALS software package",
)

dials_scale_citation = SoftwareAndCitation(
    "Acta Cryst. D",
    76,
    4,
    385,
    399,
    2020,
    "Scaling diffraction data in the DIALS software package: algorithms and new approaches for multi-crystal scaling",
    "DIALS",
    dials_version,
    "program",
    "data scaling",
    "Data scaling and merging within the DIALS software package",
)

xds_citation = SoftwareAndCitation(
    "Acta Cryst. D",
    66,
    2,
    125,
    132,
    2010,
    "XDS",
    "XDS",
    "",
    "package",
    "data reduction",
    "Data integration with the XDS package",
)

ssx_citation = SoftwareAndCitation(
    "Methods Enzymol.",
    709,
    "",
    207,
    244,
    2024,
    "Processing serial synchrotron crystallography diffraction data with DIALS",
    "DIALS",
    dials_version,
    "program",
    "data reduction",
    "Data integration of still-shot data within the DIALS package",
)
