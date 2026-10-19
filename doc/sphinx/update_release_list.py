from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime
from typing import Any, NamedTuple

import requests

DEFAULT_REPOSITORY = "dials/dials"

# Download kinds, in the order that the matrix columns are presented. Assets
# whose extension is not listed here are left out of the matrix entirely, which
# is how non-installer assets like bootstrap.py are dropped.
KINDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Graphical installer", (".pkg", ".exe", ".msi", ".dmg")),
    ("Shell installer", (".sh",)),
    ("Archive", (".tar.gz", ".tar.xz", ".tar.bz2", ".zip")),
)

# Extensions, longest first, so that e.g. ".tar.gz" wins over ".gz"
_EXTENSIONS = sorted(
    (ext for _, extensions in KINDS for ext in extensions), key=len, reverse=True
)
_KIND_FOR_EXTENSION = {ext: kind for kind, extensions in KINDS for ext in extensions}

# Both the historical "dials-v3-3-0-macosx.tar.gz" and the current
# "DIALS-3.29.0-MacOSX-arm64.pkg" asset naming
_ASSET_RE = re.compile(r"^dials-v?\d+[.-]\d+[.-]\d+-(?P<platform>.+)$", re.IGNORECASE)

_OPERATING_SYSTEMS = {
    "linux": ("Linux", 0),
    "macosx": ("macOS", 1),
    "macos": ("macOS", 1),
    "osx": ("macOS", 1),
    "windows": ("Windows", 2),
    "win": ("Windows", 2),
}


class Download(NamedTuple):
    platform: str
    kind: str
    extension: str
    url: str
    updated: datetime
    order: tuple[int, int]


def github_api_get(url: str, **kwargs) -> Any:
    """
    Fetch a GitHub API url, with authentication if available.

    If an environment variable GITHUB_TOKEN is set, thi will be used
    when accessing the releases API.
    """
    args = dict(kwargs)
    headers = {
        "Accept": "application/vnd.github.v3+json",
        **args.pop("headers", {}),
    }
    # If we've set an environment variable, use the token
    if token := os.environ.get("GITHUB_TOKEN", None):
        headers["Authorization"] = f"token {token}"
    req = requests.get(url, headers=headers, **args)
    req.raise_for_status()
    return req.json()


def _sort_tag_version(release_data: dict) -> tuple:
    """Sort a list of releases by tag version.

    Any non-numeric part of a tag is ignored rather than being an error, so
    that a repository with e.g. release-candidate tags can still be read.
    """
    parts = []
    for part in release_data["tag_name"].lstrip("v").split("."):
        number = re.match(r"\d+", part)
        parts.append(int(number.group()) if number else 0)
    return tuple(parts)


def _split_extension(filename: str) -> tuple[str, str] | None:
    """Split a known download extension off the end of an asset filename."""
    for extension in _EXTENSIONS:
        if filename.lower().endswith(extension):
            return filename[: -len(extension)], extension
    return None


def _architecture(operating_system: str, architecture: str) -> tuple[str, int]:
    """Convert an asset's architecture into a label, and an ordering for it."""
    architecture = architecture.lower()
    if architecture in {"x86_64", "amd64", "x64"}:
        # Apple describe this as "Intel", and it is the least common Mac now
        if operating_system == "macOS":
            return "Intel", 1
        return "x86-64", 0
    if architecture in {"arm64", "aarch64"}:
        if operating_system == "macOS":
            return "Apple silicon", 0
        return "ARM64", 1
    if not architecture:
        return "", 0
    return architecture, 2


def _parse_asset(asset: dict) -> Download | None:
    """Work out what platform and kind of download an asset is, if any."""
    split = _split_extension(asset["name"])
    if not split:
        return None
    stem, extension = split
    match = _ASSET_RE.match(stem)
    if not match:
        return None
    system, _, architecture = match.group("platform").partition("-")
    if system.lower() not in _OPERATING_SYSTEMS:
        return None
    operating_system, system_order = _OPERATING_SYSTEMS[system.lower()]
    architecture, architecture_order = _architecture(operating_system, architecture)

    return Download(
        platform=f"{operating_system} ({architecture})"
        if architecture
        else operating_system,
        kind=_KIND_FOR_EXTENSION[extension],
        extension=extension,
        url=asset["browser_download_url"],
        updated=datetime.strptime(asset["updated_at"], "%Y-%m-%dT%H:%M:%SZ"),
        order=(system_order, architecture_order),
    )


def _release_downloads(release_info: dict) -> list[Download]:
    """Extract the downloadable assets from a release, latest upload winning."""
    downloads: dict[tuple[str, str, str], Download] = {}
    for asset in release_info.get("assets", []):
        download = _parse_asset(asset)
        if not download:
            continue
        key = (download.platform, download.kind, download.extension)
        if key not in downloads or downloads[key].updated < download.updated:
            downloads[key] = download
    return list(downloads.values())


def _cell(downloads: list[Download]) -> list[str]:
    """Render the buttons for a single cell of the download matrix."""
    lines = []
    for download in sorted(downloads, key=lambda d: d.extension):
        if lines:
            lines.append("")
        lines.extend(
            [
                ".. button::",
                f"   :text: {download.extension}",
                f"   :link: {download.url}",
            ]
        )
    return lines


def _download_matrix(downloads: list[Download]) -> str:
    """Generate the reST markup for an OS/kind matrix of download buttons."""
    platforms = sorted({(d.order, d.platform) for d in downloads})
    kinds = [kind for kind, _ in KINDS if any(d.kind == kind for d in downloads)]
    if not platforms or not kinds:
        return ""

    rows = [[["Platform"]] + [[kind] for kind in kinds]]
    for _, platform in platforms:
        row = [[platform]]
        for kind in kinds:
            cell = [d for d in downloads if d.platform == platform and d.kind == kind]
            print(f"  {platform:<24} {kind:<20} {' '.join(d.extension for d in cell)}")
            row.append(_cell(cell))
        rows.append(row)

    table = [".. list-table::", "   :class: download-matrix", "   :header-rows: 1", ""]
    for row in rows:
        for index, cell in enumerate(row):
            # Cell contents all start in column 7, so that any continuation
            # lines can be indented consistently
            marker = "   * - " if index == 0 else "     - "
            table.append((marker + (cell[0] if cell else "")).rstrip())
            table.extend(("       " + line).rstrip() for line in cell[1:])
        table.append("")
    return "\n".join(table)


def run(args: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate the stable release download matrix for the "
        "installation page, from the assets attached to a Github release."
    )
    parser.add_argument(
        "--repo",
        default=DEFAULT_REPOSITORY,
        metavar="OWNER/NAME",
        help=f"Repository to read releases from. Default: {DEFAULT_REPOSITORY}",
    )
    parser.add_argument(
        "--release",
        metavar="TAG",
        help="Use this release instead of the highest-numbered one, matched "
        "against the release tag or name. Lets you preview the page for a "
        "release before it is the current one. Unpublished draft releases "
        "are only visible if GITHUB_TOKEN is set.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=os.path.join(os.path.dirname(__file__), "installation.stable_release"),
        metavar="FILE",
        help="Where to write the generated reST. Pass an alternative path to "
        "preview a release without overwriting the checked-in page.",
    )
    options = parser.parse_args(args)

    print(f"Checking DIALS release status ({options.repo})")
    releases_info = sorted(
        github_api_get(f"https://api.github.com/repos/{options.repo}/releases"),
        reverse=True,
        key=_sort_tag_version,
    )
    if not releases_info:
        sys.exit(f"Error: {options.repo} has no releases")

    if options.release:
        matching = [
            x for x in releases_info if options.release in {x["tag_name"], x["name"]}
        ]
        if not matching:
            sys.exit(
                f"Error: No release '{options.release}' in {options.repo}. Available: "
                + ", ".join(x["tag_name"] for x in releases_info)
            )
        release_info = matching[0]
    else:
        release_info = releases_info[0]

    # A draft release does not necessarily have a name yet
    name = release_info["name"] or release_info["tag_name"]
    print(f"Using release: {name} ({release_info['tag_name']})")

    with open(options.output, "w") as release:
        caption = f"Stable Release: {name}"
        release.write(
            f"""
{caption}
{"=" * len(caption)}

The current stable release can be downloaded from `Github <{release_info["html_url"]}>`_,
where you can also find further release notes.

""".lstrip()
        )
        release.write(_download_matrix(_release_downloads(release_info)))
    print(f"Written to {options.output}")


if __name__ == "__main__":
    run()
