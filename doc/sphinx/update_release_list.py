from __future__ import annotations

import os
import re
from datetime import datetime
from typing import Any

import requests


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


def _download_button(text, version, link):
    print(f"  {version} {text} -> {link}")
    return f".. button::\n   :text: DIALS {version} {text}\n   :link: {link}\n\n"


def _sort_tag_version(release_data: dict) -> tuple:
    """Sort a list of releases by tag version."""
    return tuple(int(x) for x in release_data["tag_name"].lstrip("v").split("."))


def _release_asset_buttons(release_info: dict) -> list[str]:
    """Generate the reST markup for the download buttons."""
    assets = {}
    for a in release_info.get("assets", []):
        tag = re.search("dials-v([^-]+)-([^-]+)-([^-]+)-(.+)", a["name"])
        if tag:
            shortname = tag.group(4)
            version = ".".join(tag.group(1, 2, 3))
            last_update = datetime.strptime(
                a["updated_at"], "%Y-%m-%dT%H:%M:%SZ"
            )  # - datetime(1970,1,1)).total_seconds()
            if shortname not in assets or assets[shortname][0] < last_update:
                assets[shortname] = (
                    last_update,
                    version,
                    a.get("browser_download_url"),
                )

    long_names = {
        "macosx.pkg": "Mac installer",
        "macosx.tar.gz": "Mac tar archive",
        "linux-x86_64.tar.xz": "Linux installer",
    }

    return [
        _download_button(long_names.get(asset, asset), _version, link)
        for asset, (_, _version, link) in assets.items()
    ]


if __name__ == "__main__":
    release_file = os.path.join(
        os.path.dirname(__file__), "installation.stable_release"
    )
    print("Checking DIALS release status")
    releases_info = sorted(
        github_api_get("https://api.github.com/repos/dials/dials/releases"),
        reverse=True,
        key=_sort_tag_version,
    )

    latest_release = releases_info[0]
    print(f"Latest release is: {latest_release['name']} ({latest_release['tag_name']})")

    with open(release_file, "w") as release:
        caption = f"Stable Release: {latest_release['name']}"
        release.write(
            f"""
{caption}
{'=' * len(caption)}

The current stable release can be downloaded from `Github <{latest_release['html_url']}>`_,
where you can also find further release notes.

""".lstrip()
        )
        release.write("".join(sorted(_release_asset_buttons(latest_release))))
