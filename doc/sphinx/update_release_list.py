import json
import os
import re
from datetime import datetime
from urllib.request import Request, urlopen


def _download_button(text, version, link):
    print(f"  {version} {text} -> {link}")
    return f".. button::\n   :text: DIALS {version} {text}\n   :link: {link}\n\n"


if __name__ == "__main__":
    release_file = os.path.join(
        os.path.dirname(__file__), "installation.stable_release"
    )
    print("Checking DIALS release status")
    url_request = Request("https://api.github.com/repos/dials/dials/releases/latest")
    with urlopen(url_request) as conn:
        release_info = json.load(conn)

    with open(release_file, "w") as release:
        caption = "Stable Release"
        if "name" in release_info:
            caption = caption + ": " + release_info["name"]
            print("Most recent major DIALS release is:", release_info["name"])
        else:
            print("Could not determine most recent major DIALS release")
        release.write(caption + "\n" + "=" * len(caption) + "\n\n")

        release.write(
            """
The current stable release can be downloaded from `Github <https://github.com/dials/dials/releases/latest>`_,
where you can also find further `release notes <https://github.com/dials/dials/releases/latest>`_.

""".lstrip()
        )

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

        buttons = [
            _download_button(long_names.get(asset, asset), _version, link)
            for asset, (_, _version, link) in assets.items()
        ]

        release.write("".join(sorted(buttons)))
