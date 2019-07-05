from __future__ import absolute_import, division, print_function

import pkg_resources
import sys

BOLD = "\033[1m"
RED = "\033[1;31m"
GREEN = "\033[32m"
NC = "\033[0m"


def read_entry_point(entry_point):
    """
    Finds all registered plugins for a given entry point.

    :return: A dictionary of entry point plugins
    """
    return {e.name: e for e in pkg_resources.iter_entry_points(entry_point)}


def installation_is_valid():
    """
    Verifies that all required plugins are present.

    :return: True if so, False otherwise
    """
    for ep, ep_dict in known_entry_points.items():
        required_plugins = set(ep_dict.get("required", []))
        if not required_plugins:
            continue
        plugins = read_entry_point(ep)
        missing_plugins = required_plugins - set(plugins)
        if missing_plugins:
            return False
    return True


known_entry_points = {
    "dxtbx.profile_model": {
        "description": "profile models",
        "required": ["gaussian_rs"],
    },
    "dxtbx.scaling_model_ext": {"description": "scaling models"},
    "dials.index.basis_vector_search_strategy": {
        "description": "Basis vector search strategies",
        "required": ["fft1d", "fft3d", "real_space_grid_search"],
    },
}

if __name__ == "__main__":
    for ep, ep_dict in known_entry_points.items():
        print(
            "{BOLD}{ep}{NC}  {ep_dict[description]}".format(
                BOLD=BOLD, NC=NC, ep=ep, ep_dict=ep_dict
            )
        )
        plugins = read_entry_point(ep)
        for p in sorted(plugins):
            print(
                " {GREEN}{p} {NC}({pp.module_name} via {BOLD}{pp.dist.project_name}{NC} {pp.dist.version})".format(
                    BOLD=BOLD, p=p, GREEN=GREEN, NC=NC, pp=plugins[p]
                )
            )
        required_plugins = set(ep_dict.get("required", []))
        missing_plugins = required_plugins - set(plugins)
        if missing_plugins:
            print(" {RED}  --- missing required plugins:{NC}".format(NC=NC, RED=RED))
            for p in missing_plugins:
                print(" {RED}{p}{NC}".format(NC=NC, p=p, RED=RED))
        print()
    sys.exit(not installation_is_valid())
