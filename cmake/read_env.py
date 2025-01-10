"""
Reads a libtbx_env file and extract information as a json string.

Some structures e.g. module dictionaries refer to the same object in multiple
places. These will be dumped multiple times. Anything that refers back to a
previous level will show as "__selfreference__". Relocatable paths are just
shown as the regular, joined path (from inspection all base off of the build
path anyway).
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from pathlib import Path
from types import ModuleType

BUILD_PATH = None
_all_relocatable_paths = set()
_all_absolute_paths = set()


def _read_obj(obj, prev=None):
    if prev is None:
        prev = []
    if obj in prev:
        return "__selfreference__"
    prev = list(prev) + [obj]

    if isinstance(obj, prop_object):
        dic = {name: _read_obj(val, prev) for name, val in obj.__dict__.items()}
        dic["__type__"] = obj._pathed_type
        return dic
    elif isinstance(obj, list):
        p = []
        for i in obj:
            p.append(_read_obj(i, prev))
        return p
    elif isinstance(obj, dict):
        return {a: _read_obj(b, prev) for a, b in obj.items()}
    else:
        return obj


class prop_object(object):
    """Object that can convert itself to a dictionary"""

    def to_dict(self):
        return _read_obj(self)


def pathed_prop_object(path):
    "Create a class that knows the path it's supposed to represent"

    class _pathed_prop_object(prop_object):
        """Object that can convert itself to a dictionary"""

        _pathed_type = path

    return _pathed_prop_object


class relocatable_path(object):
    def __repr__(self):
        _all_relocatable_paths.add(self)
        return os.path.normpath(os.path.join(str(self._anchor), self.relocatable))


class absolute_path(object):
    def __init__(self, path):
        # This init not used by unpickle - only for rewriting in here
        self._path = str(path)

    def __repr__(self):
        _all_absolute_paths.add(self)
        return self._path


def plainlify(thing):
    if isinstance(thing, (str, int, float, complex)):
        return thing
    if thing in (None, True, False):
        return thing
    if isinstance(thing, tuple):
        return tuple(map(plainlify, thing))
    if isinstance(thing, list):
        return list(map(plainlify, thing))
    if isinstance(thing, dict):
        return {plainlify(key): plainlify(value) for key, value in thing.items()}
    if isinstance(thing, set):
        return [plainlify(item) for item in thing]
    return str(thing)


def new_module(name, doc=None):
    """Create a new module and inject it into sys.modules"""
    m = ModuleType(name, doc)
    m.__file__ = name + ".py"
    sys.modules[name] = m
    return m


# Create the fake libtbx environment
libtbx = new_module("libtbx")
libtbx.env_config = new_module("libtbx.env_config")
libtbx.path = new_module("libtbx.path")
libtbx.env_config.environment = pathed_prop_object("libtbx.env_config.environment")
libtbx.env_config.build_options = pathed_prop_object("libtbx.env_config.build_options")
libtbx.env_config.module = pathed_prop_object("libtbx.env_config.module")
libtbx.path.relocatable_path = relocatable_path
libtbx.path.absolute_path = absolute_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read information from a libtbx_env file"
    )
    parser.add_argument("libtbx_env", type=Path)
    parser.add_argument(
        "--build-path",
        type=Path,
        help=(
            "The actual build path. If this disagrees with the environment, "
            "relative paths will be rewritten. This is to help match the "
            "libtbx behaviour under transplanted installs (e.g. conda)."
        ),
    )
    parser.add_argument(
        "--sys-prefix",
        type=Path,
        help=(
            "The sys.prefix for the distribution the libtbx_env sits within. "
            "If set, and the libtbx_env is an installed one, will rewrite internal prefix."
        ),
    )
    parser.add_argument(
        "--windows",
        action="store_true",
        help=(
            "Apply libtbx-installed-on-windows rewrite-rules, if the libtbx "
            "environment is an installed one. "
            "On windows, libtbx transplants the system root for some paths."
        ),
    )
    args = parser.parse_args()
    if not args.libtbx_env.is_file():
        sys.exit(f"Error: {args.libtbx_env} is not a file")
    # Use the libtbx_env path as the real build_path, for rewriting paths
    BUILD_PATH = args.build_path
    # Load the environment dump
    env = pickle.loads(args.libtbx_env.read_bytes())
    # Because of pickle, we need to __repr__ everything in order to register it
    plainlify(env.to_dict())
    # Rewrite the build folder, if we've been provided one
    if args.build_path:
        orig_build_path = str(env.build_path)
        new_build_path = str(args.build_path.resolve())
        # The relocatable path uses an "absolute path" as a parent. So,
        # make sure we translate all absolute paths that exist, and this
        # means that we've moved all relocatable paths.
        for a_path in _all_absolute_paths:
            if a_path._path == orig_build_path:
                a_path._path = new_build_path

    if env.installed and args.sys_prefix:
        # If this is an installed libtbx_env, then we have rules about rewriting it
        # this is... non-ideal, but since release libtbx_env files are broken on
        # windows, this is the best place to deal with it.
        new_prefix = absolute_path(args.sys_prefix.resolve())
        if args.windows:
            new_prefix = absolute_path(args.sys_prefix.resolve() / "library")
        # Replace the anchor for each of the repository paths
        for path in env.repository_paths:
            path._anchor = new_prefix
        env.bin_path._anchor = new_prefix
        env.exe_path._anchor = new_prefix
        env.include_path._anchor = new_prefix
        env.lib_path._anchor = new_prefix
        env.path_utility._anchor = new_prefix
        # Rewrite all of the module dist paths
        if args.windows:
            for module in env.module_list:
                for dist_path in module.dist_paths:
                    if dist_path is not None:
                        dist_path._anchor = new_prefix

    # Regenerate the plain dictionary now we've rewritten the paths
    d = plainlify(env.to_dict())
    print(json.dumps(d, indent=4))
