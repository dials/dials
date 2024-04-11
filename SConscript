import libtbx.load_env
import os
import platform
from libtbx.env_config import get_boost_library_with_python_version
from pathlib import Path

Import("env_etc")
Import("env_no_includes_boost_python_ext")


def _get_node_child(node, *path):
    """Walk a SCons node tree to a particular named child"""
    if not path:
        return node
    return _get_node_child(
        [x for x in node.children() if x.name == path[0]][0], *path[1:]
    )


# Problem: gltbx fails compilation on RHEL8 conda-forge compilers.
#
# Reason:
# - The conda-forge sysroot stdlib system uses #include_next in order
#   to chain to a wrapped stdlib.h.
# - gltbx adds the system /usr/include and /usr/local/include multiple
#   times to the include search path.
# - When the c-f sysroot stdlib tries to #include_next, instead of the
#   REAL sysroot stdlib.h, it picks up the system stdlib.h.
# - This worked on RHEL7- the system stdlib was old enough that it
#   apparently didn't use this mechanism.
# - Since the RHEL8 system stdlib is recent enough to use this mechanism,
#   it thinks it is being #included twice, and so hits the usual guards.
#
# This "System include path" was added originally to work around the fact
# that OpenGL is a complicated dependency and thus you need to use the
# library from the system. This problem can also be solved using the
# Conda-forge CDT system dependencies.
#
# So, IF we have the CDT OpenGL packages, attempt to fix the gltbx
# build. We think this is safe here because:
# - We never use the layered tbx environment any more (ran into issues)
#   so either we are building CCTBX, or we aren't using SCons.
# - These paths won't get included on other platforms.
#
gltbx_env = _get_node_child(
    env_no_includes_boost_python_ext.fs.Top, "lib", "gltbx_fonts_ext.so"
).env
while "/usr/include" in gltbx_env["CPPPATH"]:
    gltbx_env["CPPPATH"].remove("/usr/include")
while "/usr/local/include" in gltbx_env["CPPPATH"]:
    gltbx_env["CPPPATH"].remove("/usr/local/include")


env_etc.dials_dist = os.path.join(libtbx.env.dist_path("dials"), "src", "dials")
env_etc.dials_include = os.path.dirname(env_etc.dials_dist)
if not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include"):
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)

    system_includes = [x for x in env_etc.conda_cpppath if x] if libtbx.env.build_options.use_conda else []
    system_includes.append(str(Path(env_etc.scitbx_dist).parent))
    env.Append(CXXFLAGS=[f"-isystem{x}" for x in system_includes])
    env.Append(SHCXXFLAGS=[f"-isystem{x}" for x in system_includes])

    include_paths = [
        env_etc.libtbx_include,
        env_etc.scitbx_include,
        env_etc.cctbx_include,
        env_etc.ccp4io_include,
        env_etc.rstbx_include,
        env_etc.boost_include,
        env_etc.boost_adaptbx_include,
        env_etc.python_include,
        env_etc.dxtbx_include,
        env_etc.dials_include,
    ]

    # Handle cctbx bootstrap builds that pull a fixed msgpack version into modules/
    msgpack = Path(libtbx.env.dist_path("dials")).parent / "msgpack-3.1.1" / "include"
    if msgpack.is_dir():
        include_paths.append(str(msgpack))

    if libtbx.env.build_options.use_conda:
        boost_python = get_boost_library_with_python_version(
            "boost_python", env_etc.conda_libpath
        )
        env.Append(LIBPATH=env_etc.conda_libpath)
        include_paths.extend(env_etc.conda_cpppath)
    else:
        boost_python = "boost_python"
    env_etc.include_registry.append(env=env, paths=include_paths)
    env.Append(
        LIBS=env_etc.libm
        + ["scitbx_boost_python", boost_python, "boost_thread", "cctbx"],
    )

    # Fix the build environment so that it doesn't break on modern C++
    for path in list(env["CPPPATH"]):
        if "msvc9.0_include" in path:
            env["CPPPATH"].remove(path)

    # Fix compilation errors on windows, caused by function redefinition
    # See: https://github.com/boostorg/system/issues/32#issuecomment-462912013
    if env_etc.compiler == "win32_cl":
        env.Append(CPPDEFINES="HAVE_SNPRINTF")

    env.SConscript("src/dials/model/SConscript", exports={"env": env})
    env.SConscript("src/dials/array_family/SConscript", exports={"env": env})
    env.SConscript("src/dials/algorithms/SConscript", exports={"env": env})
    env.SConscript("src/dials/pychef/SConscript", exports={"env": env})
    env.SConscript("src/dials/viewer/SConscript", exports={"env": env})
    # env.SConscript('src/dials/nexus/SConscript', exports={ 'env' : env })
    env.SConscript("src/dials/util/SConscript", exports={"env": env})

    autocomplete_scons = os.path.join(
        libtbx.env.under_build(os.path.join("dials", "autocomplete")), "SConscript"
    )
    if not any(platform.win32_ver()) and os.path.isfile(autocomplete_scons):
        env.SConscript(autocomplete_scons, exports={"env": env})

    #
    # NOTE: This must go at the bottom. The LIBS are replaced with an empty
    # list. This is done because errors occur when building the tests if it
    # isn't done. Replacing the libs afterwards still results in those errors.
    #
    env.SConscript("tests/SConscript", exports={"env": env})
