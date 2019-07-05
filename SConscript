from __future__ import absolute_import, division, print_function

import libtbx.load_env
import os
import platform

Import("env_etc")

env_etc.dials_dist = libtbx.env.dist_path("dials")
env_etc.dials_include = os.path.dirname(env_etc.dials_dist)
if not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include"):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)
    include_paths = [
        env_etc.libtbx_include,
        env_etc.scitbx_include,
        env_etc.cctbx_include,
        env_etc.rstbx_include,
        env_etc.boost_include,
        env_etc.boost_adaptbx_include,
        env_etc.python_include,
        env_etc.dxtbx_include,
        env_etc.dials_include,
    ]
    msgpack = os.path.join(env_etc.dials_include, "msgpack-3.1.1", "include")
    if os.path.exists(str(msgpack)):
        include_paths.append(msgpack)
    else:
        print("msgpack header files not installed, please run libtbx.install msgpack")
    env_etc.include_registry.append(env=env, paths=include_paths)
    env.Append(
        LIBS=env_etc.libm
        + ["scitbx_boost_python", "boost_python", "boost_thread", "cctbx"]
    )
    env.SConscript("model/SConscript", exports={"env": env})
    env.SConscript("array_family/SConscript", exports={"env": env})
    env.SConscript("algorithms/SConscript", exports={"env": env})
    env.SConscript("pychef/SConscript", exports={"env": env})
    env.SConscript("viewer/SConscript", exports={"env": env})
    # env.SConscript('nexus/SConscript', exports={ 'env' : env })
    env.SConscript("util/SConscript", exports={"env": env})

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
    env.SConscript("test/SConscript", exports={"env": env})
