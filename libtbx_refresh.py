from __future__ import absolute_import, division, print_function

import dials.precommitbx.nagger
import libtbx.pkg_utils

libtbx.pkg_utils.define_entry_points(
    {
        "dxtbx.profile_model": [
            "gaussian_rs = dials.extensions.gaussian_rs_profile_model_ext:GaussianRSProfileModelExt"
        ],
        "dxtbx.scaling_model_ext": [
            "physical = dials.algorithms.scaling.model.scaling_model_ext:PhysicalScalingModelExt",
            "KB = dials.algorithms.scaling.model.scaling_model_ext:KBScalingModelExt",
            "array = dials.algorithms.scaling.model.scaling_model_ext:ArrayScalingModelExt",
        ],
        "dials.index.basis_vector_search_strategy": [
            "fft1d = dials.algorithms.indexing.basis_vector_search.strategies:FFT1D",
            "fft3d = dials.algorithms.indexing.basis_vector_search.strategies:FFT3D",
            "real_space_grid_search = dials.algorithms.indexing.basis_vector_search.strategies:RealSpaceGridSearch",
        ],
    }
)


try:
    from dials.util.version import dials_version

    print(dials_version())
except Exception:
    pass

dials.precommitbx.nagger.nag()


def _install_dials_autocompletion():
    """generate bash.sh and SConscript file in /build/dials/autocomplete"""
    import libtbx.load_env
    import os  # required due to cctbx weirdness

    # Find the dials source directory
    dist_path = libtbx.env.dist_path("dials")

    # Set the location of the output directory
    output_directory = libtbx.env.under_build(os.path.join("dials", "autocomplete"))
    try:
        os.makedirs(output_directory)
    except OSError:
        pass

    # Build a list of autocompleteable commands
    commands_dir = os.path.join(dist_path, "command_line")
    command_list = []
    for filename in sorted(os.listdir(commands_dir)):
        if not filename.startswith("_") and filename.endswith(".py"):
            # Check if this file marks itself as completable
            with open(os.path.join(commands_dir, filename)) as f:
                if "DIALS_ENABLE_COMMAND_LINE_COMPLETION" in f.read():
                    command_name = "dials.%s" % filename[:-3]
                    command_list.append(command_name)
    print("Identified autocompletable commands: " + " ".join(command_list))

    # Generate the autocompletion SConscript.
    with open(os.path.join(output_directory, "SConscript"), "w") as builder:
        builder.write(
            """import os.path
import libtbx.load_env\n
Import("env")\n\n
def dispatcher_outer(name):
    return os.path.join(libtbx.env.under_build("bin"), name)\n\n
def dispatcher_inner(name):
    return os.path.join(
        libtbx.env.dist_path("dials"), "command_line", "%s.py" % name.partition(".")[2]
    )\n\n
env.Append(
    BUILDERS={{
        "AutoComplete": Builder(action="-$SOURCE --export-autocomplete-hints > $TARGET")
    }}
)
env["ENV"]["DIALS_NOBANNER"] = "1"
for cmd in [
{}
]:
    ac = env.AutoComplete(cmd, [dispatcher_outer(cmd), dispatcher_inner(cmd)])
    Requires(ac, Dir(libtbx.env.under_build("lib")))
    Depends(ac, os.path.join(libtbx.env.dist_path("dials"), "util", "options.py"))
    Depends(ac, os.path.join(libtbx.env.dist_path("dials"), "util", "autocomplete.sh"))
""".format(
                "\n".join(['    "{}",'.format(cmd) for cmd in command_list])
            )
        )

    # Generate a bash script activating command line completion for each relevant command
    with open(os.path.join(output_directory, "bash.sh"), "w") as script:
        script.write("type compopt &>/dev/null && {\n")
        for cmd in command_list:
            script.write(" complete -F _dials_autocomplete %s\n" % cmd)
        script.write("}\n")
        script.write("type compopt &>/dev/null || {\n")
        for cmd in command_list:
            script.write(" complete -o nospace -F _dials_autocomplete %s\n" % cmd)
        script.write("}\n")

    # Find the dials build directory
    build_path = abs(libtbx.env.build_path)

    # Permanently install the autocompletion script into setpaths-scripts.
    print("Installing autocompletion script into:", end=" ")
    for filename in os.listdir(build_path):
        if filename.startswith("setpath") and filename.endswith(".sh"):
            with open(os.path.join(build_path, filename)) as f:
                original_file = f.read()
            if "DIALS_ENABLE_COMMAND_LINE_COMPLETION" not in original_file:
                marker = "\nexport PATH\n"
                original_position = original_file.find(marker)
                if original_position >= 0:
                    print(filename, end=" ")
                    insert_position = original_position + len(marker)
                    added_script = (
                        "# DIALS_ENABLE_COMMAND_LINE_COMPLETION\n"
                        '[ -n "$BASH_VERSION" ] && {\n'
                        " source $(libtbx.find_in_repositories dials/util/autocomplete.sh) && source %s || echo dials command line completion not available\n"
                        "}\n"
                        % (
                            os.path.join(
                                "$LIBTBX_BUILD", "dials", "autocomplete", "bash.sh"
                            )
                        )
                    )
                    with open(os.path.join(build_path, filename), "w") as script:
                        script.write(
                            original_file[:insert_position]
                            + added_script
                            + original_file[insert_position:]
                        )
    print()


_install_dials_autocompletion()
