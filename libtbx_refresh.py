from __future__ import absolute_import, division, print_function

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
    }
)


try:
    from dials.util.version import dials_version

    print(dials_version())
except Exception:
    pass


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

    commands_dir = os.path.join(dist_path, "command_line")
    command_list = []
    print("Identifying autocompletable commands:", end=" ")
    for f in sorted(os.listdir(commands_dir)):
        if not f.startswith("_") and f.endswith(".py"):
            if (
                "DIALS_ENABLE_COMMAND_LINE_COMPLETION"
                in open(os.path.join(commands_dir, f)).read()
            ):
                command_name = "dials.%s" % f[:-3]
                print(command_name, end=" ")
                command_list.append(command_name)
    print()

    # Generate the autocompletion SConscript.
    with open(os.path.join(output_directory, "SConscript"), "w") as builder:
        builder.write(
            """Import("env")
import os.path
import libtbx.load_env
def dispatcher_outer(name):
  return os.path.join(libtbx.env.under_build('bin'), name)
def dispatcher_inner(name):
  return os.path.join(libtbx.env.dist_path('dials'), 'command_line', '%%s.py' %% name.partition('.')[2])
env.Append( BUILDERS={'AutoComplete': Builder(action='-$SOURCE --export-autocomplete-hints > $TARGET')} )
env['ENV']['DIALS_NOBANNER'] = '1'
for cmd in [%s]:
  ac = env.AutoComplete(cmd, [dispatcher_outer(cmd), dispatcher_inner(cmd)])
  Requires(ac, Dir(libtbx.env.under_build('lib')))
  Depends(ac, os.path.join(libtbx.env.dist_path('dials'), 'util', 'options.py'))
  Depends(ac, os.path.join(libtbx.env.dist_path('dials'), 'util', 'autocomplete.sh'))
"""
            % ", ".join(["'%s'" % cmd for cmd in command_list])
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
    for f in os.listdir(build_path):
        if f.startswith("setpath") and f.endswith(".sh"):
            original_file = open(os.path.join(build_path, f)).read()
            if not "DIALS_ENABLE_COMMAND_LINE_COMPLETION" in original_file:
                marker = "\nexport PATH\n"
                original_position = original_file.find(marker)
                if original_position >= 0:
                    print(f, end=" ")
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
                    with open(os.path.join(build_path, f), "w") as script:
                        script.write(
                            original_file[:insert_position]
                            + added_script
                            + original_file[insert_position:]
                        )
    print()


_install_dials_autocompletion()
