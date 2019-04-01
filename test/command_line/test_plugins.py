from __future__ import absolute_import, division, print_function


def test_plugin_setup_is_valid():
    import dials.command_line.plugins

    assert dials.command_line.plugins.installation_is_valid()
