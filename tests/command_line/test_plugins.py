from __future__ import annotations


def test_plugin_setup_is_valid():
    import dials.command_line.plugins

    assert dials.command_line.plugins.installation_is_valid()
