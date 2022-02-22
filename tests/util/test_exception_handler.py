from __future__ import annotations

import os
import sys

import pytest

import dials.util


class _SomeError(RuntimeError):
    pass


def test_crash_with_plain_text():
    with pytest.raises(_SomeError) as exc:
        with dials.util.show_mail_on_error():
            raise _SomeError("This is a string")
    assert "report this error" in str(exc.value)
    assert "This is a string" in str(exc.value)


def test_crash_with_unicode():
    with pytest.raises(_SomeError) as exc:
        with dials.util.show_mail_on_error():
            raise _SomeError("This is a 👹♔  𝓕คｎｃ¥ ｕηเᶜόＤ𝔼  🏆🍔 string")
    assert "report this error" in str(exc.value)
    assert "This is a" in str(exc.value)


def test_crash_with_bytestring():
    with pytest.raises(_SomeError) as exc:
        with dials.util.show_mail_on_error():
            raise _SomeError(b"This is a byte string")
    assert "report this error" in repr(exc.value)
    assert "This is a" in repr(exc.value)


def test_coloured_exit(monkeypatch):
    with pytest.raises(SystemExit) as e:
        with dials.util.make_sys_exit_red():
            sys.exit("Ohno")
    assert e.value.code == "Ohno"

    class _AttySysErr:
        def isatty(self):
            return True

    # Force colouring on
    monkeypatch.setattr(sys, "stderr", _AttySysErr())

    with pytest.raises(SystemExit) as e:
        with dials.util.make_sys_exit_red():
            sys.exit("Ohno")
    assert str(e.value.code).startswith("\033[")

    # and off again...
    monkeypatch.setitem(os.environ, "NO_COLOR", "1")
    with pytest.raises(SystemExit) as e:
        with dials.util.make_sys_exit_red():
            sys.exit("Ohno")
    assert e.value.code == "Ohno"
