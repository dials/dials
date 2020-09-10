# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import pytest
import six

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
            raise _SomeError(u"This is a 👹♔  𝓕คｎｃ¥ ｕηเᶜόＤ𝔼  🏆🍔 string")
    assert u"report this error" in six.text_type(exc.value)
    assert u"This is a" in six.text_type(exc.value)


def test_crash_with_bytestring():
    with pytest.raises(_SomeError) as exc:
        with dials.util.show_mail_on_error():
            raise _SomeError(b"This is a byte string")
    assert "report this error" in repr(exc.value)
    assert "This is a" in repr(exc.value)
