# This file contains Python 3-only syntax, so must not be imported by a
# Python 2 pytest instance


def test_ignore_python3_syntax_files_with_python2_pytest():
    x = 1
    assert f"{x}" == "1"
