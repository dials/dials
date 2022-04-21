from __future__ import annotations

import io
from unittest import mock

import pytest

import boost_adaptbx.boost.python

from dials.util.ext import ostream, streambuf

ext = boost_adaptbx.boost.python.import_ext("dials_util_streambuf_test_ext")


class io_test_case:
    phrase = b"Coding should be fun"
    #          01234567890123456789

    def run(self):
        m = streambuf.default_buffer_size
        for n in range(50, 0, -1):
            streambuf.default_buffer_size = n
            self.exercise_read_failure()
            self.exercise_write_failure()
            self.exercise_read()
            self.exercise_write()
            self.exercise_seek_and_read()
            self.exercise_partial_read()
            self.exercise_write_and_seek()
        streambuf.default_buffer_size = m

    def exercise_read(self):
        self.create_file_object(mode="rb")
        words = ext.read_word(streambuf(self.file_object))
        assert words == b"Coding, should, be, fun, [ fail, eof ]"
        self.file_object.close()

    def exercise_partial_read(self):
        self.create_file_object(mode="rb")
        words = ext.partial_read(streambuf(self.file_object))
        assert words == b"Coding, should, "
        trailing = self.file_object.read()
        assert trailing == b" be fun"
        self.file_object.close()

    def exercise_read_failure(self):
        self.create_file_object(mode="rb")
        self.file_object.close()
        with pytest.raises(ValueError):
            ext.read_word(streambuf(self.file_object))
        self.file_object.close()

    def exercise_write(self):
        self.create_file_object(mode="wb")
        report = ext.write_word(ostream(self.file_object))
        assert report == b""
        assert self.file_content() == b"2 times 1.6 equals 3.2"
        self.file_object.close()

    def exercise_seek_and_read(self):
        self.create_file_object(mode="rb")
        instrumented_file = mock.Mock(spec=self.file_object, wraps=self.file_object)
        words = ext.read_and_seek(streambuf(instrumented_file))
        assert words == b"should, should, uld, ding, fun, [ eof ]"
        n = streambuf.default_buffer_size
        soughts = instrumented_file.seek.call_args_list
        # stringent tests carefully crafted to make sure the seek-in-buffer
        # optimisation works as expected
        # C.f. the comment in the C++ function actual_write_test
        assert soughts[-1] == mock.call(-4, 2)
        soughts = soughts[:-1]
        if n >= 14:
            assert soughts == []
        else:
            assert soughts[0] == mock.call(6, 0)
            soughts = soughts[1:]
            if 8 <= n <= 13:
                assert len(soughts) == 1 and self.only_seek_cur(soughts)
            elif n == 7:
                assert len(soughts) == 2 and self.only_seek_cur(soughts)
            else:
                assert soughts[0] == mock.call(6, 0)
                soughts = soughts[1:]
                assert self.only_seek_cur(soughts)
                if n == 4:
                    assert len(soughts) == 1
                else:
                    assert len(soughts) == 2
        self.file_object.close()

    def exercise_write_and_seek(self):
        self.create_file_object(mode="wb")
        instrumented_file = mock.Mock(spec=self.file_object, wraps=self.file_object)
        report = ext.write_and_seek(ostream(instrumented_file))
        assert report == b""
        expected = b"1000 times 1000 equals 1000000"
        assert self.file_content() == expected
        assert self.file_object.tell() == 9
        if streambuf.default_buffer_size >= 30:
            assert instrumented_file.write.call_count == 1
        self.file_object.close()

    @staticmethod
    def only_seek_cur(seek_calls):
        return all(call == mock.call(mock.ANY, 1) for call in seek_calls)


class bytesio_test_case(io_test_case):
    def exercise_write_failure(self):
        pass

    def create_file_object(self, mode):
        if mode == "rb":
            self.file_object = io.BytesIO(self.phrase)
        elif mode == "wb":
            self.file_object = io.BytesIO()
        else:
            raise NotImplementedError("Internal error in the test code")

    def file_content(self):
        return self.file_object.getvalue()


class mere_file_test_case(io_test_case):
    def exercise_write_failure(self):
        self.create_file_object(mode="rb")
        with pytest.raises(IOError):
            ext.write_word(streambuf(self.file_object))
        self.file_object.close()

    def create_file_object(self, mode):
        f = open("tmp_tst_python_streambuf", "wb")
        if mode.find("rb") > -1:
            f.write(self.phrase)
        f.close()
        self.file_object = open(f.name, mode)

    def file_content(self):
        self.file_object.flush()
        with open(self.file_object.name, "rb") as fh:
            result = fh.read()
        return result


def test_with_bytesio():
    bytesio_test_case().run()


@pytest.mark.xfail(
    "os.name == 'nt'", reason="crashes python process on Windows", run=False
)
def test_with_file(tmpdir):
    with tmpdir.as_cwd():
        mere_file_test_case().run()
