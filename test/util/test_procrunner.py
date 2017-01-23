from __future__ import absolute_import, division

import dials.util.procrunner
import mock
import pytest

@pytest.mark.skipif(dials.util.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('dials.util.procrunner._NonBlockingStreamReader')
@mock.patch('dials.util.procrunner.time')
@mock.patch('dials.util.procrunner.subprocess')
@mock.patch('dials.util.procrunner.Queue')
def test_run_command_aborts_after_timeout(mock_queue, mock_subprocess, mock_time, mock_streamreader):
  mock_process = mock.Mock()
  mock_process.returncode = None
  mock_subprocess.Popen.return_value = mock_process
  task = ['___']

  with pytest.raises(RuntimeError):
    dials.util.procrunner.run_process(task, -1, False)

  assert mock_subprocess.Popen.called
  assert mock_process.terminate.called
  assert mock_process.kill.called


@pytest.mark.skipif(dials.util.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('dials.util.procrunner._NonBlockingStreamReader')
@mock.patch('dials.util.procrunner.subprocess')
def test_run_command_runs_command_and_directs_pipelines(mock_subprocess, mock_streamreader):
  (mock_stdout, mock_stderr) = (mock.Mock(), mock.Mock())
  mock_stdout.get_output.return_value = mock.sentinel.proc_stdout
  mock_stderr.get_output.return_value = mock.sentinel.proc_stderr
  (stream_stdout, stream_stderr) = (mock.sentinel.stdout, mock.sentinel.stderr)
  mock_process = mock.Mock()
  mock_process.stdout = stream_stdout
  mock_process.stderr = stream_stderr
  mock_process.returncode = 99
  command = ['___']
  def streamreader_processing(*args, **kwargs):
    return {(stream_stdout,): mock_stdout, (stream_stderr,): mock_stderr}[args]
  mock_streamreader.side_effect = streamreader_processing
  mock_subprocess.Popen.return_value = mock_process

  expected = {
    'stderr': mock.sentinel.proc_stderr,
    'stdout': mock.sentinel.proc_stdout,
    'exitcode': mock_process.returncode,
    'command': command,
    'runtime': mock.ANY,
    'timeout': False,
    'time_start': mock.ANY,
    'time_end': mock.ANY
  }

  actual = dials.util.procrunner.run_process(command, 0.5, False)

  assert mock_subprocess.Popen.called
  mock_streamreader.assert_has_calls([mock.call(stream_stdout, output=mock.ANY, debug=mock.ANY, notify=mock.ANY),
                                      mock.call(stream_stderr, output=mock.ANY, debug=mock.ANY, notify=mock.ANY)], any_order=True)
  assert not mock_process.terminate.called
  assert not mock_process.kill.called
  assert actual == expected


def test_nonblockingstreamreader_can_read():
  import time
  class _stream:
    def __init__(self):
      self.data = []
      self.closed = False
    def write(self, string):
      self.data.append(string)
    def readline(self):
      while (len(self.data) == 0) and not self.closed:
        time.sleep(0.2)
      return self.data.pop(0) if len(self.data) > 0 else ''
    def close(self):
      self.closed=True

  teststream = _stream()
  testdata = ['a', 'b', 'c']

  streamreader = dials.util.procrunner._NonBlockingStreamReader(teststream, output=False)
  for d in testdata:
    teststream.write(d)
  assert not streamreader.has_finished()

  teststream.close()
  time.sleep(0.4)

  assert streamreader.has_finished()
  assert streamreader.get_output() == ''.join(testdata)

