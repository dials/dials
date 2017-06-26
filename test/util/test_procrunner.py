from __future__ import absolute_import, division

import dials.util.procrunner
import mock
import pytest

@pytest.mark.skipif(dials.util.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('dials.util.procrunner._NonBlockingStreamReader')
@mock.patch('dials.util.procrunner.time')
@mock.patch('dials.util.procrunner.subprocess')
@mock.patch('dials.util.procrunner.Pipe')
def test_run_command_aborts_after_timeout(mock_pipe, mock_subprocess, mock_time, mock_streamreader):
  mock_pipe.return_value = mock.Mock(), mock.Mock()
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

  actual = dials.util.procrunner.run_process(command, 0.5, False,
               callback_stdout=mock.sentinel.callback_stdout, callback_stderr=mock.sentinel.callback_stderr)

  assert mock_subprocess.Popen.called
  mock_streamreader.assert_has_calls([mock.call(stream_stdout, output=mock.ANY, debug=mock.ANY, notify=mock.ANY, callback=mock.sentinel.callback_stdout),
                                      mock.call(stream_stderr, output=mock.ANY, debug=mock.ANY, notify=mock.ANY, callback=mock.sentinel.callback_stderr)],
                                     any_order=True)
  assert not mock_process.terminate.called
  assert not mock_process.kill.called
  assert actual == expected

@mock.patch('dials.util.procrunner.select')
def test_nonblockingstreamreader_can_read(mock_select):
  import time
  class _stream:
    def __init__(self):
      self.data = ""
      self.closed = False
    def write(self, string):
      self.data = self.data + string
    def read(self, n):
      if self.closed:
        return ""
      if self.data == "":
        time.sleep(0.01)
        return ""
      if (len(self.data) < n):
        data = self.data
        self.data = ""
      else:
        data = self.data[:n]
        self.data = self.data[n:]
      return data
    def close(self):
      self.closed=True
  teststream = _stream()

  def select_replacement(rlist, wlist, xlist, timeout):
    assert teststream in rlist
    if teststream.closed:
      return ([teststream], [], [])
    if teststream.data == "":
      return ([], [], [])
    return ([teststream], [], [])
  mock_select.select = select_replacement

  streamreader = dials.util.procrunner._NonBlockingStreamReader(teststream, output=False)
  assert not streamreader.has_finished()
  time.sleep(0.1)
  testdata = "abc\n" * 1024
  teststream.write(testdata)
  time.sleep(0.2)
  teststream.close()
  time.sleep(0.1)

  assert streamreader.has_finished()
  output = streamreader.get_output()
  assert len(output) == len(testdata)
  assert output == testdata

def test_lineaggregator_aggregates_data():
  callback = mock.Mock()
  aggregator = dials.util.procrunner._LineAggregator(callback=callback)

  aggregator.add('some')
  aggregator.add('string')
  callback.assert_not_called()
  aggregator.add("\n")
  callback.assert_called_once_with('somestring')
  callback.reset_mock()
  aggregator.add('more')
  aggregator.add('stuff')
  callback.assert_not_called()
  aggregator.flush()
  callback.assert_called_once_with('morestuff')
