# coding: utf-8

"""
Customisable progressbar decorator for iterators.

Utilizes the progress bar from the tqdm package, with dials
customizations:
  - Specify a minimum total time with min_time - if the estimated total
    time for completion is below this period, the progress bar will not
    be shown. Defaults to 10 seconds.
  - By default, resize when terminal is resized (dynamic-columns)
  - By default, disables the progress bar for non-tty output
  - By default, the progress bar will be removed after completion

Usage:
  >>> from dials.util import progress
  >>> for i in progress(range(10)):
  ...     ...

See https://github.com/tqdm/tqdm for more in-depth usage and options.
"""

from __future__ import print_function, division, absolute_import

from tqdm import tqdm


class progress(tqdm):
  def __init__(self, *args, **kwargs):
    """
    Parameters
    ----------
    iterable  : iterable, optional
        Iterable to decorate with a progressbar.
        Leave blank to manually manage the updates.
    desc  : str, optional
        Prefix for the progressbar.
    total  : int, optional
        The number of expected iterations. If unspecified,
        len(iterable) is used if possible. As a last resort, only basic
        progress statistics are displayed (no ETA, no progressbar).
        If `gui` is True and this parameter needs subsequent updating,
        specify an initial arbitrary large positive integer,
        e.g. int(9e9).
    leave  : bool, optional
        If True, keeps all traces of the progressbar
        upon termination of iteration [default: False]
    file  : `io.TextIOWrapper` or `io.StringIO`, optional
        Specifies where to output the progress messages
        (default: sys.stderr). Uses `file.write(str)` and `file.flush()`
        methods.
    ncols  : int, optional
        The width of the entire output message. If specified,
        dynamically resizes the progressbar to stay within this bound.
        If unspecified, attempts to use environment width. The
        fallback is a meter width of 10 and no limit for the counter and
        statistics. If 0, will not print any meter (only stats).
    mininterval  : float, optional
        Minimum progress display update interval, in seconds [default: 0.1].
    maxinterval  : float, optional
        Maximum progress display update interval, in seconds [default: 10].
        Automatically adjusts `miniters` to correspond to `mininterval`
        after long display update lag. Only works if `dynamic_miniters`
        or monitor thread is enabled.
    miniters  : int, optional
        Minimum progress display update interval, in iterations.
        If 0 and `dynamic_miniters`, will automatically adjust to equal
        `mininterval` (more CPU efficient, good for tight loops).
        If > 0, will skip display of specified number of iterations.
        Tweak this and `mininterval` to get very efficient loops.
        If your progress is erratic with both fast and slow iterations
        (network, skipping items, etc) you should set miniters=1.
    min_time : float, optional
        Longest estimated total time to remain hidden for. The progress
        bar will remain hidden while the estimate for total running time
        is shorter than this value. Once the estimate goes above this,
        then the progress will remain visible even if it drops again.
        [default: 10]
    ascii  : bool, optional
        If unspecified or False, use unicode (smooth blocks) to fill
        the meter. The fallback is to use ASCII characters `1-9 #`.
    disable  : bool, optional
        Whether to disable the entire progressbar wrapper
        [default: None]. If set to None, disable on non-TTY.
    unit  : str, optional
        String that will be used to define the unit of each iteration
        [default: it].
    unit_scale  : bool or int or float, optional
        If 1 or True, the number of iterations will be reduced/scaled
        automatically and a metric prefix following the
        International System of Units standard will be added
        (kilo, mega, etc.) [default: False]. If any other non-zero
        number, will scale `total` and `n`.
    dynamic_ncols  : bool, optional
        If set, constantly alters `ncols` to the environment (allowing
        for window resizes) [default: True].
    smoothing  : float, optional
        Exponential moving average smoothing factor for speed estimates
        (ignored in GUI mode). Ranges from 0 (average speed) to 1
        (current/instantaneous speed) [default: 0.3].
    bar_format  : str, optional
        Specify a custom bar string formatting. May impact performance.
        [default: '{l_bar}{bar}{r_bar}'], where
        l_bar='{desc}: {percentage:3.0f}%|' and
        r_bar='| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, '
          '{rate_fmt}{postfix}]'
        Possible vars: l_bar, bar, r_bar, n, n_fmt, total, total_fmt,
          percentage, rate, rate_fmt, rate_noinv, rate_noinv_fmt,
          rate_inv, rate_inv_fmt, elapsed, remaining, desc, postfix.
        Note that a trailing ": " is automatically removed after {desc}
        if the latter is empty.
    initial  : int, optional
        The initial counter value. Useful when restarting a progress
        bar [default: 0].
    position  : int, optional
        Specify the line offset to print this bar (starting from 0)
        Automatic if unspecified.
        Useful to manage multiple bars at once (eg, from threads).
    postfix  : dict or *, optional
        Specify additional stats to display at the end of the bar.
        Calls `set_postfix(**postfix)` if possible (dict).
    unit_divisor  : float, optional
        [default: 1000], ignored unless `unit_scale` is True.
    gui  : bool, optional
        WARNING: internal parameter - do not use.
        Use tqdm_gui(...) instead. If set, will attempt to use
        matplotlib animations for a graphical output [default: False].
    Returns
    -------
    out  : decorated iterator.
    """
    # Handle keyword arguments for our wrapper
    self.min_time = kwargs.pop("min_time", 10)
    # Handle dials defaults that we want different from tqdm
    # Remove the progress bar when done
    if "leave" not in kwargs:
      kwargs["leave"] = False
    # Dynamically resize to the console width
    if "dynamic_ncols" not in kwargs:
      kwargs["dynamic_ncols"] = True
    # Disable non-tty by default
    if "disable" not in kwargs:
      kwargs["disable"] = None

    # Intercept the fetching of the printing function. Since some
    # methods cache this (e.g. iteration) return something that never
    # changes but instead proxies the call
    self._tqdm_status_printer = self.status_printer
    # A variable to hold the real status printer function
    self._redirect_printer = lambda x: None
    # Called with file to get the real status printing function
    self.status_printer = lambda _: lambda x: self._redirect_printer(x)

    # Set up the tqdm instance same as usual
    super(progress, self).__init__(*args, **kwargs)

    # Skip everything else if we are outright disabled
    if self.disable:
      return

    # Save what the *real* sp is so we can turn on if we need to
    self._tqdm_sp = self._tqdm_status_printer(self.fp)

    # If we are turning on from the start, just do it
    if self.min_time is None or self.min_time <= 0.05:
      self.on()

  def on(self):
    self._redirect_printer = self._tqdm_sp
    self.refresh()

  def __iter__(self):
    # Faster loop if disabled
    if self.disable:
      for value in super(progress, self).__iter__():
        yield value
      return
    # Loop and check if we've turned on each iteration
    for value in super(progress, self).__iter__():
      # Do these calculations inline to save speed on short loops
      # - iterator mode tends to be used for tighter loops
      if self._redirect_printer is not self._tqdm_sp:
        elapsed = self._time() - self.start_t
        # Calculate the rate
        if self.avg_time:
          rate = 1 / self.avg_time
        else:
          rate = self.n / elapsed
        # Use this to calculate the estimate of total loop time
        remaining = ((self.total - self.n) / rate) if rate else 0
        total_time = elapsed + remaining
        # Activate if we're expected to take longer than the min
        if total_time > self.min_time:
          self._redirect_printer = self._tqdm_sp
      yield value

  def update(self, count):
    super(progress, self).update(count)
    self._handle_turn_on()

  def _handle_turn_on(self):
    # We don't want to turn on if we're disabled, obviously
    if self.disable:
      return
    # If we're on already, do nothing further here
    if self._redirect_printer is self._tqdm_sp:
      return
    if self.total_time > self.min_time:
      elapsed = self._time() - self.start_t
      # Hide for two seconds, or the min_time - whatever is less
      if elapsed > min(self.min_time, 2):
        self.on()

  @property
  def rate(self):
    """Calculate rate using the same algorithm as tqdm"""
    if self.avg_time:
      return 1 / self.avg_time
    elapsed = self._time() - self.start_t
    return self.n / elapsed

  @property
  def total_time(self):
    """Return the current estimate of the total time required"""
    elapsed = self._time() - self.start_t
    rate = self.rate
    remaining = ((self.total - self.n) / rate) if rate else 0
    return elapsed + remaining
