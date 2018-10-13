from __future__ import absolute_import, division, print_function

from optparse import SUPPRESS_HELP, OptionParser
import sys

if __name__ == '__main__':
  parser = OptionParser(
      usage="dials.fetch_test_data [-h | --help] [-d destination] [test group]",
      description="This program downloads files used for regression tests. "
                  "These files are not required to run dials or xia2, and are only used in tests.",
  )
  parser.add_option("-?", action="help", help=SUPPRESS_HELP)
  parser.add_option("-d", "--destination", dest="destination", default=None,
                    help="Target directory for files, will be created if not present."
                         "Defaults to <build>/regression_data.")
  parser.add_option("--td", dest="download_threads", type="int", default=8, help="Number of download threads (8)")
  parser.add_option("--tv", dest="verify_threads", type="int", default=8, help="Number of file verification threads (8)")
  parser.add_option("-r", "--retry", dest="retry", type="int", default=3, help="Number of times downloads are retried (3)")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Be more verbose")

  (options, args) = parser.parse_args(sys.argv[1:])

  import dials.util.regression_data
  if not options.destination:
    import libtbx.load_env
    options.destination = libtbx.env.under_build('regression_data')
  print("Downloading regression data into directory %s\n" % options.destination)
  try:
    if not args:
      args = [None]
    for group in args:
      if group:
        print("Downloading files for test group %s" % group)
      success = dials.util.regression_data.fetch_test_data(
          options.destination,
          retry_limit=options.retry, verbose=options.verbose,
          verify_threads=options.verify_threads, download_threads=options.download_threads,
          pre_scan=False, file_group=group,
      )
      if success:
        print("Done.")
      else:
        sys.exit(1)
  except KeyboardInterrupt:
    print("\n\nInterrupted.")
    sys.exit(1)
