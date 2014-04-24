from __future__ import division
import os
import glob
import shutil
from libtbx import easy_run

import iotbx.phil
from libtbx.phil import command_line

master_phil_scope = iotbx.phil.parse("""
template = None
  .type = path
  .multiple = True
find_spots_phil = None
  .type = path
index_phil = None
  .type = path
nproc = 1
  .type = int(value_min=1)
run_xds = False
  .type = bool
run_mosflm = False
  .type = bool
xds {
  include_resolution_range = (40, 0)
    .type = floats(size=2)
  command = *xds xds_par
    .type = choice
}
""")


def run(args):
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()
  if params.find_spots_phil is not None:
    params.find_spots_phil = os.path.abspath(params.find_spots_phil)
    assert os.path.isfile(params.find_spots_phil)
  if params.index_phil is not None:
    params.index_phil = os.path.abspath(params.index_phil)
    assert os.path.isfile(params.index_phil)

  templates = params.template
  print templates

  args = []

  filenames = []

  for t in templates:
    filenames.extend(glob.glob(t))
  from dxtbx.imageset import ImageSetFactory
  imagesets = ImageSetFactory.new(filenames, check_headers=False)
  for i, imageset in enumerate(imagesets):
    print imageset.get_template()
    args.append((imageset.paths(), i+1, params))

  # sort based on the first filename of each imageset
  args.sort(key=lambda x: x[0][0])

  from libtbx import easy_mp

  nproc = params.nproc
  results = easy_mp.parallel_map(
    func=run_once,
    iterable=args,
    processes=nproc,
    method="multiprocessing",
    preserve_order=True,
    asynchronous=True,
    preserve_exception_message=True,
  )

def run_once(args):
  filenames, sweep_id, params = args
  #print filenames

  orig_dir = os.path.abspath(os.curdir)

  sweep_dir = os.path.join(orig_dir, "sweep_%03i" %sweep_id)
  print sweep_dir
  if not os.path.exists(sweep_dir):
    os.mkdir(sweep_dir)
  os.chdir(sweep_dir)

  log = open('%s/sweep_%03i.log' %(sweep_dir, sweep_id), 'wb')
  print >> log, filenames

  cmd = " ".join(["dials.import"] + sorted(filenames))
  print >> log, cmd
  result = easy_run.fully_buffered(command=cmd)
  result.show_stdout(out=log)
  result.show_stderr(out=log)
  args = ["dials.find_spots",
                  "datablock.json",
                  "--nproc=1"
                  ]
  if params.find_spots_phil is not None:
    args.append(params.find_spots_phil)
  cmd = " ".join(args)
  print >> log, cmd
  result = easy_run.fully_buffered(command=cmd)
  result.show_stdout(out=log)
  result.show_stderr(out=log)

  args = ["dials.index",
          "datablock.json",
          "strong.pickle",
          ]

  if params.index_phil is not None:
    args.append(params.index_phil)
  cmd = " ".join(args)
  print >> log, cmd
  result = easy_run.fully_buffered(command=cmd)
  result.show_stdout(out=log)
  result.show_stderr(out=log)

  if params.run_xds:
    g = glob.glob('experiments.json')
    if len(g) == 0:
      return

    cmd = " ".join(["dials.export_xds", "experiments.json"])
    print >> log, cmd
    result = easy_run.fully_buffered(command=cmd)
    result.show_stdout(out=log)
    result.show_stderr(out=log)

    g = glob.glob('xds')
    if len(g) == 0:
      g = glob.glob('xds_[0-9]*')

    g = [os.path.abspath(p) for p in g]

    sweep_dir_full = os.path.abspath('.')

    for xds_dir in g:
      os.chdir(xds_dir)

      os.mkdir("run_1")
      os.chdir("run_1")

      shutil.copyfile("../XDS.INP", "XDS.INP")
      shutil.copyfile("../XPARM.XDS", "XPARM.XDS")

      no_scale = True

      # only refine crystal parameters since we probably know the detector,
      # beam and rotation axis parameters much more accurately from the
      # reference dataset
      with open("XDS.INP", "ab") as f:
        print >> f, "REFINE(INTEGRATE)= ORIENTATION CELL"
        print >> f, "REFINE(CORRECT)= ORIENTATION CELL"
        print >> f, "INCLUDE_RESOLUTION_RANGE= %.1f %.1f" %tuple(
          params.xds.include_resolution_range)

        #if no_scale:
          #print >> f, "MINIMUM_I/SIGMA=50"
          #print >> f, "CORRECTIONS="
          #print >> f, "NBATCH=1"

      result = easy_run.fully_buffered(command=params.xds.command)
      with open("xds.log", "wb") as xds_log:
        result.show_stdout(out=xds_log)
        result.show_stderr(out=xds_log)

      os.chdir("../")

      if os.path.exists("run_1/GXPARM.XDS"):

        os.mkdir("run_2")
        shutil.copyfile("XDS.INP", "run_2/XDS.INP")
        shutil.copyfile("run_1/GXPARM.XDS", "run_2/XPARM.XDS")
        os.chdir("run_2")

        # don't refine anything more the second time
        with open("XDS.INP", "ab") as f:
          print >> f, "REFINE(INTEGRATE)="
          print >> f, "REFINE(CORRECT)="
          print >> f, "INCLUDE_RESOLUTION_RANGE= %.1f %.1f" %tuple(
          params.xds.include_resolution_range)

          if no_scale:
            print >> f, "MINIMUM_I/SIGMA=50"
            print >> f, "CORRECTIONS="
            print >> f, "NBATCH=1"

        result = easy_run.fully_buffered(command=params.xds.command)
        with open("xds.log", "wb") as xds_log:
          result.show_stdout(out=xds_log)
          result.show_stderr(out=xds_log)

  elif params.run_mosflm:
    g = glob.glob('experiments.json')
    if len(g) == 0:
      return

    cmd = " ".join(["dials.export_mosflm", "experiments.json"])
    print >> log, cmd
    result = easy_run.fully_buffered(command=cmd)
    result.show_stdout(out=log)
    result.show_stderr(out=log)

    g = glob.glob('mosflm')
    if len(g) == 0:
      g = glob.glob('mosflm_[0-9]*')

    g = [os.path.abspath(p) for p in g]

    sweep_dir_full = os.path.abspath('.')

    for mosflm_dir in g:
      os.chdir(mosflm_dir)

      with open("mosflm.in", "ab") as mosflm_in:
        print >> mosflm_in, """\
MOSAIC 0.2
refinement residual 15.0
refinement include partials

postref fix all
postref nosegment
process 1 1
go
  """ %(i, i)

      os.chdir("mosflm")
      cmd = "ipmosflm < mosflm.in"
      print >> log, cmd
      result = easy_run.fully_buffered(cmd)
      result.show_stdout(out=log)
      result.show_stderr(out=log)

      os.chdir("../")

  log.close()

  os.chdir(orig_dir)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
