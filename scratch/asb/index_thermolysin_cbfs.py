from __future__ import division
from libtbx import easy_run
import sys, os

num_procs = 1
def do_work(path):
  root        = os.path.dirname(path)
  basename    = os.path.splitext(os.path.basename(path))[0]
  datablock   = os.path.join(root, basename + "_datablock.json")
  mask        = os.path.join(root, basename + "_mask.p")
  spots       = os.path.join(root, basename + "_strong.pickle")
  indexed     = os.path.join(root, basename + "_indexed.pickle")
  experiments = os.path.join(root, basename + "_experiments.json")
  refined     = os.path.join(root, basename + "_refined.json")
  integrated  = os.path.join(root, basename + "_integrated.pickle")

  print "Preparing to index", basename

  cmd = "dials.import %s output=%s"%(path,datablock)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(datablock): return

  cmd = "dials.generate_mask %s border=1 output.mask=%s"%(datablock, mask)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(mask): return

  cmd = "dials.find_spots input.datablock=%s threshold.xds.sigma_strong=15 min_spot_size=3 lookup.mask=%s output.reflections=%s"%(datablock, mask, spots)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(spots): return

  cmd = "dials.index input.reflections=%s input.datablock=%s method=fft1d beam.fix=all detector.fix=all known_symmetry.unit_cell=93,93,130,90,90,120 known_symmetry.space_group=P6122 n_macro_cycles=5 d_min_final=0.5 refinement.reflections.outlier.algorithm=tukey refinement.reflections.weighting_strategy.override=stills refinement.reflections.weighting_strategy.delpsi_constant=1000000  output.experiments=%s output.reflections=%s"%(spots, datablock, experiments, indexed)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(experiments): return
  if not os.path.exists(indexed): return

  cmd = "dials.refine weighting_strategy.override=stills weighting_strategy.delpsi_constant=1000000 reflections.outlier.algorithm=tukey %s %s  output.experiments_filename=%s"%(experiments, indexed, refined)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(refined): return

  cmd = "dials.integrate outlier.algorithm=null integration.intensity.algorithm=sum %s %s output.reflections=%s"%(refined, indexed, integrated)
  print cmd
  easy_run.call(cmd)
  if not os.path.exists(integrated): return

  print basename, "indexed succesfully"

  results.append((integrated,refined))

if num_procs > 1:
  from multiprocessing import Manager, Process

  mgr = Manager()
  results = mgr.list()
  def worker():
    for item in iter( q.get, None ):
      do_work(item)
      q.task_done()
    q.task_done()

  q = mgr.JoinableQueue()
  procs = []
  for i in range(num_procs):
    procs.append( Process(target=worker) )
    procs[-1].daemon = True
    procs[-1].start()

  for path in sys.argv[1:]:
    q.put(path)

  q.join()

  for p in procs:
    q.put( None )

  q.join()

  for p in procs:
    p.join()
else:
  results = []
  for path in sys.argv[1:]:
    do_work(path)

if len(results) > 0:
  print "Writing phil input to metrology refiner..."

  f = open("all.phil", 'w')
  for (indexed, experiments) in results:
    f.write("input {\n")
    f.write("  experiments = %s\n"%experiments)
    f.write("  reflections = %s\n"%indexed)
    f.write("}\n")
  f.close()
else:
  print "Nothing indexed"

print "Done"
