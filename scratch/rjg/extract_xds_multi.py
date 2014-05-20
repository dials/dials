from __future__ import division
import os
import glob

from libtbx import easy_run

import iotbx.phil
from libtbx.phil import command_line

master_phil_scope = iotbx.phil.parse("""
resolve_indexing_ambiguity = False
  .type = bool
  .help = "Attempt to resolve an indexing ambiguity using the "
          "Brehm-Diederichs algorithm."
overlaps {
  find_overlaps = False
    .type = bool
  max_overlap_fraction = 0.0
    .type = float(value_min=0.0)
  max_overlap_pixels = 0
    .type = int(value_min=0)
}
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
nproc = 1
  .type = int(value_min=1)
aimless {
  command = None
    .type = str
    .multiple = True
  executable = aimless
    .type = path
}
""")

def run(args):
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
    args=args, custom_processor="collect_remaining")
  working_phil.show()
  params = working_phil.extract()

  files = args

  from cctbx import crystal
  from iotbx.reflection_file_reader import any_reflection_file

  file_name_dict = {}

  wedge_id = -1
  wedge_number = -1

  wedge_number_to_wedge_id = {}

  assert params.space_group is not None
  assert params.unit_cell is not None
  space_group = params.space_group.group()
  unit_cell = params.unit_cell
  crystal_symmetry = crystal.symmetry(
    unit_cell=unit_cell, space_group=space_group)

  for file_name in files:
    file_name = os.path.abspath(file_name)
    print file_name
    wedge_number_ = None
    for s in file_name.split(os.path.sep):
      if s.startswith('sweep_'):
        wedge_number_ = int(os.path.splitext(s)[0][-3:])
        print "wedge_number:", wedge_number_
        break
    if wedge_number_ is not None:
      wedge_number = wedge_number_
    else:
      wedge_number += 1
    lattice_id = 1
    for s in file_name.split(os.path.sep):
      if s.startswith('lattice_'):
        lattice_id = int(os.path.splitext(s)[0].split('_')[-1])
        print "lattice_id:", lattice_id
        break
    wedge_id += 1
    print "wedge_id: %i, wedge_number: %i, lattice_id: %i" %(
      wedge_id, wedge_number, lattice_id)
    wedge_number_to_wedge_id.setdefault(wedge_number, [])
    wedge_number_to_wedge_id[wedge_number].append(wedge_id)

    #if not intensities.crystal_symmetry().is_similar_symmetry(
      #crystal_symmetry, relative_length_tolerance=0.1):
      #continue

    file_name_dict[wedge_id] = file_name

  if params.overlaps.find_overlaps:
    # figure out the overlapping reflections and save the miller indices
    # for later on
    reject_hkl = {}
    for wedge_n, wedge_ids in wedge_number_to_wedge_id.iteritems():
      print "Wedge", wedge_n
      if len(wedge_ids) > 1:
        for wedge_id in wedge_ids:
          args = ["dials.import_xds",
                  os.path.split(file_name_dict[wedge_id])[0],
                  "--output='experiments_%i.json'" %wedge_id]
          cmd = " ".join(args)
          print cmd
          result = easy_run.fully_buffered(cmd).raise_if_errors()

          args = ["dials.import_xds",
                  file_name_dict[wedge_id],
                  "experiments_%i.json" %wedge_id,
                  "--input=reflections",
                  "--output='integrate_hkl_%i.pickle'" %wedge_id]
          cmd = " ".join(args)
          print cmd
          result = easy_run.fully_buffered(cmd).raise_if_errors()

        from dials.command_line import find_overlaps
        args = ['experiments_%i.json' %wedge_id for wedge_id in wedge_ids]
        args.extend(['integrate_hkl_%i.pickle' %wedge_id for wedge_id in wedge_ids])
        args.append("nproc=%s" %params.nproc)
        args.append("max_overlap_fraction=%f" %params.overlaps.max_overlap_fraction)
        args.append("max_overlap_pixels=%f" %params.overlaps.max_overlap_pixels)
        args.append("save_overlaps=False")
        overlaps = find_overlaps.run(args)
        miller_indices = overlaps.overlapping_reflections['miller_index']
        overlapping = [
          miller_indices.select(
            overlaps.overlapping_reflections['id'] == i_lattice)
          for i_lattice in range(len(wedge_ids))]
        for wedge_id, overlaps in zip(wedge_ids, overlapping):
          reject_hkl[wedge_id] = overlaps

  for wedge_n, wedge_ids in wedge_number_to_wedge_id.iteritems():
    for wedge in wedge_ids:
      cmd = """\
pointless -copy xdsin %s hklout integrate_hkl_%03.f.mtz << EOF
SPACEGROUP %s
EOF
""" %(file_name_dict[wedge], wedge, space_group.type().lookup_symbol())
      log = open('pointless_%03.f.log' %wedge, 'wb')
      print >> log, cmd
      result = easy_run.fully_buffered(command=cmd)
      result.show_stdout(out=log)
      result.show_stderr(out=log)

      if params.overlaps.find_overlaps:
        from cctbx import miller
        from iotbx import mtz
        m = mtz.object(file_name="integrate_hkl_%03.f.mtz" %wedge)
        orig_indices = m.extract_original_index_miller_indices()
        overlaps = reject_hkl.get(wedge)
        if overlaps is not None and len(overlaps) > 0:
          matches = miller.match_multi_indices(overlaps, orig_indices)
          before = m.n_reflections()
          print "before: %i reflections" %m.n_reflections()
          for i_ref in sorted(matches.pair_selection(1).iselection(), reverse=True):
            m.delete_reflection(i_ref)
          after = m.n_reflections()
          print "after: %i reflections" %m.n_reflections()
          m.add_history("Removed %i overlapping reflections" %len(overlaps))
          m.write("integrate_hkl_%03.f.mtz" %wedge)

  g = glob.glob("integrate_hkl_*.mtz")

  if params.resolve_indexing_ambiguity:
    from cctbx.command_line import brehm_diederichs
    args = g
    args.append("asymmetric=1")
    args.append("save_plot=True")
    args.append("show_plot=False")
    brehm_diederichs.run(args)
  g = glob.glob("integrate_hkl_*_reindexed.mtz")

  for file_name in g:
    wedge_number = int(os.path.splitext(
      os.path.basename(file_name))[0].replace('_reindexed', '')[-3:])
    #print wedge_number, wedge_number
    result = any_reflection_file(file_name)
    mtz_object = result.file_content()
    #if not mtz_object.crystals()[0].crystal_symmetry().is_similar_symmetry(
      #crystal_symmetry, relative_length_tolerance=0.1):
      #continue
    for batch in mtz_object.batches():
      batch.set_num(batch.num() + 1000 * wedge_number)
    batches = mtz_object.get_column('BATCH')
    batches.set_values(batches.extract_values() + 1000*wedge_number)
    mtz_object.write("rebatch-%i.mtz" %(wedge_number))

  g = glob.glob("rebatch-*.mtz")

  cmd = """\
pointless -copy hklin %s hklout pointless.mtz << EOF
ALLOW OUTOFSEQUENCEFILES
TOLERANCE 4
SPACEGROUP %s
EOF
""" %(" ".join(g), space_group.type().lookup_symbol())

  log = open('pointless_all.log', 'wb')
  print >> log, cmd
  result = easy_run.fully_buffered(command=cmd)
  result.show_stdout(out=log)
  result.show_stderr(out=log)

  cmd = """\
aimless pointless.mtz << EOF
OUTPUT UNMERGED TOGETHER
%s
EOF
""" %("\n".join(params.aimless.command))

  log = open('aimless.log', 'wb')
  print >> log, cmd
  result = easy_run.fully_buffered(command=cmd)
  result.show_stdout(out=log)
  result.show_stderr(out=log)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
