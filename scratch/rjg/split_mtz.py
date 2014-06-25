from __future__ import division
import iotbx.phil
from libtbx.phil import command_line
from cctbx.array_family import flex
from iotbx.reflection_file_reader import any_reflection_file
from iotbx import merging_statistics
import iotbx.mtz

master_phil_scope = iotbx.phil.parse("""
batch_multiplier = 1000
  .type = int
split_at_batch = None
  .type = int
export_unmerged = False
  .type = bool
""")


def run(args):
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
    args=args, custom_processor="collect_remaining")
  working_phil.show()
  params = working_phil.extract()
  batch_multiplier = params.batch_multiplier
  split_at_batch = params.split_at_batch
  assert split_at_batch is not None

  assert len(args) == 1
  file_name = args[0]
  reader = any_reflection_file(file_name)
  assert reader.file_type() == 'ccp4_mtz'

  as_miller_arrays = reader.as_miller_arrays(merge_equivalents=False)
  mtz_object = reader.file_content()
  intensities = [ma for ma in as_miller_arrays
                 if ma.info().labels == ['I', 'SIGI']][0]
  intensities = intensities.customized_copy(
    indices=mtz_object.extract_original_index_miller_indices()).set_info(
      intensities.info())
  batch_ids = [ma for ma in as_miller_arrays
             if ma.info().labels == ['BATCH']][0]
  batch_ids = batch_ids.customized_copy(
    indices=mtz_object.extract_original_index_miller_indices()).set_info(
      batch_ids.info())
  intensities = intensities.customized_copy(anomalous_flag=True).set_info(
    intensities.info())
  intensities.set_observation_type_xray_intensity()

  run_id = flex.int()
  run_batch_id = flex.int()
  for b in batch_ids.data():
    r_id, b_id = divmod(b, batch_multiplier)
    run_id.append(r_id)
    run_batch_id.append(b_id)

  sel = run_batch_id < split_at_batch

  for negate in (False, True):
    if not negate:
      outfile = "split_1.mtz"
    else:
      outfile = "split_2.mtz"

    intensities_sel = intensities.select(sel, negate=negate).set_info(
      intensities.info())

    result = iotbx.merging_statistics.dataset_statistics(intensities_sel)
    result.show()

    intensities_merged_anom = intensities_sel.merge_equivalents().array()
    intensities_merged = intensities_sel.customized_copy(
      anomalous_flag=False).merge_equivalents().array()

    dataset = intensities_merged_anom.as_mtz_dataset(
      "I", wavelength=intensities.info().wavelength)
    dataset.add_miller_array(
      miller_array=intensities_merged, column_root_label="IMEAN")
    m = dataset.mtz_object()
    m.write(outfile)

  if params.export_unmerged:

    #for negate in (False, True):
      #if not negate:
        #isel = sel.iselection()
        #outfile = "unmerged_1.mtz"
      #else:
        #isel = (~sel).iselection()
        #outfile = "unmerged_2.mtz"

      #m = iotbx.mtz.object()
      #m.set_title(mtz_object.title())
      #m.set_space_group_info(mtz_object.space_group_info())

      #batches = mtz_object.batches()
      #batch_sel = flex.bool(batches.size(), False)
      #for i, b in enumerate(batches):
        #r_id, b_id = divmod(b.num(), batch_multiplier)
        #if (((not negate) and (b_id < split_at_batch)) or
            #(negate and (b_id >= split_at_batch))):
          #o = m.add_batch()
          #o.set_num(b.num())
          #o.set_nbsetid(b.nbsetid())
          #o.set_ncryst(b.ncryst())
          #o.set_time1(b.time1())
          #o.set_time2(b.time2())
          #o.set_title(b.title())
          #o.set_ndet(b.ndet())
          #o.set_theta(b.theta())
          #o.set_lbmflg(b.lbmflg())
          #o.set_alambd(b.alambd())
          #o.set_delamb(b.delamb())
          #o.set_delcor(b.delcor())
          #o.set_divhd(b.divhd())
          #o.set_divvd(b.divvd())
          #o.set_so(b.so())
          #o.set_bbfac(b.bbfac())
          #o.set_bscale(b.bscale())
          #o.set_sdbfac(b.sdbfac())
          #o.set_sdbscale(b.sdbscale())
          #o.set_nbscal(b.nbscal())
          #o.set_cell(b.cell())
          #o.set_lbcell(b.lbcell())
          #o.set_umat(b.umat())
          #o.set_crydat(b.crydat())
          #o.set_lcrflg(b.lcrflg())
          #o.set_datum(b.datum())
          #o.set_detlm(b.detlm())
          #o.set_dx(b.dx())
          #o.set_e1(b.e1())
          #o.set_e2(b.e2())
          #o.set_e3(b.e3())
          #o.set_gonlab(b.gonlab())
          #o.set_jsaxs(b.jsaxs())
          #o.set_ngonax(b.ngonax())
          #o.set_phixyz(b.phixyz())
          #o.set_phistt(b.phistt())
          #o.set_phirange(b.phirange())
          #o.set_phiend(b.phiend())
          #o.set_scanax(b.scanax())
          #o.set_misflg(b.misflg())
          #o.set_jumpax(b.jumpax())
          #o.set_ldtype(b.ldtype())


      #for x in m.crystals():
        #x_ = m.add_crystal(x.name(), x.project_name(), x.unit_cell_parameters())
        #for d in x.datasets():
          #d_ = x_.add_dataset(d.name(), d.wavelength())

          #nref = isel.size()
          #m.adjust_column_array_sizes(nref)
          #m.set_n_reflections(nref)

          #for column in d.columns():
            #d_.add_column(column.label(), column.type()).set_values(
              #column.extract_values())
      #print

      #m.write(outfile)

    for negate in (False, True):
      if negate:
        outfile = "split_unmerged_1.mtz"
      else:
        outfile = "split_unmerged_2.mtz"
      reader = any_reflection_file(file_name)
      mtz_object = reader.file_content()
      if negate:
        isel = (~sel).iselection()
      else:
        isel = sel.iselection()
      mtz_object.delete_reflections(isel)
      mtz_object.write(outfile)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
