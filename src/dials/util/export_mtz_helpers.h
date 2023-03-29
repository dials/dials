#ifndef DIALS_UTIL_EXPORT_MTZ_HELPERS_H
#define DIALS_UTIL_EXPORT_MTZ_HELPERS_H

#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/constants.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/uctbx.h>
#include <iotbx/mtz/object.h>
#include <scitbx/array_family/misc_functions.h>

#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/mtz.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/symmetry.hpp>

#include <string>

namespace dials { namespace util {

  using cctbx::uctbx::unit_cell;
  using iotbx::mtz::object;
  using scitbx::mat3;
  using scitbx::af::min_index;

  void add_dials_batches(iotbx::mtz::object* mtz,
                         int dataset_id,
                         af::tiny<int, 2> image_range,
                         int batch_offset,
                         float wavelength,
                         float mosaic,
                         af::shared<float> phi_start,
                         af::shared<float> phi_range,
                         af::const_ref<float, af::c_grid<2> > cell_array,
                         af::const_ref<float, af::c_grid<2> > umat_array,
                         af::tiny<int, 2> panel_size,
                         float panel_distance,
                         af::shared<float> axis,
                         af::shared<float> source) {
    // static assertions at the start
#if defined(__DECCXX_VER)
    BOOST_STATIC_ASSERT(sizeof(float) == sizeof(int));
#else
    IOTBX_ASSERT(sizeof(float) == sizeof(int));
#endif

    CMtz::MTZ* ptr = mtz->ptr();
    CMtz::MTZBAT* p = ptr->batch;
    CMtz::MTZBAT* p_tail = p;
    int max_batch_number = 0;
    while (p != 0) {
      max_batch_number = std::max(max_batch_number, p->num);
      p_tail = p;
      p = p->next;
    }
    int i_batch;
    int n_batches = image_range[1] - image_range[0] + 1;
    batch_offset += image_range[0] - 1;
    if (max_batch_number > batch_offset) {
      batch_offset = max_batch_number;
    }
    for (i_batch = 0; i_batch < n_batches; i_batch++) {
      int batch_num = batch_offset + i_batch + 1;

      // original code from add_batch below - probably excessive in this
      // case but seems harmless at this point - moved sizeof() assert out

      boost::scoped_array<float> buf(new float[NBATCHINTEGERS + NBATCHREALS]);
      std::fill_n(buf.get(), NBATCHINTEGERS + NBATCHREALS, static_cast<float>(0));
      IOTBX_ASSERT(CMtz::ccp4_lwbat(ptr, 0, batch_num, buf.get(), "") == 1);
      p = (p_tail == 0 ? ptr->batch : p_tail->next);
      IOTBX_ASSERT(p != 0);
      IOTBX_ASSERT(p->next == 0);
      IOTBX_ASSERT(p->num == batch_num);

      // now we need to write in all the numbers - just use direct data
      // access in the data structure...

      p->nbsetid = dataset_id;
      p->ncryst = 1;
      p->time1 = 0.0;
      p->time2 = 0.0;
      sprintf((char*)&(p->title), "Batch %d", batch_num);
      p->ndet = 1;
      p->theta[0] = 0.0;
      p->theta[1] = 0.0;
      p->lbmflg = 0;
      p->alambd = wavelength;
      p->delamb = 0.0;
      p->delcor = 0.0;
      p->divhd = 0.0;
      p->divvd = 0.0;

      for (int j = 0; j < 3; j++) {
        p->so[j] = source[j];
      }

      p->source[0] = 0;
      p->source[1] = 0;
      p->source[2] = 0;
      // Find element closest to -1 for ideal beam direction
      std::size_t i_min = scitbx::af::min_index(source.const_ref());
      p->source[i_min] = -1;

      p->bbfac = 0.0;
      p->bscale = 1.0;
      p->sdbfac = 0.0;
      p->sdbscale = 0.0;
      p->nbscal = 0;

      // cell constants and flags to indicate what was refined - these should
      // probably be != -1 :-/
      for (int j = 0; j < 6; j++) {
        p->cell[j] = cell_array(i_batch, j);
        p->lbcell[j] = -1;
      }

      for (int j = 0; j < 9; j++) {
        p->umat[j] = umat_array(i_batch, j);
      }

      p->crydat[0] = mosaic;

      for (int j = 1; j <= 11; j++) {
        p->crydat[j] = 0.0;
      }

      p->lcrflg = 0;
      p->datum[0] = 0.0;
      p->datum[1] = 0.0;
      p->datum[2] = 0.0;

      // detector limits
      p->detlm[0][0][0] = 0;
      p->detlm[0][0][1] = panel_size[0];
      p->detlm[0][1][0] = 0;
      p->detlm[0][1][1] = panel_size[1];
      p->detlm[1][0][0] = 0;
      p->detlm[1][0][1] = 0;
      p->detlm[1][1][0] = 0;
      p->detlm[1][1][1] = 0;

      p->dx[0] = panel_distance;
      p->dx[1] = 0;

      // goniometer axes and names, and scan axis number, and num axes, missets
      // should we be using this to unroll the setting matrix etc? at this time
      // assume single-axis goniometer

      p->e1[0] = axis[0];
      p->e1[1] = axis[1];
      p->e1[2] = axis[2];

      p->e2[0] = 0.0;
      p->e2[1] = 0.0;
      p->e2[2] = 0.0;
      p->e3[0] = 0.0;
      p->e3[1] = 0.0;
      p->e3[2] = 0.0;

      strcpy((char*)&(p->gonlab[0]), "AXIS");
      strcpy((char*)&(p->gonlab[1]), "");
      strcpy((char*)&(p->gonlab[2]), "");

      p->jsaxs = 1;
      p->ngonax = 1;

      // missets all 0 - encoded in U matrix
      for (int j = 0; j < 3; j++) {
        p->phixyz[0][j] = 0.0;
        p->phixyz[1][j] = 0.0;
      }

      p->phistt = phi_start[i_batch];
      p->phirange = phi_range[i_batch];
      p->phiend = phi_start[i_batch] + phi_range[i_batch];

      for (int j = 0; j < 3; j++) {
        p->scanax[j] = axis[j];
      }

      p->misflg = 0;
      p->jumpax = 0;
      p->ldtype = 2;

      p_tail = p;
    }
  }

  mat3<double> ub_to_mosflm_u(const mat3<double> UB, unit_cell uc) {
    scitbx::af::double6 p = uc.parameters();
    scitbx::af::double6 rp = uc.reciprocal_parameters();

    double d2r = atan(1.0) / 45.0;

    scitbx::mat3<double> mosflm_B(rp[0],
                                  rp[1] * cos(d2r * rp[5]),
                                  rp[2] * cos(d2r * rp[4]),
                                  0,
                                  rp[1] * sin(d2r * rp[5]),
                                  -rp[2] * sin(d2r * rp[4]) * cos(d2r * p[3]),
                                  0,
                                  0,
                                  1.0 / p[2]);

    scitbx::mat3<double> mosflm_U = UB * mosflm_B.inverse();

    return mosflm_U;
  }

  class GemmiMtzObject {
  public:
    GemmiMtzObject() {}

    void set_title(const char* title) {
      mtz_.title = title;
    }

    void add_history(const char* history) {
      mtz_.history.emplace_back(history);
    }

    void set_space_group_by_name(const char* symbol) {
      mtz_.spacegroup = gemmi::find_spacegroup_by_name(symbol);
    }

    int get_max_batch_number() {
      int max_batch_number = 0;
      for (gemmi::Mtz::Batch b : mtz_.batches)
        max_batch_number = std::max(max_batch_number, b.number);
      return max_batch_number;
    }

    void add_batch(gemmi::Mtz::Batch batch) {
      mtz_.batches.push_back(batch);
    }

    void add_dataset(const char* project_name,
                     const char* crystal_name,
                     const char* dataset_name,
                     af::double6 const& unit_cell_parameters,
                     float wavelength) {
      std::array<double, 6> params;
      for (std::size_t i = 0; i < unit_cell_parameters.size(); i++)
        params[i] = unit_cell_parameters[i];
      gemmi::UnitCell cell(params);
      mtz_.datasets.push_back({++current_data_set_id_,
                               project_name,
                               crystal_name,
                               dataset_name,
                               cell,
                               wavelength});
    }

    void add_column(const char* column_name, const char column_type) {
      mtz_.add_column(column_name, column_type, current_data_set_id_, -1, false);
    }

    void add_column_data(af::const_ref<float> const& values) {
      column_data_.push_back(values);
    }

    void set_n_reflections(int n_reflections) {
      mtz_.nreflections = n_reflections;
    }

  private:
    gemmi::Mtz mtz_;
    int current_data_set_id_ = 0;
    std::vector<af::const_ref<float> > column_data_;
  };

  void add_dials_batches_gemmi(GemmiMtzObject& mtz,
                               int dataset_id,
                               af::tiny<int, 2> image_range,
                               int batch_offset,
                               float wavelength,
                               float mosaic,
                               af::shared<float> phi_start,
                               af::shared<float> phi_range,
                               af::const_ref<float, af::c_grid<2> > cell_array,
                               af::const_ref<float, af::c_grid<2> > umat_array,
                               af::tiny<int, 2> panel_size,
                               float panel_distance,
                               af::shared<float> axis,
                               af::shared<float> source) {
    // static assertions at the start
#if defined(__DECCXX_VER)
    BOOST_STATIC_ASSERT(sizeof(float) == sizeof(int));
#else
    IOTBX_ASSERT(sizeof(float) == sizeof(int));
#endif

    int max_batch_number = mtz.get_max_batch_number();
    int i_batch;
    int n_batches = image_range[1] - image_range[0] + 1;
    batch_offset += image_range[0] - 1;
    if (max_batch_number > batch_offset) {
      batch_offset = max_batch_number;
    }

    for (i_batch = 0; i_batch < n_batches; i_batch++) {
      int batch_num = batch_offset + i_batch + 1;

      // Prepare a batch header
      gemmi::Mtz::Batch batch;
      batch.set_dataset_id(dataset_id);
      batch.title = "Batch " + std::to_string(batch_num);

      // Set the batch ints first
      // ints[0] to ints[2] are set already by the Batch struct
      // ints[3] is batch.number
      batch.number = batch_num;
      // We don't set lbcell (refinement flags for unit cell),
      // because it's probably not used by any program anyway.
      // batch.ints[4] to [9] = left unset
      // Skip jumpax, which is defined as:
      // reciprocal axis closest to principle goniostat axis E1
      // batch.ints[11] = left unset
      // We assume one crystal, one goniostat axis, one detector.
      batch.ints[12] = 1;  // ncryst
      // batch.ints[13] lcrflg (mosaicity model: 0 = isotropic, 1 = anisotropic)
      batch.ints[14] = 2;  // ldtype 3D
      batch.ints[15] = 1;  // jsaxs - goniostat scan axis number
      // batch.ints[16] nbscal (number of batch scales & Bfactors)
      batch.ints[17] = 1;  // ngonax - number of goniostat axes
      // batch.ints[18] lbmflg (flag for type of beam info: 0 for alambd, delamb; 1 also
      // delcor, divhd, divvd)
      batch.ints[19] = 1;  // ndet
      // ints[20] set by set_dataset_id

      // Set the floats
      // floats[0] to floats[5] is the cell
      for (int j = 0; j < 6; j++) {
        batch.floats[j] = cell_array(i_batch, j);
      }
      // floats[6] to floats[14] is the Umat
      for (int j = 0; j < 9; j++) {
        batch.floats[6 + j] = umat_array(i_batch, j);
      }
      // floats[15] to floats[20] are missetting angles, here all zero as they are
      // encoded in the U matrix floats[21] is crydat[0]. All other 11 crydat elements
      // are left as zero
      batch.floats[21] = mosaic;
      // floats[33] to floats[35] are datum values of goniostat axes, here all zero
      batch.floats[36] = phi_start[i_batch];  // phistt
      batch.floats[37] = phi_start[i_batch] + phi_range[i_batch];
      ;  // phiend
      // floats[38] to floats[40] is scanax, the rotation axis in lab frame
      for (int j = 0; j < 3; j++) {
        batch.floats[38 + j] = axis[j];
      }
      // floats[41] and floats[42] are time1 and time2, left as 0.0
      batch.floats[43] = 1.0;  // bscale (batch scale)
      // floats[44] is bbfac, left as zero
      // floats[45] is sdbscale, left as zero
      // floats[46] is sdbfac, left as zero
      batch.floats[47] = phi_range[i_batch];  // phirange
      // floats[59] to floats[61] is e1
      for (int j = 0; j < 3; j++) {
        batch.floats[59 + j] = axis[j];
      }
      // floats[62] to floats[64] is e2, all left as zero
      // floats[65] to floats[67] is e3, all left as zero
      // floats[80] to floats[82] is source, the idealised source vector
      // Find element closest to -1 for ideal beam direction
      std::size_t i_min = scitbx::af::min_index(source.const_ref());
      batch.floats[80 + i_min] = -1.f;

      // floats[83] to floats[85] is so the source vector (including tilts)
      for (int j = 0; j < 3; j++) {
        batch.floats[83 + j] = source[j];
      }
      batch.floats[86] = wavelength;  // alambd

      // floats[87] is delamb (dispersion (deltalambda / lambda), left as zero
      // floats[88] is delcor (correlated component), left as zero
      // floats[89] is divhd (horizontal beam divergence), left as zero
      // floats[90] is divvd (vertical beam divergence), left as zero

      // floats[111] and floats[112] is dx (xtal to detector distance). Only set the
      // first
      batch.floats[111] = panel_distance;

      // floats[113] to floats[120] are detlm (detector limits)
      batch.floats[114] = panel_size[0];  // NX
      batch.floats[116] = panel_size[1];  // NY

      batch.axes.push_back("AXIS");  // gonlab[0]

      mtz.add_batch(batch);
    }
  }

}}  // namespace dials::util

#endif
