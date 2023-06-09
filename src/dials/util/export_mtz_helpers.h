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

  class BatchArrays {
  public:
    BatchArrays(int max_batch_number,
                int dataset_id,
                af::tiny<int, 2> image_range,
                int batch_offset,
                float wavelength,
                float mosaic,
                af::const_ref<float> phi_start,
                af::const_ref<float> phi_range,
                af::const_ref<float, af::c_grid<2> > cell_array,
                af::const_ref<float, af::c_grid<2> > umat_array,
                af::tiny<int, 2> panel_size,
                float panel_distance,
                af::const_ref<float> axis,
                af::const_ref<float> source) {
      int n_batches = image_range[1] - image_range[0] + 1;
      batch_offset += image_range[0] - 1;
      if (max_batch_number > batch_offset) {
        batch_offset = max_batch_number;
      }

      // TODO should resize the arrays properly here and then just add elements by
      // index. What to do about the elements that are set by the Batch struct?

      for (int i_batch = 0; i_batch < n_batches; i_batch++) {
        int batch_num = batch_offset + i_batch + 1;

        // Set the batch ints first
        // ints[0] to ints[2] are set already by the Batch struct
        // ints[3] is batch.number
        _ints.push_back(batch_num);
        // We don't set lbcell (refinement flags for unit cell),
        // because it's probably not used by any program anyway.
        // batch.ints[4] to [9] = left unset
        // Skip jumpax, which is defined as:
        // reciprocal axis closest to principle goniostat axis E1
        // batch.ints[11] = left unset
        // We assume one crystal, one goniostat axis, one detector.
        _ints.push_back(1);  // ncryst
        // batch.ints[13] lcrflg (mosaicity model: 0 = isotropic, 1 = anisotropic)
        _ints.push_back(2);  // ldtype 3D
        _ints.push_back(1);  // jsaxs - goniostat scan axis number
        // batch.ints[16] nbscal (number of batch scales & Bfactors)
        _ints.push_back(1);  // ngonax - number of goniostat axes
        // batch.ints[18] lbmflg (flag for type of beam info: 0 for alambd, delamb; 1
        // also delcor, divhd, divvd)
        _ints.push_back(1);  // ndet
        // ints[20] set by set_dataset_id

        // Set the floats
        // floats[0] to floats[5] is the cell
        for (int j = 0; j < 6; j++) {
          _floats.push_back(cell_array(i_batch, j));
        }
        // floats[6] to floats[14] is the Umat
        for (int j = 0; j < 9; j++) {
          _floats.push_back(umat_array(i_batch, j));
        }
        // floats[15] to floats[20] are missetting angles, here all zero as they are
        // encoded in the U matrix floats[21] is crydat[0]. All other 11 crydat elements
        // are left as zero
        _floats.push_back(mosaic);
        // floats[33] to floats[35] are datum values of goniostat axes, here all zero
        _floats.push_back(phi_start[i_batch]);  // phistt
        _floats.push_back(phi_start[i_batch] + phi_range[i_batch]);
        ;  // phiend
        // floats[38] to floats[40] is scanax, the rotation axis in lab frame
        for (int j = 0; j < 3; j++) {
          _floats.push_back(axis[j]);
        }
        // floats[41] and floats[42] are time1 and time2, left as 0.0
        _floats.push_back(1.0);  // bscale (batch scale)
        // floats[44] is bbfac, left as zero
        // floats[45] is sdbscale, left as zero
        // floats[46] is sdbfac, left as zero
        _floats.push_back(phi_range[i_batch]);  // phirange
        // floats[59] to floats[61] is e1
        for (int j = 0; j < 3; j++) {
          _floats.push_back(axis[j]);
        }
        // floats[62] to floats[64] is e2, all left as zero
        // floats[65] to floats[67] is e3, all left as zero
        // floats[80] to floats[82] is source, the idealised source vector
        // Find element closest to -1 for ideal beam direction
        std::size_t i_min = scitbx::af::min_index(source);
        for (int j = 0; j < 3; j++) {
          if (j == i_min) {
            _floats.push_back(axis[j]);
          } else {
            _floats.push_back(0.0);
          }
        }
        _floats.push_back(-1.f);

        // floats[83] to floats[85] is so the source vector (including tilts)
        for (int j = 0; j < 3; j++) {
          _floats.push_back(source[j]);
        }
        _floats.push_back(wavelength);  // alambd

        // floats[87] is delamb (dispersion (deltalambda / lambda), left as zero
        // floats[88] is delcor (correlated component), left as zero
        // floats[89] is divhd (horizontal beam divergence), left as zero
        // floats[90] is divvd (vertical beam divergence), left as zero

        // floats[111] and floats[112] is dx (xtal to detector distance). Only set the
        // first
        _floats.push_back(panel_distance);

        // floats[113] to floats[120] are detlm (detector limits)
        _floats.push_back(panel_size[0]);  // NX
        _floats.push_back(panel_size[1]);  // NY
      }
    }

    af::shared<float> get_floats() {
      return _floats;
    }

    af::shared<int> get_ints() {
      return _ints;
    }

  private:
    af::shared<float> _floats;
    af::shared<int> _ints;
  };

}}  // namespace dials::util

#endif
