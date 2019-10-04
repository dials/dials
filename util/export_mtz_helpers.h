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

namespace dials { namespace util {

  using cctbx::uctbx::unit_cell;
  using iotbx::mtz::object;
  using scitbx::mat3;

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
                         af::shared<float> s0n) {
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
        p->so[j] = s0n[j];
      }

      p->source[0] = 0;
      p->source[1] = 0;
      p->source[2] = -1;

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

  mat3<double> dials_u_to_mosflm(const mat3<double> dials_U, unit_cell uc) {
    scitbx::af::double6 p = uc.parameters();
    scitbx::af::double6 rp = uc.reciprocal_parameters();
    scitbx::mat3<double> dials_B = uc.fractionalization_matrix().transpose();
    scitbx::mat3<double> dials_UB = dials_U * dials_B;

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

    scitbx::mat3<double> mosflm_U = dials_UB * mosflm_B.inverse();

    return mosflm_U;
  }
}}  // namespace dials::util

#endif
