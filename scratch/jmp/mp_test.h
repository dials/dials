
#include <omptbx/omp_or_stubs.h>
#include <dials/algorithms/integration/fitrs/fit.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformSpec;
  using dials::algorithms::profile_model::gaussian_rs::transform::Forward;


  void test(
      const Beam &beam,
      const Detector &detector,
      const Goniometer &goniometer,
      const Scan &scan,
      double sigma_b,
      double sigma_m,
      std::size_t grid_size,
      af::reflection_table data) {

    DIALS_ASSERT(data.contains("s1"));
    DIALS_ASSERT(data.contains("xyzcal.mm"));
    DIALS_ASSERT(data.contains("shoebox"));

    af::const_ref< vec3<double> > s1 = data["s1"];
    af::const_ref< vec3<double> > xyzmm = data["xyzcal.mm"];
    af::const_ref< Shoebox<> > sbox = data["shoebox"];

    TransformSpec<> spec(
      beam,
      detector,
      goniometer,
      scan,
      sigma_b,
      sigma_m,
      5.0,
      grid_size);

    #pragma omp parallel for
    for (std::size_t i = 0; i < s1.size(); ++i) {

      std::size_t panel = sbox[i].panel;
      int6 bbox = sbox[i].bbox;
      int2 image_size = detector[panel].get_image_size();
      DIALS_ASSERT(sbox[i].is_consistent());
      if (bbox[0] < 0 || bbox[1] > image_size[0] ||
          bbox[2] < 0 || bbox[3] > image_size[1]) {
        continue;
      }

      Forward<> transform(
        spec,
        s1[i],
        xyzmm[i][2],
        sbox[i]);
    }
  }

}}
