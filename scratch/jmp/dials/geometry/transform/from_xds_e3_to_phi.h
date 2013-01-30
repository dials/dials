
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_XDS_E3_TO_PHI_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_XDS_E3_TO_PHI_H

namespace dials { namespace geometry { namespace transform {

/** A class to transform from XDS e3 coord to the rotation angle, phi */
class FromXdsE3ToPhi {

public:

    /** Default constructor */
    FromXdsE3ToPhi()
        : zeta_(0.0),
          phi_(0.0) {}
          
    /** 
     * Initialise the class 
     * @param zeta The xds zeta parameters
     * @param phi The rotation angle
     */
    FromXdsE3ToPhi(double zeta, double phi)
        : zeta_(zeta),
          phi_(phi) {}
    
    /**
     * Apply the transform
     * @param e3 The XDS e3 coordinate
     * @returns The rotation angle phi'
     */
    double apply(double e3) const {
        return e3 / zeta_ + phi_;
    }

private:

    double zeta_;
    double phi_;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_XDS_E3_TO_PHI_H
