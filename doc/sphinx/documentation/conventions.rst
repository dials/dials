Conventions
===========

Coordinate frames
-----------------

The diffractometer equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the vector :math:`\vec{h}` to describe a position in *fractional
reciprocal space* in terms of the reciprocal lattice basis vectors :math:`\vec{a^*},
\vec{b^*}` and :math:`\vec{c^*}`.

.. math::
   :label: miller_index

   \vec{h} = \begin{pmatrix}
   h \\
   k \\
   l \\
   \end{pmatrix} = h \vec{a^*} + k \vec{b^*} + l \vec{c^*}


The special positions at which h, k and l are integer define the *reciprocal
lattice points* for which (hkl) are the *Miller indices*.

The basic diffractometer equation relates a position :math:`\vec{h}` to a
position :math:`\vec{r_\phi}` in *Cartesian reciprocal space*. This space is
defined so that its axes coincide with the axes of the *laboratory frame*. The
distinction is necessary because distances in reciprocal space are measured in
units of :math:`\unicode{x212B}^{-1}`. However, for convenience it is often acceptable to
refer to either Cartesian reciprocal space or the real space laboratory frame as
the "lab frame", when the correct choice is clear by context. The diffractometer
equation is

.. math::
   :label: diffractometer

   \vec{r_\phi} = \mathbf{R} \mathbf{A} \vec{h}

where :math:`\mathbf{R}` is the *goniostat rotation matrix* and
:math:`\mathbf{A}` is the *crystal setting matrix*, while its inverse
:math:`\mathbf{A}^{-1}` is referred to as the *indexing matrix*. The product
:math:`\mathbf{A} \vec{h}` may be written as :math:`\vec{r_0}`, which is a
position in the :math:`\phi`-axis frame, a Cartesian frame that coincides with
the laboratory frame at a rotation angle of :math:`\phi=0`. This makes clear
that the setting matrix does not change during the course of a rotation
experiment (notwithstanding small "misset" rotations --- see
`Orientation matrix`.

For an experiment performed using the rotation method we use here :math:`\phi`
to refer to the angle about the actual axis of rotation, even when this is
effected by a differently labelled axis on the sample positioning equipment
(such as an :math:`\omega` axis of a multi-axis goniometer). Only in code
specifically dealing with sample positioning equipment might we need to redefine
the labels of axes.  Outside of such modules, the rotation angle is :math:`\phi`
and the axis of rotation is :math:`\vec{e}`, which together with the definition
of the laboratory frame determine the rotation matrix :math:`\mathbf{R}`.

Orthogonalisation convention
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following [#Busing1967]_ we may decompose the setting matrix :math:`\mathbf{A}`
into the product of two matrices, conventionally labelled :math:`\mathbf{U}` and
:math:`\mathbf{B}`. We name :math:`\mathbf{U}` the *orientation matrix* and
:math:`\mathbf{B}` the *reciprocal space orthogonalisation matrix*. These names
are in common, but not universal use. In particular, some texts (for example
[#Paciorek1999]_ refer to the product (i.e. our setting matrix) as the
"orientation matrix".

Of these two matrices, :math:`\mathbf{U}` is a pure rotation matrix and is
dependent on the definition of the lab frame, whilst :math:`\mathbf{B}` is not
dependent on this definition. :math:`\mathbf{B}` does depend however on a choice
of orthogonalisation convention, which relates :math:`\vec{h}` to a position in
the *crystal-fixed Cartesian system*. The basis vectors of this orthogonal
Cartesian frame are fixed to the reciprocal lattice *via* this convention.

There are infinitely many ways that :math:`\mathbf{A}` may be decomposed into a
pair :math:`\mathbf{U} \mathbf{B}`. The symbolic expression of
:math:`\mathbf{B}` is simplified when the crystal-fixed Cartesian system is chosen
to be aligned with crystal real or reciprocal space axes. For example,
[#Busing1967]_ use a frame in which the basis vector :math:`\vec{i}` is parallel
to reciprocal lattice vector :math:`\vec{a^*}`, while :math:`\vec{j}` is chosen
to lie in the plane of :math:`\vec{a^*}` and :math:`\vec{b^*}`. Unfortunately,
this convention is then disconnected from the standard *real space*
orthogonalisation convention, usually called the *PDB convention* [#PDB1992]_.
This standard is essentially universal in crystallographic software for the
transformation of fractional crystallographic coordinates to positions in
orthogonal space, with units of :math:`\AA`. In particular, it is the convention
used in the cctbx [#GrosseKunstleve2002]_. The convention states that the
orthogonal coordinate :math:`x` is determined from a fractional coordinate
:math:`u` by

.. math::
   :label: realspaceortho

   \vec{x} = \mathbf{O} \vec{u}

where the matrix :math:`O` is the *real space orthogonalisation matrix*. This
matrix transforms to a crystal-fixed Cartesian frame that is defined such that
its basis vector :math:`\vec{i}` is parallel to the real space lattice vector
:math:`\vec{a}`, while :math:`\vec{j}` lies in the :math:`(\vec{a}, \vec{b})`
plane. The elements of this matrix made explicit in a compact form are

.. math::
   :label: realspaceorthomatrix

   \mathbf{O} =
   \begin{pmatrix}
   a & b\cos{\gamma} &  c\cos{\beta} \\
   0 & b\sin{\gamma} & -c\sin{\beta}\cos{\alpha^*} \\
   0 & 0             &  c\sin{\beta}\sin{\alpha^*} \\
   \end{pmatrix}

It is desirable to specify our *reciprocal space* orthogonalisation convention
in terms of this real space orthogonalisation convention.  [#Giacovazzo2002]_
derives relationships between real and reciprocal space. Of particular interest
from that text we have

.. math::
   :label: realreciprocaltransforms
   :nowrap:

   \begin{eqnarray}
   \vec{x} & = & \mathbf{M}^\mathsf{T} \vec{x}^\prime \nonumber \\
   \vec{x^*} & = & \mathbf{M}^{-1} \vec{x^*}^\prime
   \end{eqnarray}

By analogy, equate :math:`\vec{x^*}^\prime` with :math:`\vec{h}` and
:math:`\mathbf{B}` with :math:`\mathbf{M}^{-1}`. Also equate
:math:`\mathbf{M}^\mathsf{T}` with :math:`\mathbf{O}` and :math:`\vec{x}^\prime`
with :math:`\vec{u}`. We then see that

.. math::
   :label: reciprocalortho

   \mathbf{B} = \left( \mathbf{O}^{-1} \right)^\mathsf{T} = \mathbf{F}^\mathsf{T}

where :math:`\mathbf{F}` is designated the *real space fractionalisation
matrix*.  This is easily obtained in cctbx by a method of a
:samp:`cctbx.uctbx.unit_cell` object.

A symbolic expression for :math:`\mathbf{F}` in terms of the real space unit cell
parameters is given by [#RuppWeb]_ from which we derive :math:`\mathbf{B}` simply:

.. math::
   :label: recipspaceorthomatrix

   \mathbf{B} =
   \begin{pmatrix}
   \frac{1}{a} &
   0 &
   0 \\
   -\frac{\cos{\gamma}}{a\sin{\gamma}} &
   \frac{1}{b\sin{\gamma}} &
   0 \\
   \frac{bc}{V}\left( \frac{\cos{\gamma} \left( \cos{\alpha} - \cos{\beta}\cos{\gamma} \right)}{\sin{\gamma}} - \cos{\beta}\sin{\gamma} \right) &
   -\frac{ac \left( \cos{\alpha} - \cos{\beta}\cos{\gamma} \right)}{V\sin{\gamma}} &
   \frac{ab\sin{\gamma}}{V} \\
   \end{pmatrix}

with :math:`V = abc \sqrt{ 1 - \cos^2{\alpha} - \cos^2{\beta} - \cos^2{\gamma} +
2 \cos{\alpha}\cos{\beta}\cos{\gamma}}`

Orientation matrix
------------------

.. \label{sec:U_matrix}

The matrix :math:`\mathbf{U}` "corrects" for the orthogonalisation convention
implicit in the choice of :math:`\mathbf{B}`. As the crystal-fixed Cartesian
system and the :math:`\phi`-axis frame are both orthonormal, Cartesian frames
with the same scale, it is clear that :math:`\mathbf{U}` must be a pure rotation
matrix. Its elements are clearly dependent on the mutual orientation of these
frames.

It is usual to think of the orientation as a fixed property of the "sweep".  In
practice the orientation is parameterised such that it becomes a function of
time, to account for crystal slippage (the true degree of this is unknown but
expected to be small; Mosflm uses crystal orientation parameters to account for
inadequacies in other aspects of the experimental description). To reconcile
these points, the current orientation may be expanded into a fixed, datum part
and a variable time-dependent part that is parameterised. That gives

.. math::

   \vec{r_\phi} = \mathbf{\Psi}\mathbf{R}\mathbf{U_0}\mathbf{B}\vec{h}

where :math:`\Psi` is the combined rotation matrix for the misset expressed as
three angles, :math:`\psi_x, \psi_y` and :math:`\psi_z` in the laboratory frame.

In Mosflm these angles are converted to their equivalents in the
:math:`\phi-` axis frame, where:

.. math::

   \vec{r_\phi} = \mathbf{R}\mathbf{\Phi}\mathbf{U_0}\mathbf{B}\vec{h}

At this stage it is unclear which set of angles are the best choice for
parameterisation of the crystal orientation.

The laboratory frame
^^^^^^^^^^^^^^^^^^^^

An important design goal of the DIALS project is that all algorithms should be
fully vectorial. By this we mean that it should be possible to change the
reference frame arbitrarily and all calculations should work appropriately in
the new frame.

Nevertheless, it is useful to adopt a particular standard frame of reference for
meaningful comparison of results, communication between components of the
software and for an agreed definition of what the laboratory consists of
(incompatible definitions can be reasonably argued for, such as that it should
be either fixed to the detector, or to the rotation axis and beam).

In the interests of standardisation, we choose to adopt the Image CIF (imgCIF)
reference frame [#Bernstein2006]_, [#Hammersley2006]_.

Summary of coordinate frames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 - :math:`\vec{h}` gives a position in *fractional reciprocal space*, fixed to
   the crystal.
 - :math:`\mathbf{B}\vec{h}` gives that position in the *crystal-fixed Cartesian
   system* (basis aligned to crystal axes by the orthogonalization convention)
 - :math:`\mathbf{UB}\vec{h}` gives the :math:`\phi`-axis frame (rotates with
   the crystal, axes aligned to lab frame at :math:`\phi=0`)
 - :math:`\mathbf{RUB}\vec{h}` gives *Cartesian reciprocal space* (fixed wrt the
   laboratory)
 - The diffraction geometry relates this to the
   direction of the scattering vector :math:`\vec{s}` in the *laboratory frame*
 - Projection along :math:`\vec{s}` impacts an *abstract sensor frame* giving a
   2D position of the reflection position on a sensor.
 - This position is converted to the *pixel position* for the 2D position on an
   image in number of pixels (starting 0,0 at the origin).


.. rubric:: References

.. [#Bernstein2006] `Bernstein, H. J. in Int. Tables Crystallogr. 199–205 (IUCr, 2006). <http://it.iucr.org/Ga/ch3o7v0001/>`_
.. [#Busing1967] Busing, W. R. & Levy, H. A. Angle calculations for 3- and 4-circle X-ray and neutron diffractometers. Acta Crystallogr. 22, 457–464 (1967).
.. [#Giacovazzo2002] Giacovazzo, C. Fundamentals of Crystallography. (Oxofrd University Press, USA, 2002).
.. [#GrosseKunstleve2002] Grosse-Kunstleve, R. W., Sauter, N. K., Moriarty, N. W. & Adams, P. D. The Computational Crystallography Toolbox: crystallographic algorithms in a reusable software framework. J. Appl. Crystallogr. 35, 126–136 (2002).
.. [#Hammersley2006] `Hammersley, A. P., Bernstein, H. J. & Westbrook, D. in Int. Tables Crystallogr. 444–458 (IUCr, 2006). <http://it.iucr.org/Ga/ch4o6v0001/>`_
.. [#Paciorek1999] Paciorek, W. A., Meyer, M. & Chapuis, G. On the geometry of a modern imaging diffractometer. Acta Crystallogr. Sect. A Found. Crystallogr. 55, 543–557 (1999).
.. [#PDB1992] `PDB. Atomic Coordinate and Bibliographic Entry Format Description. Brookhaven Natl. Lab. (1992). <http://www.wwpdb.org/docs/documentation/file-format/PDB_format_1992.pdf>`_
.. [#RuppWeb] `Rupp, B. Coordinate system transformation. <http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm>`_

