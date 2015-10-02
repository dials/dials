+++++
About
+++++

The DIALS framework is being developed in a fully open-source, collaborative
environment. In the spirit of the open source movement, we welcome
collaboration from anyone who wishes to contribute to the project. The
framework development is currently a joint effort between developers from
`Diamond Light Source`_ and `CCP4`_, based in the UK, and developers at
`Lawrence Berkeley National Laboratory`_, USA. Our common interests
allow us to distribute effort and combine expertise effectively, greatly
enhancing our joint productivity and allowing us to tackle the sizeable task of
writing a new data processing package within a reasonable time frame.

To avoid "reinventing the wheel" as much as possible, the DIALS project builds
on knowledge accumulated over many decades in the field of data processing for
MX. We benefit greatly from the altruism of experts who contribute their ideas
and advice either directly or via their detailed publications on existing
algorithms and packages. At the heart of the DIALS framework lies a design
philosophy of hardware abstraction and a generalised model of the experiment
that is inspired directly by material published on the seminal workshops on
position sensitive detector software [#Lure]_. Continuing in the spirit of these
workshops we hold regular meetings, with talks from invited speakers, and code
camps in which specific problems are addressed by intensive effort across the
collaboration. Summaries of these meetings and copies of slides given as
presentations are available :doc:`here </links>`. We have begun by reproducing
published spot finding and integration algorithms used by integration packages
such as XDS [#XDS]_ and MOSFLM [#MOSFLM]_. We expect new research to lead to more
advanced or case-specific algorithms in future, which may either be implemented
by us or by other groups who choose to work with the toolkit.

DIALS is written using Python and C++, making heavy use of the `cctbx`_ [#RWGK]_
for core crystallographic calculations and much infrastructure including a
complete build system. Seamless interaction between the C++ and Python
components of this *hybrid system* is enabled by `Boost.Python`_. Python provides
a useful ground for rapid prototyping, after which core algorithms and data
structures may be transferred over to C++ for speed. High level interfaces of
the hybrid system remain in Python, facilitating further development and code
reuse both within DIALS and by third parties.


Development Teams
=================

DIALS East
----------

Development of DIALS in the UK is funded by the `BioStruct-X`_ EU grant,
`Diamond Light Source`_ and `CCP4`_, and led by `Dr Gwyndaf Evans`_.
Developers include Luis Fuentes-Montero, Markus Gerstel, Richard Gildea,
James Parkhurst and Graeme Winter at `Diamond Light Source`_, and David Waterman
at `CCP4`_.

DIALS West
----------

Development of DIALS at `Lawrence Berkeley National Laboratory`_, USA is led by
`Dr Nicholas Sauter`_ and supported by `National Institutes of Health`_ /
`National Institute of General Medical Sciences`_ grant R01-GM095887: *Realizing
new horizons in X-ray crystallography data processing*. Work at LBNL is
performed under `Department of Energy`_ contract DE-AC02-05CH11231.
Developers include Aaron Brewster, Tara Michels-Clark and Iris Young.

Acknowledgements
================

We are grateful to the following people who have contributed to the development
of DIALS either in the form of code, or through numerous intellectual discussions:

Muhamed Amin,
Alun Ashton,
Gleb Bourenkov,
Gerard Bricogne,
Phil Evans,
Nat Echols,
Johan Hattne,
Andrew Leslie,
Nigel Moriarty,
Garib Murshadov,
Takanori Nakane,
Jim Pflugrath,
Harry Powell,
Ian Rees,
Jon Schuermann
and
Matthew Webber.

.. [#Lure] `Bricogne, G. (1987). Proceedings of the CCP4 Daresbury Study Weekend, pp. 120-145.`
.. [#XDS] `Kabsch, W. (2010). Acta Cryst. D66, 125-132.`
.. [#MOSFLM] `Leslie, A. G. W. and Powell H. R. (2007), Evolving Methods for Macromolecular Crystallography, 245, 41-51. ISBN 978-1-4020-6314-5.`
.. [#RWGK] `Grosse-Kunstleve, R. W., Sauter, N. K., Moriarty, N. W., & Adams, P. D. (2002). Journal of Applied Crystallography. 35, 126–136.`

.. _`BioStruct-X`: http://www.biostruct-x.org/
.. _`Boost.Python`: http://www.boost.org/doc/libs/1_59_0/libs/python/doc/index.html
.. _`cctbx`: http://cctbx.sourceforge.net/
.. _`CCP4`: http://www.ccp4.ac.uk/
.. _`Diamond Light Source`: http://www.diamond.ac.uk/Home.html
.. _`Dr Gwyndaf Evans`: http://www.diamond.ac.uk/Beamlines/Mx/VMXm/Staff/Evans.html
.. _`Dr Nicholas Sauter`: http://pbd.lbl.gov/scientists/nicholas-sauter/
.. _`Lawrence Berkeley National Laboratory`: http://www.lbl.gov/
.. _`National Institutes of Health`: http://www.nih.gov/
.. _`National Institute of General Medical Sciences`: http://www.nigms.nih.gov/
.. _`Department of Energy`: http://www.energy.gov/
