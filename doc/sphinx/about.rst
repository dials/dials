+++++
About
+++++

Atomic structures are vital for scientific research.  From understanding
the mechanisms of biology to the design of new materials for industrial applications,
knowing the placement of atoms in a molecule assists in the control of disease progression,
drug design, the development of therapies such as vaccines, building new compounds
for better batteries, and finding new methods of carbon sequestration.

Crystallography is one of the prime technologies used to generate high quality structures,
generated from crystallographic data collected at large facilities like synchrotrons and
X-ray free electron lasers, on large micro electron diffraction instruments, at neutron
sources, or at home in a lab.  At these facilities and instruments, crystals are either
rotated in a beam of photons, neutrons, or electrons, or are exposed one at a time in a
serial experiment, and diffraction data from detectors are collected and stored in image
files. The DIALS software package (Diffraction Integration for Advanced Light Sources) is
used throughout the world for processing crystallographic data in automated pipelines or
through graphical user interfaces.

The DIALS software is developed in a fully open-source, collaborative
environment. The main development teams are based at `Diamond Light Source`_
and `CCP4`_, in the UK, and at `Lawrence Berkeley National Laboratory`_, USA.
However, in the spirit of the open source movement, we welcome
collaboration from anyone who wishes to contribute to the project.

To avoid "reinventing the wheel" as much as possible, the DIALS project builds
on knowledge accumulated over many decades in the field of crystallographic
data processing. We benefit greatly from the altruism of experts who contribute their ideas
and advice, either directly or via their detailed publications on existing
algorithms and packages such as XDS [#XDS]_ and MOSFLM [#MOSFLM]_. At the heart
of the DIALS framework lies a design philosophy of hardware abstraction and a
generalised model of the experiment
that is inspired directly by material published on the seminal workshops on
position sensitive detector software [#Lure]_. Continuing in the spirit of these
workshops we held our own series of meetings, with talks from invited speakers, and code
camps in which specific problems are addressed by intensive effort across the
collaboration. Summaries of these meetings and copies of slides given as
presentations are available :doc:`here </workshops/index>`.

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

DIALS UK
----------

Development of DIALS in the UK is funded by the `Wellcome Trust`_,
`Diamond Light Source`_ and `CCP4`_, and led by `Dr Gwyndaf Evans`_.

DIALS US
----------

Development of DIALS at `Lawrence Berkeley National Laboratory`_, USA is led by
`Dr Aaron Brewster`_ and supported by `National Institutes of Health`_ /
`National Institute of General Medical Sciences`_ grant R24GM154040: *DIALS:
supporting structural biology through open source diffraction processing software*.
Work at LBNL is performed under `Department of Energy`_ contract
DE-AC02-05CH11231.

Acknowledgements
================

We are grateful to all those who have contributed to the development
of DIALS.

.. literalinclude:: ../../AUTHORS
    :language: rst

In addition, we acknowledge guidance and ideas gained through numerous
intellectual discussions with the following:

Alun Ashton,
Gleb Bourenkov,
Gerard Bricogne,
Phil Evans,
Andrew Leslie,
Nigel Moriarty,
Garib Murshudov,
Jim Pflugrath,
Harry Powell,
Jon Schuermann
and
Matthew Webber.

.. [#Lure] `Bricogne, G. (1987). Proceedings of the CCP4 Daresbury Study Weekend, pp. 120-145.`
.. [#XDS] `Kabsch, W. (2010). Acta Cryst. D66, 125-132.`
.. [#MOSFLM] `Leslie, A. G. W. and Powell H. R. (2007), Evolving Methods for Macromolecular Crystallography, 245, 41-51. ISBN 978-1-4020-6314-5.`
.. [#RWGK] `Grosse-Kunstleve, R. W., Sauter, N. K., Moriarty, N. W., & Adams, P. D. (2002). Journal of Applied Crystallography. 35, 126–136.`

.. _`Wellcome Trust`: https://wellcome.ac.uk/
.. _`Boost.Python`: http://www.boost.org/doc/libs/1_59_0/libs/python/doc/index.html
.. _`cctbx`: https://github.com/cctbx
.. _`CCP4`: http://www.ccp4.ac.uk/
.. _`Diamond Light Source`: http://www.diamond.ac.uk/Home.html
.. _`Dr Gwyndaf Evans`: http://www.diamond.ac.uk/Beamlines/Mx/VMXm/Staff/Evans.html
.. _`Dr Aaron Brewster`: https://biosciences.lbl.gov/profiles/aaron-brewster/
.. _`Lawrence Berkeley National Laboratory`: http://www.lbl.gov/
.. _`National Institutes of Health`: http://www.nih.gov/
.. _`National Institute of General Medical Sciences`: http://www.nigms.nih.gov/
.. _`Department of Energy`: http://www.energy.gov/
