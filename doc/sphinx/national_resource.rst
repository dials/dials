++++++++++++++++++++++++++
US DIALS National Resource
++++++++++++++++++++++++++

.. _national_resource:

In the US, DIALS is supported through a R24 National Resource grant (R24GM154040)
from the `National Institutes of Health`_ / `National Institute of General Medical Sciences`_
at `Lawrence Berkeley National Laboratory`_.  The goal of the resource is to support both
users and facilities in using DIALS to analyse their diffraction data. This is part of a
larger international software project that includes support from staff in the US, UK, and
worldwide.

.. container:: twocol

   .. container:: leftside-nih

        .. image:: https://upload.wikimedia.org/wikipedia/commons/f/fb/NIH_2012_logo_arrow.svg
           :alt: NIH NIGMS
           :target: `NIGMS`_
           :height: 70px

   .. container:: leftside-nih

        | National Institute of
        | General Medical Sciences

.. rst-class:: clear-both

Goals
=====

The goals of the Resource are two-fold: 1) support users, instrument scientists,
and beamline scientists in using DIALS to analyze their data and 2) maintain and upgrade DIALS
as a software product.

Community support activities
----------------------------

DIALS has many use cases. Just a few of a variety of possible examples:

-  The individual scientist processing data on their laptop
-  The synchrotron beamline processing thousands of datasets a week an a local cluster
-  The electron diffraction team at a university solving structures using FIB milling single crystal methods
-  The XFEL researcher processing tens of thousands of images a minute on a nation-leading supercomputer

In each case, different support is may be needed.  The individual scientist may need help with parameters for indexing.  The synchroton may need to minimize memory load.  The electron micrscope may require different spotfinder paramters than default.  The XFEL data may need new super-computing submission arguments.  For all of this and more, the Resource can help.  Simply email a staff member or join the community!

Additionally, workshops and training for DIALS will be available and direct support can be arranged by way of site visits by Resource staff.  Please be in touch!

Software activities
-------------------

DIALS is built on cctbx, and in turn many packages depend on DIALS.  The Resource conducts the following activities to maintain and upgrade DIALS:

- Codebase optimization: as processing needs increase and data volumes accelerate, DIALS staff refactors and optimizes the codebase, profiling code, bringing critical algoritms into C++ or onto GPUs, and increasing code stability.
- Codebase maintenance: the Resource tracks modern coding standards and upgrades the code to support them, be they Python, C++, or HDF5 standards.  Additionally as external libraries change and evolve, the Resource continually runs regression tests to monitor for changes that affect code function and performance
- Provide user interfaces: through GUIs and command line programs, DIALS staff build tools that allow data visualization and statistical analysis.

Resource administration
=======================

Staff
-----

- Aaron Brewster: Director, Principal Investigator, representative of DIALS to standards bodies and user communities
- Daniel Paley: senior staff and lead for outreach, maintenance, and optimization
- David Mittan-Moreau: synchrotron and MicroED integration
- Herbert Bernstein: CBF and NeXus support and maintenance

Board
-----

The DIALS advisory board will include crystallographic software developers, beamline facility staff, and representatives of standards bodies to ensure the needs of all aspects of the crystallographic data collection and processing community are being addressed in the DIALS development efforts.  Annual board meetings will include discussions of metrics, improvements, and assessment of progress, and will charge developers to move towards completing specific milestones and deliverables. The board is currently being assembled.

.. _`Lawrence Berkeley National Laboratory`: http://www.lbl.gov/
.. _`National Institutes of Health`: http://www.nih.gov/
.. _`National Institute of General Medical Sciences`: https://www.nigms.nih.gov/
.. _`NIGMS`: https://www.nigms.nih.gov/
