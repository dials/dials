DIALS 3.28.0 (2026-02-23)
=========================

No significant changes.


DIALS 3.27.0 (2026-02-23)
=========================

Features
--------

- ``dials.index``: Two new parameters are introduced, ``xy_rmsd_threshold=``, and ``n_indexed_threshold=``, to exclude obvious bad indexing solutions. The default ``xy_rmsd_threshold`` is 6.0 pixels for rotation data and 2.0 pixels for stills, while ``n_indexed_threshold=5`` in either case. (`#3032 <https://github.com/dials/dials/issues/3032>`_)
- Add psana2 support, which allows reading data from LCLS-II file formats (`#3038 <https://github.com/dials/dials/issues/3038>`_)
- New tool: ``dials.tof_integrate`` with 1D and 3D profile fitting methods, for processing time-of-flight data. (`#3039 <https://github.com/dials/dials/issues/3039>`_)
- Add the ``refstat`` algorithm, a utility to check space group symmetry based on Olex2 functionality. An associated command-line program, ``dev.dials.refstat_symmetry_analysis`` is also added, but this is intended to be temporary as the feature will eventually be incorporated into ``dials.symmetry``. (`#3040 <https://github.com/dials/dials/issues/3040>`_)
- ``dials.export``: Add output of hklf2 files. (`#3047 <https://github.com/dials/dials/issues/3047>`_)
- ``dials.index``: Change behaviour of ``joint_index=`` option. Multi-sweep data now requests an explicit ``joint_index=true`` or ``=false`` to specify how to handle multi-sweep data. (`#3051 <https://github.com/dials/dials/issues/3051>`_)
- ``dials.cluster_unit_cell``: Add optional HTML and JSON dendrogram output. (`#3055 <https://github.com/dials/dials/issues/3055>`_)
- ``dials.correlation_matrix``: Use silhouette score for quantifying cluster solution quality, and include dataset labels on clustering plots. (`#3060 <https://github.com/dials/dials/issues/3060>`_)
- ``dials.image_viewer``: Add a new ellipse tool, to aid the calculation of elliptical distortion parameters. (`#3067 <https://github.com/dials/dials/issues/3067>`_)


Bugfixes
--------

- Fix failure to build on Visual Studio 2022, due to excessively large obj files. (`#3019 <https://github.com/dials/dials/issues/3019>`_)
- ``dials.index``: Return best pinkIndexer orientation matrices when ``target_lattices`` > ``max_refine``. Deprecate the min_lattices parameter in favour of target_lattices. (`#3023 <https://github.com/dials/dials/issues/3023>`_)
- ``dials.cosym``: Fix crash if run with ``lattice_group=`` parameter. (`#3048 <https://github.com/dials/dials/issues/3048>`_)
- Reduce memory footprint of the ``Importer`` used by the DIALS ``ArgumentParser`` (`#3072 <https://github.com/dials/dials/issues/3072>`_)
- ``dials.integrate``: Avoid ``FutureWarning`` from deprecated usage of ``functools.partial``. (`#3073 <https://github.com/dials/dials/issues/3073>`_)
- ``dials.integrate``: Include overhead in the calculation of memory required per process, report memory calculations in MB rather than GB, and reduce the default ``max_memory_usage`` parameter to 0.8 rather than 0.9. (`#3074 <https://github.com/dials/dials/issues/3074>`_)
- ``dials.index``: the ``ffbidx`` indexer is no longer installed by default, to reduce the size of DIALS builds. (`#3079 <https://github.com/dials/dials/issues/3079>`_)
- ``dials.rs_mapper``: Extended to work with double precision raw image data. (`#3082 <https://github.com/dials/dials/issues/3082>`_)
- ``dials.image_viewer``: Various improvements to the unit cell and ellipse tools (`#3084 <https://github.com/dials/dials/issues/3084>`_)
- ``dials.rs_mapper``: Print help and exit, when provided empty input. (`#3085 <https://github.com/dials/dials/issues/3085>`_)
- ``dials.image_viewer``: Fix issues with the display of HKL labels on rings made by the unit cell tool. (`#3088 <https://github.com/dials/dials/issues/3088>`_)
- ``dials.split_still_data``: Fix bug in splitting logic for nxs/h5 data in recent versions of dials, where all images were assigned to the first group. (`#3094 <https://github.com/dials/dials/issues/3094>`_)
- Fix tests of the ``refstat`` code on Windows. (`#3095 <https://github.com/dials/dials/issues/3095>`_)
- ``dials.ssx_integrate``: Fix crash if both ``unit_cell.fix=True`` and ``orientation.fix=True`` (`#3096 <https://github.com/dials/dials/issues/3096>`_)
- Fix missing residual calculation in ToF profile1d fitting causing failures on MacOS. (`#3098 <https://github.com/dials/dials/issues/3098>`_)
- Fix autodoc documentation for the reflection table. (`#3104 <https://github.com/dials/dials/issues/3104>`_)
- Fix cosym using a read-only array with pandas 3.0+. (`#3107 <https://github.com/dials/dials/issues/3107>`_)
- ``dials.ssx_index``: Less noisy output to the terminal in Windows. (`#3109 <https://github.com/dials/dials/issues/3109>`_)
- ``dials.correlation_matrix``: Fix bug with OPTICS classification, if single cluster present and no large gradient within allowed max_eps range. (`#3113 <https://github.com/dials/dials/issues/3113>`_)


Improved Documentation
----------------------

- Fix "Show/Hide Log" dropdown buttons in tutorials (`#3093 <https://github.com/dials/dials/issues/3093>`_)


Deprecations and Removals
-------------------------

- Remove IOTA and Prime from installation. These are both deprecated. (`#3030 <https://github.com/dials/dials/issues/3030>`_)


Misc
----

- `#3046 <https://github.com/dials/dials/issues/3046>`_, `#3057 <https://github.com/dials/dials/issues/3057>`_, `#3059 <https://github.com/dials/dials/issues/3059>`_, `#3063 <https://github.com/dials/dials/issues/3063>`_, `#3065 <https://github.com/dials/dials/issues/3065>`_, `#3069 <https://github.com/dials/dials/issues/3069>`_, `#3077 <https://github.com/dials/dials/issues/3077>`_, `#3091 <https://github.com/dials/dials/issues/3091>`_, `#3100 <https://github.com/dials/dials/issues/3100>`_, `#3101 <https://github.com/dials/dials/issues/3101>`_, `#3102 <https://github.com/dials/dials/issues/3102>`_, `#3106 <https://github.com/dials/dials/issues/3106>`_
