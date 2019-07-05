from __future__ import absolute_import, division, print_function


class _Extension(object):
    """A base class for extension groups.
    This contains a common lookup mechanism and phil scope generator.
    """

    @classmethod
    def load(cls, name):
        """Get the requested extension class by name.

        :param name: The name of the extension
        :returns: The extension class

        """
        for e in cls.extensions():
            if e.name == name:
                return e

    @classmethod
    def phil_scope(cls):
        """Get the phil scope for the interface or extension.

        :returns: The phil scope for the interface or extension

        """
        from libtbx.phil import parse

        if cls == _Extension:
            raise RuntimeError("Extension has no phil parameters")
        doc = "\n".join('"%s"' % d for d in cls.__doc__.strip().splitlines())
        master_scope = parse("%s .help=%s {}" % (cls.name, doc))
        main_scope = master_scope.get_without_substitution(cls.name)
        assert len(main_scope) == 1
        main_scope = main_scope[0]
        if "phil" in cls.__dict__:
            main_scope.adopt_scope(cls.phil())

        def ext_names(extensions):
            names = []
            default_index = -1
            for ext in extensions:
                if "default" in ext.__dict__:
                    default_index = len(names)
                names.append(ext.name)
            if default_index < 0:
                default_index = 0
            if names:
                names[default_index] = "*" + names[default_index]
            return names

        exts = cls.extensions()
        if exts:
            algorithm = parse(
                """
        algorithm = %s
          .help = "The choice of algorithm"
          .type = choice
      """
                % " ".join(ext_names(exts))
            )
            main_scope.adopt_scope(algorithm)
            for ext in exts:
                if hasattr(ext, "phil"):
                    help_str = "\n".join(
                        ['"%s"' % line for line in ext.__doc__.split()]
                    )
                    ext_master_scope = parse("%s .help=%s {}" % (ext.name, help_str))
                    ext_phil_scope = ext_master_scope.get_without_substitution(ext.name)
                    assert len(ext_phil_scope) == 1
                    ext_phil_scope = ext_phil_scope[0]
                    ext_phil_scope.adopt_scope(ext.phil())
                    main_scope.adopt_scope(ext_master_scope)
        return master_scope


class SpotFinderThreshold(_Extension):
    """ Extensions for threshold algorithms to be used in spot finding. """

    scope = "spotfinder"
    name = "threshold"

    @classmethod
    def extensions(cls):
        from dials.extensions.dispersion_spotfinder_threshold_ext import (
            DispersionSpotFinderThresholdExt,
        )
        from dials.extensions.dispersion_extended_spotfinder_threshold_ext import (
            DispersionExtendedSpotFinderThresholdExt,
        )
        from dials.extensions.helen_spotfinder_threshold_ext import (
            HelenSpotFinderThresholdExt,
        )
        from dials.extensions.global_spotfinder_threshold_ext import (
            GlobalSpotFinderThresholdExt,
        )

        return [
            DispersionSpotFinderThresholdExt,
            DispersionExtendedSpotFinderThresholdExt,
            HelenSpotFinderThresholdExt,
            GlobalSpotFinderThresholdExt,
        ]


class ProfileModel(_Extension):
    """
    The interface definition for a profile model.

    """

    scope = "profile"
    name = "profile"

    @classmethod
    def extensions(cls):
        from dials.extensions.gaussian_rs_profile_model_ext import (
            GaussianRSProfileModelExt,
        )

        return [GaussianRSProfileModelExt]


class Centroid(_Extension):
    """ Extensions for centroid algorithms. """

    scope = "integration"
    name = "centroid"

    @classmethod
    def extensions(cls):
        from dials.extensions.simple_centroid_ext import SimpleCentroidExt

        return [SimpleCentroidExt]


class Background(_Extension):
    """ Extensions for background algorithms. """

    scope = "integration"
    name = "background"

    @classmethod
    def extensions(cls):
        from dials.extensions.glm_background_ext import GLMBackgroundExt
        from dials.extensions.gmodel_background_ext import GModelBackgroundExt
        from dials.extensions.simple_background_ext import SimpleBackgroundExt
        from dials.extensions.null_background_ext import NullBackgroundExt
        from dials.extensions.median_background_ext import MedianBackgroundExt

        return [
            GLMBackgroundExt,
            GModelBackgroundExt,
            SimpleBackgroundExt,
            NullBackgroundExt,
            MedianBackgroundExt,
        ]
