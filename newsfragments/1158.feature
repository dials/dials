array_family.flex interface has changed: background and centroid algorithms are
set via public properties. Instead of flex.strategy use functools.partial with
the same signature. as_miller_array() raises KeyError instead of Sorry.
.extract_shoeboxes() lost its verbosity parameter, use log levels instead.
