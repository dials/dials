``dials.import``: ``geometry.*`` overrides are now also passed to the Format classes as the
images are read, so that a Format class which declares metadata it cannot read from the file
(via dxtbx's ``missing_metadata``) is filled in immediately rather than having to invent dummy
values. If such a Format class is used without the overrides it needs, the import stops with an
error naming exactly which parameters to supply, instead of silently processing the data with
fake geometry.
