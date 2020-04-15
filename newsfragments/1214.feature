Crystal model now has a new recalculated_unit_cell attribute. This allows it to store
a post-refined unit cell (e.g. from dials.two_theta_refine) in addition to that from
traditional geometry refinement (which was used for prediction). Downstream programs
such as dials.scale and dials.export will now use the recalculated unit cell 
where appropriate.
