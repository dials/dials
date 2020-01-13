Introduce `dials.rescale_diamond_anvil_cell` to correct the absorption of the incident and diffracted X-ray beam by the diamond anvils in a pressure cell.
Call `dials.rescale_diamond_anvil_cell` on the output of `dials.integrate` and then proceed to use post-integration tools as normal, just as though the sample had been measured in air.
