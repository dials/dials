Introduce `dials.anvil_correction` to correct the absorption of the incident and diffracted X-ray beam by the diamond anvils in a pressure cell.
Call `dials.anvil_correction` on the output of `dials.integrate` and then proceed to use post-integration tools as normal, just as though the sample had been measured in air.
