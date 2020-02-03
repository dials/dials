dials.anvil_correction
======================

Introduction
------------

.. python_string:: dials.command_line.anvil_correction.help_message

Full parameter definitions
--------------------------

.. phil:: dials.command_line.anvil_correction.phil_scope
   :expert-level: 2
   :attributes-level: 2

Details
-------

The path lengths :math:`l_0` and :math:`l_1` of the incident and diffracted beams through the diamond anvils of the pressure cell are illustrated in the schematic below.
The anvils are assumed to have equal thickness :math:`t` and to be perfectly parallel.
The cell is fixed to the goniometer, so its orientation depends on the goniometer rotation :math:`R`.
When the goniometer is at the zero datum, the anvils' normal is :math:`\mathbf{\hat{n}}`, so in general, the anvils' normal is :math:`R\mathbf{\hat{n}}`.

.. image:: /figures/diamond_anvil_cell.svg
  :width: 600
  :alt: Schematic of a diamond anvil cell diffraction measurement.

Since the magnitude of the incident and diffracted beam vectors :math:`\mathbf{s}_0` and :math:`\mathbf{s}_1` is simply :math:`\left|\mathbf{s}_i\right| = 1/\lambda`, the path lengths :math:`l_0` and :math:`l_1` are

.. math::
  l_i = \frac{t}{\left|\cos{\left(\alpha_i\right)}\right|} = \frac{t}{\left|\lambda\mathbf{s}_i \cdot R\mathbf{\hat{n}}\right|} \,\text.

As a result of absorption in the anvil material, the intensity of each beam is reduced by a factor :math:`\exp{\left(-\mu l_i\right)}`, where :math:`\mu` is the linear absorption coefficient.
Hence each diffraction spot has intensity attenuated by a factor

.. math::
  e^{-\mu\left(l_0 + l_1\right)}\text.

With knowledge of the density :math:`\rho` of the diamond anvils, :math:`\mu` can be calculated from tabulated values of the mass attenuation coefficient for carbon :math:`\left(\mu/\rho\right)_\text{C}`.
The mass attenuation coefficient is taken from data collated by the US National Institute of Standards and Technology [NIST]_.

.. [NIST] Hubbell, J.H. and Seltzer, S.M. (2004), *Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption Coefficients* (version 1.4). [Online] Available: http://physics.nist.gov/xaamdi [2020-01-31]. National Institute of Standards and Technology, Gaithersburg, MD, USA.

After integrating the observed diffraction spot intensities, we can recover an approximation of the intensity that might be expected in the absence of X-ray attenuation by the anvil material.
This is simply achieved by multiplying each of the profile-fitted integrated intensities, the summation integrated intensities and their respective standard deviations by a factor

.. math::
  e^{\left(\mu/\rho\right)_\text{C}\,\rho\,\left(l_0 + l_1\right)}\text.

Note that in the case of the standard deviations, this correction may subtly contradict certain assumptions in the error model of your chosen scaling utility.
The effect is not anticipated to be very significant in most cases and no attempt is made to account for it at this stage.
