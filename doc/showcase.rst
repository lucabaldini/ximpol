.. _showcase:

Showcase
========

Here are a few source models and simulation outputs illustrating some of the
basic ximpol capabilities.


.. _showcase_casa

Cas A
-----

Full source model definition in `ximpol/config/casa.py
<https://github.com/lucabaldini/ximpol/blob/master/ximpol/config/casa.py>`_.

**Input model**

The spectral model is taken (by hand) from Figure 5 of E.A. Helder and J. Vink,
*"Characterizing the non-thermal emission of Cas A"*, Astrophys. J. 686 (2008)
1094--1102 (`arxiv link <http://arxiv.org/abs/0806.3748>`_). The spectrum of
Cas A is a complex superposition of thermal and non thermal emission, and
for our purposes, we call *thermal* anything that is making up for the lines
and *non-thermal* all the rest, as illustrated in the figure below.

.. image:: figures/showcase/casa_model_spectrum.png
   :width: 75%
   :align: center

The morphology of the source is energy-dependent in a non trivial way.
We start from two Chandra images (in the 1.5--3 keV and 4--6 keV energy ranges,
respectively) and associate the former to the thermal spectral component
and the latter to the non-thermal one (note the absence of spectral lines
between 4 and 6 keV).

.. image:: figures/showcase/casa_model_le_image.png
   :width: 49.6%
.. image:: figures/showcase/casa_model_he_image.png
   :width: 49.6%

For the polarization, we assume that the thermal component is unpolarized,
while for the non-thermal component we use a simple geometrical, radially
symmetric model (loosely inspired from radio observations) where the
polarization angle is tangential and the polarization degree is zero at the
center of the source and increases toward the edges (see figure below).
 
.. image:: figures/showcase/casa_model_he_polmap.png
   :width: 75%
   :align: center
           
Our total model of the region of interest is therefore the superposition of
two indipendent components, with different spectral, morphological and
polarimetric properties. Crude as it is, it's a good benchmark for the
observation simulator.


**Simulation output**

Below is a binned count map of a 250 ks simulated XIPE observation of Cas A,
based on the model described above.

.. image:: figures/showcase/casa_cmap.png
   :width: 75%
   :align: center

When the entire source is analyzed at once, most of the polarization averages
out and even in the high-energy band, where the emission is predominantly
non-thermal, the residual polarization degree resulting from the averaging
of the different emission regions is of the order of 5%.

.. image:: figures/showcase/casa_mod_le.png
   :width: 49.6%
.. image:: figures/showcase/casa_mod_he.png
   :width: 49.6%

On the other hand, spatially- and energy-resolved polarimetry would in this
case reveal much of the richness in the original polarization pattern.
Below is an example of the azimuthal distributions in the two energy bands
for the circular region of interest indicated by the white circle in the left
plot. (The green and blue lines in the ROI indicate the reconstructed
polarization angle.) The comparison with the previous, spatially averaged
distributions is striking.

.. image:: figures/showcase/casa_reg0009_mcube.png
   :width: 100%

By mapping the entire field ov view with suitable regions of interest we can
in fact (at least qualitatively) recover the input polarization pattern,
as shown in the figre below. (Note that at the center of the image the
polarization is close to zero and the arrows have little meaning.)

.. image:: figures/showcase/casa_reg_all.png
   :width: 75%
   :align: center
           

The Crab pulsar
---------------


The Crab complex
----------------



