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
   :width: 49.5%
.. image:: figures/showcase/casa_model_he_image.png
   :width: 49.5%

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

Below are some outputs of a 250 ks simulated observation of Cas A, based on the
model described above. When the entire source is selected (correponding to the
case where the polarimeter has no imaging capabilities), most of the
polarization averages out.


On the other hand, spatially- and energy-resolved polarimetry could in this
case reveal much of the richness in the original polarization pattern.





The Crab pulsar
---------------


The Crab complex
----------------



