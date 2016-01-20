Defining source models
======================

.. warning:: This need to be though through and documented.

While we don't have all the machinery in place to define source models, yet,
we do have a fully functional (though rudimentary) configuration file
that illustrate how we are planning on implementing complex model with
non trivial energy- and time-dependence of the underlying parameters.

The file `srcmodel/config/stationary_point_pl.py
<https://github.com/lucabaldini/ximpol/blob/master/ximpol/srcmodel/config/stationary_point_pl.py>`_
defines a stationary source with a power-law spectrum and energy- and
time-independent polarization degree and angle, but in a way that is suggestive
how we'll go about defining more complex models.

.. literalinclude:: ../ximpol/srcmodel/config/stationary_point_pl.py
   :language: python
   :start-after: pass

In the mid term, the basic data structures that we envisage to define
arbitrary models is something along the lines of

* a ``xSourceComponent`` class, encapsulating all the morphological, spectral
  and polarization properties of a given component;
* a ``xSource`` class, containing a list of ``xSpectralComponents``;
* a ``xSource List`` class, containing a list of sources (i.e., our models).

(Names might change.)
