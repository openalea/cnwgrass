
.. _elongwheat:

.. module:: openalea.elongwheat

Elong-Wheat documentation
############################

Module description
==================

.. topic:: Overview

    Elong-Wheat is a model of leaf elongation. Individual leaves are represented by a hiddenzone
    compartment (encompassing all growing and mature tissues of the growing leaves located inside the pseudostem) and
    any exposed tissues of the growing leaves.

    Two main approaches are available in Elongwheat depending on the coupling to :mod:`openalea.turgorgrowth`:

    * The rates of leaf elongation are totally regulated by the C-N concentration in the hiddenzone and coordination rules between successive leaves.
      Leaf water status is not taking into account.

    * The rates of leaf elongation is co-regulated by (i) the leaf water status and (ii) the C-N concentration in the hiddenzone and coordination rules between successive leaves.

    

Documentation
=============

.. toctree::
    :maxdepth: 2

    User Guide<user.rst>   
    Reference Guide<ref.rst>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |elongwheat| replace:: :mod:`elongwheat`

