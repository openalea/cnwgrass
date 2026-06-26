
.. _morphogenesis_user:

Morphogenesis User Guide
#########################

.. contents::

Introduction
============

Each growing leaves is represented by a set of compartments:

- a hiddenzone: encompasses growing and mature tissues of the growing leaf that are
  enclosed in the pseudo stem.
- any exposed (emerged) part of the lamina (photosynthetically active)
- any exposed part of the sheath (photosynthetically active).

Model functioning depends on the desired version:

- When leaf status is not taken into account (hydraulics=False)
  Tissue emergence is calculated according to the length of the pseudostem which depends on
  the above sheath and internodes of each growing leaf. Leaf elongation rate is regulated
  by the amino acids and sucrose concentration inside the hiddenzone. The dynamics of leaf elongation
  is coordinated with the emergence of the previous leaf:
    - before leaf n-1 emergence (Phase I), leaf n follows an exponential-like function with a strong metabolic regulation
    - after leaf n-1 emergence (Phase II), leaf n follows a predefined sigmoidal kinetics with low metabolic control

- When leaf status is taken into account (hydraulics=True)
  Leaf elongation rate in Phase I is co-regulated by metabolite concentration in the hidden zone and water potential of the xylem.
  Phase II is almost deactivated, as leaf elongation rate is calculated in :mod:`openalea.hydraulics`:
  In this version, the separation between the lamina and the sheath is only occurring at the end of leaf elongation.


At leaf n emergence, its maximal width and surfacic mass are defined according to the sucrose concentration
of the hiddenzone averaged during to phyllochrons.

Morphogenesis also simulates some functions of the shoot apical meristem : leaf primordia emission and floral transition
Morphogenesis has on option to force the final dimensions of each leaves to those defined in :mod:`openalea.morphogenesis.parameters` module

Inputs of Morphogenesis
========================

- Initial dimensions (m) and structural masses (g) of each hiddenzone, laminae, sheaths and internodes
- Initial content of the metabolites calculated by CN-Metabolism (µmol)
- air and soil temperature (°C) : used to calculate the temperature of the hiddenzones and shoot apical meristem
- the xylem water potential is also needed if Morphogenesis is coupled with a hydraulic model.

Details on each inputs are given in the docstring.

Outputs of Morphogenesis
=========================
Updated leaf dimensions and shoot apical meristem status

Details on each outputs are given in the docstring.


Package architecture
=====================

Morphogenesis is a Python package which consists of several Python modules:

* :mod:`openalea.morphogenesis.model`: the state and the equations of the model, two classes available (:class:`MorphogenesisModel <openalea.morphogenesis.model.MorphogenesisModel>` and :class:`MorphogenesisModelHydraulics <openalea.morphogenesis.model.MorphogenesisModelHydraulics>`)
* :mod:`openalea.morphogenesis.parameters`: the parameters of the model, two classes available (:class:`Parameters <openalea.morphogenesis.parameters.Parameters>` and :class:`ParametersHydraulics <openalea.morphogenesis.parameters.ParametersHydraulics>`)
* :mod:`openalea.morphogenesis.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.morphogenesis.converter`: functions to convert Morphogenesis inputs/outputs to/from Pandas dataframes.