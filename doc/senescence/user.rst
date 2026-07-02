
.. _senescence_user:

Senescence User Guide
#########################

.. contents::

Introduction
============

Senescence simulates leaf senescence according to :
 - the ratio between the amount of proteins in the leaf at a given timestep and the maximal protein amount recorded in that leaf.
   When the ratio drops below a threshold, the senescence is triggered. The dynamics of N is provided by CN-Metabolism.
 - a maximal age of the leaf. Leaf age is expressed in time compensated for the effects of temperature.

When the senescence is triggered, a small fraction of the leaf dies, reducing its green area.
A fixed proportion of C-N metabolites are remobilised from the death to living tissues.
Root senescence is assumed to occur only after flor l transition.

Inputs of Senescence
========================

- dimensions, green area, structural masses of each shoot organ
- amounts of C and N metabolites of each shoot organs
- structural mass and C -N metabolites for the root compartment


Outputs of Senescence
=========================

Updated structural masses, green areas, dimensions and C-N amounts


Package architecture
=====================

Senescence is a Python package which consists of several Python modules:

* :mod:`openalea.cnwgrass.senescence.model`: the state and the equations of the model,
* :mod:`openalea.cnwgrass.senescence.parameters`: the parameters of the model,
* :mod:`openalea.cnwgrass.senescence.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.cnwgrass.senescence.converter`: functions to convert Senescence inputs/outputs to/from Pandas dataframes.