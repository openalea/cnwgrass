
.. _senescwheat_user:

Senesc-Wheat User Guide
#########################

.. contents::

Introduction
============

Senesc-Wheat simulates leaf senescence according to :
 - the ratio between the amount of proteins in the leaf at a given timestep and the maximal protein amount recorded in that leaf.
   When the ratio drops below a threshold, the senescence is triggered. The dynamics of N is provided by CN-Wheat.
 - a maximal age of the leaf. Leaf age is expressed in time compensated for the effects of temperature.

When the senescence is triggered, a small fraction of the leaf dies, reducing its green area.
A fixed proportion of C-N metabolites are remobilised from the death to living tissues.
Root senescence is assumed to occur only after flor l transition.

Inputs of Senesc-Wheat
========================

- dimensions, green area, structural masses of each shoot organ
- amounts of C and N metabolites of each shoot organs
- structural mass and C -N metabolites for the root compartment


Outputs of Senesc-Wheat
=========================

Updated structural masses, green areas, dimensions and C-N amounts


Package architecture
=====================

Senesc-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.senescwheat.model`: the state and the equations of the model,
* :mod:`openalea.senescwheat.parameters`: the parameters of the model,
* :mod:`openalea.senescwheat.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.senescwheat.converter`: functions to convert Senesc-Wheat inputs/outputs to/from Pandas dataframes.