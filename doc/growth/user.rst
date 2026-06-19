
.. _growthwheat_user:

Growth-Wheat User Guide
#########################

.. contents::

Introduction
============

Growth-Wheat simulates leaf growth in mass according to leaf elongation calculated in Elong-Wheat.
For non-emerged leaves, a empiric and constant relation between leaf length and its structural mass is used.
For emerged leaves, the structural mass growth is computed from their green area (given by a geometrical model such as Adel-Wheat) and
a surfacic mass defined by Elong-Wheat. According to the computed increase in structural mass, Growth-Wheat
calculates the related consumption of sucrose and amino acids. Growth-Wheat also handles the transfer of CN metabolites between the
hiddenzone and the newly emerged part of the leaf if any. When Growth-Wheat is coupled to :mod:`openalea.turgorgrowth` (hydraulics=True),
then leaf mass increase is completely attributed to the lamina during its growth. At the end of leaf elongation, its mass and solutes are distributed
between the lamina and the sheath.
Root growth and the related metabolite consumption are regulated by the local concentration of
sucrose. The maximal rate of root growth also depends on the reproductive status of the plant.
When Growth-Wheat is coupled to :mod:`openalea.turgorgrowth` (hydraulics=True), root growth also depends on the xylem water potential.


Inputs of Growth-Wheat
========================

- Initial structural masses (g), dimensions (m) of each hiddenzone + the green area for laminae, sheaths and internodes
- Initial content of the metabolites calculated by CN-Wheat (µmol)
- 'xylem_water_potential': the water potential of the xylem (MPa)

Details on each inputs are given in the docstring.

Outputs of Growth-Wheat
=========================

Updated structural masses for hiddenzones, shoot organs and roots
Updated CN metabolites contents

Details on each outputs are given in the docstring.


Package architecture
=====================

Growth-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.growthwheat.model`: the state and the equations of the model, two classes available (:class:`GrowthWheatModel <openalea.growthwheat.model.GrowthWheatModel>` and :class:`GrowthWheatModelHydraulics <openalea.growthwheat.model.GrowthWheatModelHydraulics>`)
* :mod:`openalea.growthwheat.parameters`: the parameters of the model,  two classes available (:class:`Parameters <openalea.growthwheat.parameters.Parameters>` and :class:`ParametersHydraulics <openalea.growthwheat.parameters.ParametersHydraulics>`)
* :mod:`openalea.growthwheat.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.growthwheat.converter`: functions to convert Growth-Wheat inputs/outputs to/from Pandas dataframes.