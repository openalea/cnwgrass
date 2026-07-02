
.. _growth_user:

Growth User Guide
#########################

.. contents::

Introduction
============

Growth simulates leaf growth in mass according to leaf elongation calculated in Morphogenesis.
For non-emerged leaves, a empiric and constant relation between leaf length and its structural mass is used.
For emerged leaves, the structural mass growth is computed from their green area (given by a geometrical model such as Adel-Wheat) and
a surfacic mass defined by Morphogenesis. According to the computed increase in structural mass, Growth
calculates the related consumption of sucrose and amino acids. Growth also handles the transfer of CN metabolites between the
hiddenzone and the newly emerged part of the leaf if any. When Growth is coupled to :mod:`openalea.cnwgrass.hydraulics` (hydraulics=True),
then leaf mass increase is completely attributed to the lamina during its growth. At the end of leaf elongation, its mass and solutes are distributed
between the lamina and the sheath.
Root growth and the related metabolite consumption are regulated by the local concentration of
sucrose. The maximal rate of root growth also depends on the reproductive status of the plant.
When Growth is coupled to :mod:`openalea.cnwgrass.hydraulics` (hydraulics=True), root growth also depends on the xylem water potential.


Inputs of Growth
========================

- Initial structural masses (g), dimensions (m) of each hiddenzone + the green area for laminae, sheaths and internodes
- Initial content of the metabolites calculated by CN-Metabolism (µmol)
- 'xylem_water_potential': the water potential of the xylem (MPa)

Details on each inputs are given in the docstring.

Outputs of Growth
=========================

Updated structural masses for hiddenzones, shoot organs and roots
Updated CN metabolites contents

Details on each outputs are given in the docstring.


Package architecture
=====================

Growth is a Python package which consists of several Python modules:

* :mod:`openalea.cnwgrass.growth.model`: the state and the equations of the model, two classes available (:class:`GrowthModel <openalea.cnwgrass.growth.model.GrowthModel>` and :class:`GrowthModelHydraulics <openalea.cnwgrass.growth.model.GrowthModelHydraulics>`)
* :mod:`openalea.cnwgrass.growth.parameters`: the parameters of the model,  two classes available (:class:`Parameters <openalea.cnwgrass.growth.parameters.Parameters>` and :class:`ParametersHydraulics <openalea.cnwgrass.growth.parameters.ParametersHydraulics>`)
* :mod:`openalea.cnwgrass.growth.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.cnwgrass.growth.converter`: functions to convert Growth inputs/outputs to/from Pandas dataframes.