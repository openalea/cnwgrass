.. _integration_user:

Integration User Guide
#########################

.. contents::

Introduction
============

Integration is a subpackage used to couple the models cn-metabolism, gas-exchange, morphogenesis, hydraulics,
growth, respiration, senescence, adel-wheat and caribu. Each model has a specific facade used to read
its inputs from a MTG object and write its outputs after running.
The MTG is provided by adel-wheat model.


Package architecture
=====================

Integration is a Python package which consists of several Python modules:

* :mod:`openalea.cnwgrass.integration.caribu_facade`: the interface between Caribu and the MTG
* :mod:`openalea.cnwgrass.integration.cnmetabolism_facade`: the interface between CN-Metabolism and the MTG
* :mod:`openalea.cnwgrass.integration.morphogenesis_facade`: the interface between Morphogenesis and the MTG
* :mod:`openalea.cnwgrass.integration.gasexchange_facade`: the interface between Gas-Exchange and the MTG
* :mod:`openalea.cnwgrass.integration.build_outputs`: the interface used to extract model results from the MTG
* :mod:`openalea.cnwgrass.integration.runner`: the interface that run the models using the different facades and user inputs
* :mod:`openalea.cnwgrass.integration.growth`: the interface between Growth and the MTG
* :mod:`openalea.cnwgrass.integration.senescence_facade`: the interface between Senescence and the MTG
* :mod:`openalea.cnwgrass.integration.tools`: some tools mainly related to the creation of graphs
* :mod:`openalea.cnwgrass.integration.hydraulics_facade`: the interface between Hydraulics and the MTG