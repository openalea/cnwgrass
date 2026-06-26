
.. _respiration_user:

Respiration User Guide
#########################

.. contents::

Introduction
============

Respiration simulates the different respiration fluxes related to the main physiological processes
described in Thornley and Cannell (2000).
The respiration fluxes described are related:
- the growth of shoot and root organs and grains
- phloem loading of sucrose by source organs
- nitrates, ammonium and other ions uptake by roots
- nitrates reduction
- N fixation for legumes only
- residual processes ((cost from protein turn-over, cell ion gradients, futile cycles...)


Inputs of Respiration
========================

Details on each inputs are given in the docstring.


Outputs of Respiration
=========================

Amount of C lost by respiration (µmol C) for each described process.

Details on each outputs are given in the docstring.


Package architecture
=====================

Respiration is a Python package which consists of several Python modules:

* :mod:`openalea.respiration.model`: the state and the equations of the model,