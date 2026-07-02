
.. _hydraulics_user:

Hydraulics User Guide
###################

.. contents::

Introduction
************

    Hydraulics is a model simulating water fluxes into wheat plants described at organ level.
    Shoot organs are represented as a set of laminae, sheaths and internodes, while the root system is ignored in the current version.
    Water influx is simulated at organ scale based on water potential differences with a common compartment shared at the axis scale : a "xylem*.
    Leaf elongation follows two distinct phases separated by the emergence of the previous leaf.
    During the initial exponential-like phase, leaf elongation is co-regulated by metabolite concentrations and xylem water potential.
    In the second phase, 3D leaf elongation is simulated using a turgor-driven growth approach,
    whereby metabolic concentrations affect osmotic potential.


Implementation and architecture
*******************************

Software design
===============

.. _package_architecture:

Package architecture
--------------------

Hydraulics is a Python package which consists of several Python modules:

* :mod:`openalea.cnwgrass.hydraulics.model`: the state and the equations of the model,
* :mod:`openalea.cnwgrass.hydraulics.parameters`: the parameters of the model,
* :mod:`openalea.cnwgrass.hydraulics.simulation`: the simulator (front-end) to run the model,
* :mod:`openalea.cnwgrass.hydraulics.postprocessing`: the post-processing and graph functions,
* :mod:`openalea.cnwgrass.hydraulics.tools`: tools to help for the validation of the outputs,
* and :mod:`openalea.cnwgrass.hydraulics.converter`: functions to convert Hydraulics inputs/outputs to/from Pandas dataframes.


Parameters, variables and equations
-----------------------------------

Hydraulics is defined at culm scale, the crop being represented as a population of
individual culms. Culms are considered as a set of botanical modules representing
the root system, each photosynthetic organ and the whole grains.

Computationally, the :class:`population <hydraulics.model.Population>` is described as a composition
of objects, organized in a multiscale tree-like structure:

* a :class:`population <hydraulics.model.Population>` contains one or several :class:`plant(s) <hydraulics.model.Plant>`,
    * each :class:`plant <hydraulics.model.Plant>` contains one or several :class:`axis(es) <hydraulics.model.Axis>`,
        * each :class:`axis <hydraulics.model.Axis>` contains:
            [* one :class:`set of roots <hydraulics.model.Roots>`,] --> NOT YET FUNCTIONAL
            * one :class:`xylem <hydraulics.model.Xylem>`,
            * and one or several :class:`phytomer(s) <hydraulics.model.Phytomer>` ; each :class:`phytomer <hydraulics.model.Phytomer>` contains:
                * and/or one :class:`lamina <hydraulics.model.Lamina>` ; each :class:`lamina <hydraulics.model.Lamina>` contains:
                    * one exposed :class:`lamina element <hydraulics.model.LaminaElement>`,
                    * and/or one enclosed :class:`lamina element <hydraulics.model.LaminaElement>`,
                * and/or one :class:`internode <hydraulics.model.Internode>` ; each :class:`internode <hydraulics.model.Internode>` contains:
                    * one exposed :class:`internode element <hydraulics.model.InternodeElement>`,
                    * and/or one enclosed :class:`internode element <hydraulics.model.InternodeElement>`,
                * and/or one :class:`sheath <hydraulics.model.Sheath>` ; each :class:`sheath <hydraulics.model.Sheath>` contains:
                    * one exposed :class:`sheath element <hydraulics.model.SheathElement>`,
                    * and/or one enclosed :class:`sheath element <hydraulics.model.SheathElement>`.


Soil water content and water potential are stored and computed in objects of type :class:`hydraulics.model.Soil`.

These objects include water related state variables, variations in which are
represented by ordinary differential equations driven by the Lockart and Ortega approach proposed in Coussement et al. (2018).
Each object is connected to a common pool, the xylem, to allow water fluxes at the axis scale.

Thus, each class of :mod:`hydraulics.model` defines:

* constants to represent the parameters of the model, 
* attributes to store the current state of the model as compartment values, 
* and methods to compute fluxes and derivatives in the system of differential equations. 



The parameters of the model are stored in module :mod:`parameters <hydraulics.parameters>`.
Module :mod:`parameters <hydraulics.parameters>` follows the same tree-like structure as module :mod:`model <hydraulics.model>`.

Front-end
---------

Module :mod:`simulation <hydraulics.simulation>` is the front-end of Hydraulics, which allows
to :meth:`initialize <hydraulics.simulation.Simulation.initialize>` and :meth:`run <hydraulics.simulation.Simulation.run>`
a :meth:`simulation <hydraulics.simulation.Simulation>`.

At :meth:`initialization step <hydraulics.simulation.Simulation.initialize>`, we first check the
consistency of the :attr:`population <hydraulics.simulation.Simulation.population>`
and :attr:`soils <hydraulics.simulation.Simulation.soils>` given by the user. Then we
set the initial conditions which will be used by the solver.  

When we :meth:`run <hydraulics.simulation.Simulation.run>` the model over 1 time step, we first
:meth:`update the initial conditions <hydraulics.simulation.Simulation._update_initial_conditions>`.
Then we call the function :func:`odeint <scipy.integrate.odeint>` of the library :mod:`SciPy <scipy>` 
to integrate the system of differential equations over 1 :attr:`time step <hydraulics.simulation.Simulation.time_step>`.
The derivatives needed by :func:`odeint <scipy.integrate.odeint>` are computed by 
method :meth:`_calculate_all_derivatives <hydraulics.simulation.Simulation._calculate_all_derivatives>`.
If no error occurs and :func:`odeint <scipy.integrate.odeint>` manages to integrate the 
system successfully, then we update the state of the model setting the attributes of :attr:`population` 
and :attr:`soils` (if any) to the compartment values returned by :func:`odeint <scipy.integrate.odeint>`, and
we compute the :meth:`integrative variables of the population <hydraulics.model.Population>`.


Module :mod:`simulation <hydraulics.simulation>` also implements :class:`exception handling <hydraulics.simulation.SimulationError>`
:mod:`logging <logging>`.

Postprocessing and tools
--------------------------------

After running the simulation over 1 or several time steps, the user can apply :mod:`postprocessing <hydraulics.postprocessing>` on
the outputs of the model. These post-processing are defined in module :mod:`postprocessing <hydraulics.postprocessing>`, and can be computed
using function :func:`postprocessing.postprocessing <hydraulics.postprocessing.postprocessing>`.


Finally, module :mod:`tools <hydraulics.tools>` defines functions to:

* set up loggers,

Module :mod:`converter <hydraulics.converter>` implements functions to convert
Hydraulics internal :attr:`population <hydraulics.simulation.Simulation.population>` and
:attr:`soils <hydraulics.simulation.Simulation.soils>` :func:`to <hydraulics.converter.to_dataframes>`
and :func:`from <hydraulics.converter.from_dataframes>` :class:`Pandas dataframes <pandas.DataFrame>`.


Constraints on use
==================

Consistency of the inputs
------------------------- 

The input :attr:`population <hydraulics.simulation.Simulation.population>` given by the user
must the following topological rules: 

* the :class:`population <hydraulics.model.Population>` contains at least one :class:`plant <hydraulics.model.Plant>`,
* each :class:`plant <hydraulics.model.Plant>` contains at least one :class:`axis <hydraulics.model.Axis>`,
* each :class:`axis <hydraulics.model.Axis>` must have:
    * one :class:`xylem <hydraulics.model.Xylem>`,
* each :class:`phytomer <hydraulics.model.Phytomer>` must have at least:
    * one photosynthetic organ, among :class:`lamina <hydraulics.model.Lamina>`, :class:`internode <hydraulics.model.Internode>`, or :class:`sheath <hydraulics.model.Sheath>`,
    * or one :class:`hiddenzone <hydraulics.model.HiddenZone>`,
* each :class:`photosynthetic organ <hydraulics.model.PhotosyntheticOrgan>` must have one enclosed element and/or one exposed element.
  Elements enclosed and exposed must be of a type derived from class :class:`PhotosyntheticOrganElement <hydraulics.model.PhotosyntheticOrganElement>`,
  that is one of :class:`lamina element <hydraulics.model.LaminaElement>`,
  :class:`internode element <hydraulics.model.InternodeElement>` or :class:`sheath element <hydraulics.model.SheathElement>`.
  An element must belong to an organ of the same type (e.g. a :class:`lamina element <hydraulics.model.LaminaElement>`
  must belong to a :class:`lamina <hydraulics.model.Lamina>`).

If no external soil model is defined, the input :attr:`soils <hydraulics.simulation.Simulation.soils>` given by the user must
supply a :class:`soil <hydraulics.model.Soil>` for each :class:`axis <hydraulics.model.Axis>`.

These rules prevent from inconsistency in the modeled system. There are checked 
automatically at :meth:`initialization step <hydraulics.simulation.Simulation.initialize>`.
If the :attr:`population <hydraulics.simulation.Simulation.population>` or the :class:`soils <hydraulics.model.Soil>`
breaks these rules, then the :class:`simulator <hydraulics.simulation.Simulation>` raises an exception
with appropriate error message.

Continuity of the model
-----------------------

To integrate the system of ordinary differential equations (ODE), the function :func:`odeint <scipy.integrate.odeint>` 
takes as first parameters a function which computes the derivatives at t0::

    dy/dt = func(y, t0, ...)

where ``y`` is a vector.

This function is also called RHS (Right Hand Side) function.

In Hydraulics, the RHS function is defined by the method :meth:`_calculate_all_derivatives <hydraulics.simulation.Simulation._calculate_all_derivatives>`.
   
If the RHS function has a discontinuity, this may lead to integration failure, the raise of an exception, and a premature end of the execution. 
A discontinuity in RHS function could be due to the use of inconsistent parameters or to bug(s) in the equations of the model. 

If you get a warning of type "ODEintWarning: Excess work done on this call (perhaps wrong Dfun type)", 
you can try to increase the value of ``ODEINT_MXSTEP`` defined in the body of the method :meth:`run <hydraulics.simulation.Simulation.run>`.
But you should first enable and check the logs to see if you can settle the problem ahead of the integration. 
Sometimes, this warning is just due to a local discontinuity which does not affect the whole result of the simulation.   

.. _inputs_and_outputs:

Inputs and outputs
******************

The inputs and the outputs of the model consist in state variables describing
the state of the population at a given step.   

All state variables are defined in the classes of the module :mod:`hydraulics.model`.
At a given step, instances of these classes stored the state parameters and state variables which represent 
the state of the system. Metabolic inputs for each compartment are defined in individual csv files located
in the inputs folder of your project.
  
See module :mod:`hydraulics.model` for a documentation on the inputs and outputs of the model.

.. _post_processing:

Post-processing
***************

The functions which compute the post-processing are defined in the module :mod:`hydraulics.postprocessing`,
by botanical object.

See module :mod:`hydraulics.postprocessing` for a documentation on the post-processing which can be applied
on the outputs of the  of the model.
