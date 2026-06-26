
.. _cnmetabolism_user:

CN-Metabolism User Guide
###################

.. contents::

Introduction
************

    CN-Metabolism is a model simulating the distribution
    of carbon and nitrogen metabolites into wheat plants described at organ level.
    Shoot organs are represented as a set of laminae, sheaths and internodes, while the root system is represented
    as a single compartment. For each organ, CN-Metabolism considers structural material, mobile metabolites (amino acids,
    nitrates and sucrose) and storage metabolites (proteins and fructans). Synthesis, degradation and transport of
    the different metabolites are simulated for each compartment according to the local metabolite concentration
    and the local environment (mainly light, temperature, CO2, soil N). CN-Metabolism represents these flows as a set
    of differential equations which are solved at each time step (hourly) assuming a constant plant architecture.
    CN fluxes related to organ growth are therefore accounted for in a supplementary models (morphogenesis, growth).
    Soil N dynamics and root N uptake can be either calculated by CN-Metabolism or forced by an external model (using external_soil_model=True)


Implementation and architecture
*******************************

Software design
===============

.. _package_architecture:

Package architecture
--------------------

CN-Metabolism is a Python package which consists of several Python modules:

* :mod:`openalea.cnmetabolism.model`: the state and the equations of the model,
* :mod:`openalea.cnmetabolism.parameters`: the parameters of the model,
* :mod:`openalea.cnmetabolism.simulation`: the simulator (front-end) to run the model,
* :mod:`openalea.cnmetabolism.postprocessing`: the post-processing and graph functions,
* :mod:`openalea.cnmetabolism.tools`: tools to help for the validation of the outputs,
* and :mod:`openalea.cnmetabolism.converter`: functions to convert CN-Metabolism inputs/outputs to/from Pandas dataframes.

.. figure:: ./image/architecture.png
   :width: 50%
   :align: center

   CN-Metabolism architecture

Parameters, variables and equations
-----------------------------------

CN-Metabolism is defined at culm scale, the crop being represented as a population of
individual culms. Culms are considered as a set of botanical modules representing
the root system, each photosynthetic organ and the whole grains.

.. figure:: ./image/botanical_description.jpeg
   :width: 25%
   :align: center

   Botanical description of the culm structure of wheat as implemented in the model (Barillot et al., 2016)

Computationally, the :class:`population <cnmetabolism.model.Population>` is described as a composition
of objects, organized in a multiscale tree-like structure:

* a :class:`population <cnmetabolism.model.Population>` contains one or several :class:`plant(s) <cnmetabolism.model.Plant>`,
    * each :class:`plant <cnmetabolism.model.Plant>` contains one or several :class:`axis(es) <cnmetabolism.model.Axis>`,
        * each :class:`axis <cnmetabolism.model.Axis>` contains:
            * one :class:`set of roots <cnmetabolism.model.Roots>`,
            * one :class:`phloem <cnmetabolism.model.Phloem>`,
            * zero or one :class:`set of grains <cnmetabolism.model.Grains>`,
            * and one or several :class:`phytomer(s) <cnmetabolism.model.Phytomer>` ; each :class:`phytomer <cnmetabolism.model.Phytomer>` contains:
                * one :class:`chaff <cnmetabolism.model.Chaff>` ; each :class:`chaff <cnmetabolism.model.Chaff>` contains:
                    * one exposed :class:`chaff element <cnmetabolism.model.ChaffElement>`,
                    * and/or one enclosed :class:`chaff element <cnmetabolism.model.ChaffElement>`,
                * and/or one :class:`peduncle <cnmetabolism.model.Peduncle>` ; each :class:`peduncle <cnmetabolism.model.Peduncle>` contains:
                    * one exposed :class:`peduncle element <cnmetabolism.model.PeduncleElement>`,
                    * and/or one enclosed :class:`peduncle element <cnmetabolism.model.PeduncleElement>`,
                * and/or one :class:`lamina <cnmetabolism.model.Lamina>` ; each :class:`lamina <cnmetabolism.model.Lamina>` contains:
                    * one exposed :class:`lamina element <cnmetabolism.model.LaminaElement>`,
                    * and/or one enclosed :class:`lamina element <cnmetabolism.model.LaminaElement>`,
                * and/or one :class:`internode <cnmetabolism.model.Internode>` ; each :class:`internode <cnmetabolism.model.Internode>` contains:
                    * one exposed :class:`internode element <cnmetabolism.model.InternodeElement>`,
                    * and/or one enclosed :class:`internode element <cnmetabolism.model.InternodeElement>`,
                * and/or one :class:`sheath <cnmetabolism.model.Sheath>` ; each :class:`sheath <cnmetabolism.model.Sheath>` contains:
                    * one exposed :class:`sheath element <cnmetabolism.model.SheathElement>`,
                    * and/or one enclosed :class:`sheath element <cnmetabolism.model.SheathElement>`.

.. figure:: ./image/population.png
   :width: 50%
   :align: center

   The multiscale tree-like structure of a population of plants

The nitrate concentration in soil is stored and computed in objects of type :class:`cnmetabolism.model.Soil`.

These objects include structural, storage and mobile materials, variations in which are
represented by ordinary differential equations driven by the main metabolic activities.
Each object consists of different metabolites and is connected to a common
pool, the phloem, to allow C–N fluxes.

Thus, each class of :mod:`cnmetabolism.model` defines:

* constants to represent the parameters of the model, 
* attributes to store the current state of the model as compartment values, 
* and methods to compute fluxes and derivatives in the system of differential equations. 



The parameters of the model are stored in module :mod:`parameters <cnmetabolism.parameters>`.
Module :mod:`parameters <cnmetabolism.parameters>` follows the same tree-like structure as module :mod:`model <cnmetabolism.model>`.

Front-end
---------

Module :mod:`simulation <cnmetabolism.simulation>` is the front-end of CN-Metabolism, which allows
to :meth:`initialize <cnmetabolism.simulation.Simulation.initialize>` and :meth:`run <cnmetabolism.simulation.Simulation.run>`
a :meth:`simulation <cnmetabolism.simulation.Simulation>`.

At :meth:`initialization step <cnmetabolism.simulation.Simulation.initialize>`, we first check the
consistency of the :attr:`population <cnmetabolism.simulation.Simulation.population>`
and :attr:`soils <cnmetabolism.simulation.Simulation.soils>` given by the user. Then we
set the initial conditions which will be used by the solver.  

When we :meth:`run <cnmetabolism.simulation.Simulation.run>` the model over 1 time step, we first
:meth:`update the initial conditions <cnmetabolism.simulation.Simulation._update_initial_conditions>`.
Then we call the function :func:`odeint <scipy.integrate.odeint>` of the library :mod:`SciPy <scipy>` 
to integrate the system of differential equations over 1 :attr:`time step <cnmetabolism.simulation.Simulation.time_step>`.
The derivatives needed by :func:`odeint <scipy.integrate.odeint>` are computed by 
method :meth:`_calculate_all_derivatives <cnmetabolism.simulation.Simulation._calculate_all_derivatives>`.
If no error occurs and :func:`odeint <scipy.integrate.odeint>` manages to integrate the 
system successfully, then we update the state of the model setting the attributes of :attr:`population` 
and :attr:`soils` (if any) to the compartment values returned by :func:`odeint <scipy.integrate.odeint>`, and
we compute the :meth:`integrative variables of the population <cnmetabolism.model.Population>`.

.. figure:: ./image/run.png
   :width: 25%
   :align: center
   
   A run of the model

Module :mod:`simulation <cnmetabolism.simulation>` also implements :class:`exception handling <cnmetabolism.simulation.SimulationError>`
:mod:`logging <logging>`.

Postprocessing and tools
--------------------------------

After running the simulation over 1 or several time steps, the user can apply :mod:`postprocessing <cnmetabolism.postprocessing>` on
the outputs of the model. These post-processing are defined in module :mod:`postprocessing <cnmetabolism.postprocessing>`, and can be computed
using function :func:`postprocessing.postprocessing <cnmetabolism.postprocessing.postprocessing>`.


Finally, module :mod:`tools <cnmetabolism.tools>` defines functions to:

* set up loggers,

Module :mod:`converter <cnmetabolism.converter>` implements functions to convert
CN-Metabolism internal :attr:`population <cnmetabolism.simulation.Simulation.population>` and
:attr:`soils <cnmetabolism.simulation.Simulation.soils>` :func:`to <cnmetabolism.converter.to_dataframes>`
and :func:`from <cnmetabolism.converter.from_dataframes>` :class:`Pandas dataframes <pandas.DataFrame>`.

Here is an activity diagram of this example:

.. figure:: ./image/example.png
   :width: 15%
   :align: center

   Dataflow of the example

Constraints on use
==================

Consistency of the inputs
------------------------- 

The input :attr:`population <cnmetabolism.simulation.Simulation.population>` given by the user
must the following topological rules: 

* the :class:`population <cnmetabolism.model.Population>` contains at least one :class:`plant <cnmetabolism.model.Plant>`,
* each :class:`plant <cnmetabolism.model.Plant>` contains at least one :class:`axis <cnmetabolism.model.Axis>`,
* each :class:`axis <cnmetabolism.model.Axis>` must have:
    * one :class:`set of roots <cnmetabolism.model.Roots>`,
    * one :class:`phloem <cnmetabolism.model.Phloem>`,
    * zero or one :class:`set of grains <cnmetabolism.model.Grains>`,
    * at least one :class:`phytomer <cnmetabolism.model.Phytomer>`,
* each :class:`phytomer <cnmetabolism.model.Phytomer>` must have at least:
    * one photosynthetic organ, among :class:`chaff <cnmetabolism.model.Chaff>`, :class:`peduncle <cnmetabolism.model.Peduncle>`,
      :class:`lamina <cnmetabolism.model.Lamina>`, :class:`internode <cnmetabolism.model.Internode>`, or :class:`sheath <cnmetabolism.model.Sheath>`,
    * or one :class:`hiddenzone <cnmetabolism.model.HiddenZone>`,
* each :class:`photosynthetic organ <cnmetabolism.model.PhotosyntheticOrgan>` must have one enclosed element and/or one exposed element.
  Elements enclosed and exposed must be of a type derived from class :class:`PhotosyntheticOrganElement <cnmetabolism.model.PhotosyntheticOrganElement>`,
  that is one of :class:`chaff element <cnmetabolism.model.ChaffElement>`, :class:`lamina element <cnmetabolism.model.LaminaElement>`,
  :class:`internode element <cnmetabolism.model.InternodeElement>`, :class:`peduncle element <cnmetabolism.model.PeduncleElement>` or
  :class:`sheath element <cnmetabolism.model.SheathElement>`.
  An element must belong to an organ of the same type (e.g. a :class:`lamina element <cnmetabolism.model.LaminaElement>`
  must belong to a :class:`lamina <cnmetabolism.model.Lamina>`).

If no external soil model is defined, the input :attr:`soils <cnmetabolism.simulation.Simulation.soils>` given by the user must
supply a :class:`soil <cnmetabolism.model.Soil>` for each :class:`axis <cnmetabolism.model.Axis>`.

These rules prevent from inconsistency in the modeled system. There are checked 
automatically at :meth:`initialization step <cnmetabolism.simulation.Simulation.initialize>`.
If the :attr:`population <cnmetabolism.simulation.Simulation.population>` or the :class:`soils <cnmetabolism.model.Soil>`
breaks these rules, then the :class:`simulator <cnmetabolism.simulation.Simulation>` raises an exception
with appropriate error message.

Continuity of the model
-----------------------

To integrate the system of ordinary differential equations (ODE), the function :func:`odeint <scipy.integrate.odeint>` 
takes as first parameters a function which computes the derivatives at t0::

    dy/dt = func(y, t0, ...)

where ``y`` is a vector.

This function is also called RHS (Right Hand Side) function.

In CN-Metabolism, the RHS function is defined by the method :meth:`_calculate_all_derivatives <cnmetabolism.simulation.Simulation._calculate_all_derivatives>`.
   
If the RHS function has a discontinuity, this may lead to integration failure, the raise of an exception, and a premature end of the execution. 
A discontinuity in RHS function could be due to the use of inconsistent parameters or to bug(s) in the equations of the model. 

If you get a warning of type "ODEintWarning: Excess work done on this call (perhaps wrong Dfun type)", 
you can try to increase the value of ``ODEINT_MXSTEP`` defined in the body of the method :meth:`run <cnmetabolism.simulation.Simulation.run>`.
But you should first enable and check the logs to see if you can settle the problem ahead of the integration. 
Sometimes, this warning is just due to a local discontinuity which does not affect the whole result of the simulation.   

.. _inputs_and_outputs:

Inputs and outputs
******************

The inputs and the outputs of the model consist in state variables describing
the state of the population at a given step.   

All state variables are defined in the classes of the module :mod:`cnmetabolism.model`.
At a given step, instances of these classes stored the state parameters and state variables which represent 
the state of the system. Metabolic inputs for each compartment are defined in individual csv files located
in the inputs folder of your project.
  
See module :mod:`cnmetabolism.model` for a documentation on the inputs and outputs of the model.

.. _post_processing:

Post-processing
***************

The functions which compute the post-processing are defined in the module :mod:`cnmetabolism.postprocessing`,
by botanical object.

See module :mod:`cnmetabolism.postprocessing` for a documentation on the post-processing which can be applied
on the outputs of the  of the model.
