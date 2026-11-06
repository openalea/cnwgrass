.. _metafspm_gas_exchange:

Gas exchange as a metafspm component
####################################

.. contents::
   :local:

``openalea.cnwgrass.metafspm.gas_exchange`` is CNW-Grass's :ref:`gas exchange <gasexchange>` as a
`metafspm <https://github.com/openalea/metafspm>`_ component, ``GasExchange``: the first of CNW-Grass's submodules
ported (the others are not, yet). It computes, for each shoot element of a plant model's DataStructure (Adel's
``StemElement`` and ``LeafElement``), its photosynthesis (Farquhar), its stomatal conductance (BWB, Leuning, Tuzet or
both of the last), its surface temperature and its transpiration, at every step of a metafspm ``Scene``, with the
legacy's equations and numerics.

Using it
========

A plant model adds it to its components; the elements' height comes from Adel's ``z1``, ``z2`` through a translator
formula (``adel_links()``); ``width``, ``diameter`` and ``green_area`` are Adel's, shared by their names; ``PARa``
is the light's (Caribu's ``CaribuLightModel``, an environment model of the scene)::

    from openalea.adel.development import ScenarioShoot
    from openalea.cnwgrass.metafspm.gas_exchange import GasExchange, adel_links

    class Wheat(ScenarioShoot):
        translator = adel_links()

        def __init__(self, data_structure, time_step):
            super().__init__(data_structure, time_step)
            self.gas_exchange = GasExchange(data_structure=data_structure, time_step=time_step)
            self.couple(self.development, self.architecture, self.gas_exchange, translator=self.translator)

        def run(self):
            super().run()
            self.gas_exchange()

    scene = Scene(planting([Planting()], Wheat, {}), environment=[CaribuLightModel], time_step=3600.,
                  forcings=meteo)            # PARi, air_temperature, relative_humidity, ambient_CO2, wind_speed

- **Outputs** per element: ``Ag``, ``An``, ``Rd`` (µmol m-2 s-1), ``gs`` (mol m-2 s-1), ``Ci`` (µmol mol-1),
  ``Ts`` (°C), ``Tr`` (mmol m-2 s-1), the area it exchanges over ``exchanging_area`` (m2) and ``transpiration``
  (``Tr`` · ``exchanging_area``, mmol s-1, extensive: its sum per Organ or per plant is the organ's or the plant's);
  per plant, ``canopy_height``; the capacities at 25 °C (``Vc_max25``, ...). The fluxes per unit area are weighted by
  ``exchanging_area``.
- **An element crossing the soil surface** (a shoot whose crown is below the soil, as Caribu's light is cut there)
  exchanges over its part above the soil only: ``emerged_links()`` gives the gas exchange each element's height, the
  centre of its triangles' parts above the Scene's soil surface (metafspm's ``CentreIn("air")``: NaN wholly below,
  no exchange), and its green area above it, the architecture's ``green_area`` times its green triangles' share
  above the soil (``FractionIn("air", weight="green")``, Adel's per-triangle green share: the senescence going down
  from the top), under a receiver-private name (``names={"green_area": "green_area_in_air"}``): the gas exchange
  reads it as its ``green_area``, Adel's ``green_area`` staying the whole element's.
- **The stomatal model** is a setting, ``stomatal_model="BWB"`` (the default), ``"Leuning"``, ``"Tuzet"`` or
  ``"hydraulics"``. Tuzet and hydraulics read the elements' ``leaf_water_potential`` (MPa), e.g. MECHA's
  ``water_potential`` (hPa at its Compartments): ``mecha_links()`` links it, its factor 1e-4, the mean over each
  element's Compartments weighted by MECHA's ``liquid_volume`` (the air spaces and pores, whose potential is the air's
  Kelvin equivalent, left out); a feedback lagged by one scene step, outside MECHA's solve. The plant model registers
  ``leaf_water_potential`` at the SubOrgans before coupling (an optional input is not registered until given).
- **Parameters** are per plant (a plant's scenario may change any), with their units, bounds and references;
  ``TOTAL_NITROGEN`` is the legacy's other parametrisation (the surfacic nitrogen with its structural part).
- **The nitrogen and carbohydrates** (``surfacic_nitrogen``, ``surfacic_NSC``) are inputs a CN model will give; until
  then their defaults (2 g m-2, the legacy's ``NA_0``, and 0). ``surfacic_photosynthetic_nitrogen``,
  ``surfacic_total_nitrogen`` and ``surfacic_nonstructural_carbohydrates`` compute them from the amounts as the
  legacy does.

Tutorials
=========

``openalea.cnwgrass.metafspm.simulation.simulate()`` computes, hour by hour, the gas exchange of a small stand of
Adel's wheat plants lit by Caribu (its plant model ``ExchangingWheat``), its results as tables. Two series of
notebooks use it, in ``example/tutorials/``:

- ``user/``, for beginners, without the formalism: a first run, changing a parameter, combining parameters, a
  sensitivity analysis (metafspm's ``variants``), a measured weather, the light, the plants and the leaf water
  potential of the models it runs with;
- ``developer/``: an overview, metafspm's component formalism on ``GasExchange`` (declarations, steps, forcings,
  settings, the stomatal models, the corrections) and its coupling in a plant model and a scene (Caribu's ``PARa``
  by the scene's identity links, its transpiration summed per organ into MECHA's ``leaf_transpiration``), and
  extending the component.

The equations
=============

.. automodule:: openalea.cnwgrass.metafspm.gas_exchange
   :no-members:

The component
=============

.. metafspm-component:: openalea.cnwgrass.metafspm.gas_exchange.GasExchange

The legacy and the port
=======================

The port is checked against the legacy ``Simulation`` (``test/test_gas_exchange.py``), with ``corrections=()``: bit
for bit on random elements for the four stomatal models and three weathers, with or without the NSC's inhibition,
with either nitrogen parametrisation, at the primitive scale, on the legacy's standalone example, and on Adel's
plants in a metafspm Scene (every element computed, two plants with their own parameters).

The corrections (``corrections="recommended"``, the default, applies them all):

.. list-table::
   :header-rows: 1

   * - correction
     - what it corrects
   * - ``vapour_pressure_slope``
     - the slope of the saturation vapour pressure at the air temperature, Tetens' derivative es(Ta) 17.4 · 239 /
       (Ta + 239)² with Ta in °C (the legacy adds 273.15 to Ta: a slope about four times too small at the fixed
       point's first iteration)
   * - ``compensation_point_units``
     - the CO2 compensation point Γ of the Leuning, Tuzet and hydraulics models in µmol mol-1, as Cs, its
       temperature in K (the legacy: Γ0 in mol mol-1 and Ta in °C minus Tref' in K, Γ about 0.002 µmol mol-1
       instead of about 27)
   * - ``primitive_aggregation``
     - at the primitive scale, every output the area-weighted mean of the primitives' (the legacy: Ag only, the
       others the last primitive's)

The legacy's other behaviours, and what the port does:

.. list-table::
   :header-rows: 1

   * - the legacy
     - the port
   * - parameters set by ``setattr`` on the module (``update_parameters``), shared by every simulation
     - parameters per plant, set by the scenarios
   * - the main stem only
     - every axis
   * - a previous Ci of 0 (its own output for an element without exchange) divides by 0 at the next step
     - the next step starts from the ambient CO2
   * - an air at 0 °C divides the relative change of Ts by 0 (when Ci converges at once, in the dark)
     - from 0 °C, Ts's change in °C
   * - a canopy height of 0 (ln(0)), a stem of null diameter (0 / 0) fail
     - no exchange without a canopy; a null diameter's boundary layer its limit, 0
   * - the facade's ``'pedoncule'`` for ``'peduncle'``: a peduncle's characteristic dimension its width
     - every non-lamina organ's its diameter
   * - the primitive scale with the default ``area_prim`` (0) fails
     - its inputs required
   * - a non-converged element printed
     - one ``RuntimeWarning`` per call, with their number
   * - a non-lamina element's Ag times ``EFFICENCY_STEM``, its An not (An ≠ Ag − Rd)
     - kept (``stem_efficiency``)
   * - the fixed point (Picard, 1 %, 30 iterations) oscillates for some elements (Ag dropping to 0 at a low Ci,
       then Ci back to Ca), most with the Leuning and Tuzet models
     - kept: the same iterations and the same last iterate
   * - mature hidden elements (CNW-Grass's facade: no height) exchange nothing, their Ts the apex's
     - the same by label (``hidden_labels``); the standalone ``Simulation`` computes any element with a height
