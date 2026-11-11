"""
CNW-Grass's gas exchange as a metafspm component: each shoot element's photosynthesis, stomatal conductance, surface
temperature and transpiration, from the light it absorbs, its nitrogen and the weather (the legacy
openalea.cnwgrass.gasexchange, its Simulation and its facade).

    GasExchange(data_structure=data_structure, time_step=3600.)       a FunctionalComponent at the SubOrgan scale (Adel's
                                                          elements: StemElement, LeafElement), on a plant model's
                                                          DataStructure; its parameters per plant
    GasExchange(..., stomatal_model="Tuzet")              the stomatal conductance limited by the leaf water
                                                          potential (an input, e.g. MECHA's)

Its steps follow the physics:
    _canopy_height    each plant's canopy height Zh: the highest of its exchanging elements' heights (the legacy's per
                      axis), 0 without any
    _capacities       each element's photosynthetic capacities at 25 °C from its surfacic nitrogen Na:
                      Vc_max25, Jmax25, TPU25, Rdark25 = S · (Na − Na_min), and the electron yield α = S_α · Na + β
    _exchange         for each exchanging element, the fixed point of the photosynthesis, the stomatal conductance,
                      the internal CO2 and the energy balance (on Ci and Ts, from Ci the previous step's and Ts = Ta);
                      Ag, An, Rd, gs, Ci, Ts, Tr, the area it exchanges over (exchanging_area: the green area it
                      reads) and its transpiration Tr · exchanging_area

The equations (units: µmol m-2 s-1 for the CO2 fluxes and PAR, µmol mol-1 for the CO2 concentrations, °C for the
temperatures unless in K, kPa for the vapour pressures):

Temperature dependence of a parameter p (Braune et al. 2009; Bernacchi et al. 2001 for Kc, Ko, Rdark), T in K,
Tref = 298.15 K, R in kJ mol-1 K-1::

    p(T) = p25 · exp(ΔHa (T − Tref) / (R Tref T)) · [1 + exp((Tref ΔS − ΔHd) / (R Tref))]
                                                   / [1 + exp((T ΔS − ΔHd) / (R T))]     (the bracket for Vc_max,
                                                                                         Jmax and TPU only)

Photosynthesis (Farquhar et al. 1980; Müller et al. 2005, Braune et al. 2009, Evers et al. 2010), at Ts and Ci::

    Ac = Vc_max (Ci − Γ*) / (Ci + Kc (1 + O2 / Ko))                          RuBisCO-limited
    J  = (Jmax + α PAR − sqrt((Jmax + α PAR)² − 4 θ α PAR Jmax)) / (2 θ)
    Aj = J (Ci − Γ*) / (4 Ci + 8 Γ*)                                         RuBP regeneration-limited
    Vomax = Vc_max Ko Γ* / (0.5 Kc O2),  Vo = Vomax O2 / (O2 + Ko (1 + Ci / Kc))
    Ap = (1 − Γ* / Ci) (3 TPU + Vo)                                          TPU-limited
    Ag = min(Ac, Aj) (1 − I(NSC))   with NSC_retroinhibition (the default), I(NSC) = min(Imax (NSC − NSCmin) /
                                    (K + NSC − NSCmin), 1) above NSCmin, else 0 (Azcón-Bieto 1983)
    Ag = min(Ac, Aj, Ap)            without
    Rd = Rdark (Rd_A + (1 − Rd_A) Rd_B^(PAR / Rd_C))                         (Müller et al. 2005, eq. 19)
    An = Ag − Rd, and Ag = An = 0 when Ag <= 0

Stomatal conductance to water vapour gs (mol m-2 s-1), at the surface CO2 Cs = Ca − An · 1.37 / gb::

    BWB         gs = g0 + m Ag RH / Cs,  m = δ1 Na^δ2                       (Ball et al. 1987; Braune et al. 2009)
    Leuning     gs = g0 + 1.6 m' Ag fD / (Cs − Γ),  fD = 1 / (1 + VPD / D0)                        (Leuning 1995)
    Tuzet       gs = g0 + 1.6 m' Ag fψ / (Cs − Γ),  fψ = 1 / (1 + (ψ / ψref)^n)     (Tuzet et al. 2003; Alabasha)
    hydraulics  gs = g0 + 1.6 m' Ag fD fψ / (Cs − Γ)                                (Leuning and Tuzet together)
    Γ = Γ0 (1 + Γ1 (T − Tref') + Γ2 (T − Tref')²)       (corrected: Γ0 in µmol mol-1 and T in K; see CORRECTIONS)

Internal CO2::

    Ci = Ca − An (1.6 / gs + 1.37 / gb)

Energy balance (Penman-Monteith; Campbell and Norman 1998 for the wind in the canopy), with the canopy height Zh,
the wind Ur at the reference height zr (2 m), the element's height z and characteristic dimension w (a blade's
width, a stem's or a sheath's diameter)::

    d = 0.7 Zh,  z0 = 0.1 Zh,  u* = κ Ur / ln((zr − d) / z0),  u = u* / κ ln((Zh − d) / z0) exp(a (z / Zh − 1))
    rbh = 154 sqrt(w / u)                    a lamina (a horizontal plane in forced convection)
    rbh = w / (1.2e-5 (u w / 1.5e-5)^0.47)   a cylinder (a stem, a sheath; 0 for a null diameter, its limit)
    ra  = ln((zr − d) / z0)² / (κ² Ur)
    Rn  = PAR · 1.48 / 4.55                                  (the absorbed global radiation, W m-2)
    s   = es(Ta) 17.4 · 239 / (Ta + 239)²  at the first iteration (Ts = Ta), else (es(Ts) − es(Ta)) / (Ts − Ta)
    E   = max(0, (s Rn + ρcp VPD / (rbh + ra)) / (λ (s + γ (rbw + ra + rs) / (rbh + ra))))     (mm s-1)
    Ts  = Ta + (rbh + ra)(Rn − λ E) / ρcp,  rbw = 0.96 rbh,  rs = Patm / (gs R (Ts + 273.15)),
    es(T) = 0.611 exp(17.4 T / (239 + T)),  VPD = es(Ta) (1 − RH)
    Tr  = E · 1e6 / 18                                       (mmol m-2 s-1)

The fixed point stops when Ci and Ts change by less than convergence_tolerance (relative, 1 %), or after
maximum_iterations (30). A non-lamina element's Ag is then multiplied by stem_efficiency (0.78; An and Rd are not,
as in the legacy). At the primitive scale (primitive_scale), each of the element's primitives (PARa_primitives,
their areas primitive_areas) is solved in turn, its Ci and Ts the previous primitive's to start, and the element's
outputs are their area-weighted means.

Which elements: those labelled exchanging_labels (Adel's StemElement and LeafElement, the legacy's LeafElement1)
with a green area and a height; an element labelled hidden_labels (inside the sheaths' tube), or exchanging without
green area or height, or with a null surfacic nitrogen, exchanges nothing (0), its Ts the apex temperature
(SAM_temperature, else the air's) and its Ci 0 (as the legacy; the next step starts it again from Ca). The other
elements (roots, in a DataStructure holding the whole plant) get 0 and Ts NaN. Every axis is computed (the legacy:
the main stem only).

Its inputs: PARa (CaribuLight), green_area, height and width (the architecture: Adel's ShootArchitecture; where an
architecture crosses the soil, an assembly gives the gas exchange the green area of the element's part above it, as
Caribu's light is cut there: a link of a receiver-private name, metafspm's Link(names={"green_area": ...}), from the
architecture's green_area and the green triangles' share above the soil, data_structure.geometry; Adel's green_area untouched), the
organs' diameter (Adel's), the surfacic nitrogen and NSC (a CN model's, later; defaults until then), the leaf water
potential leaf_water_potential (MPa; for Tuzet and hydraulics: MECHA's water_potential, hPa at its Compartments, the
mean of its liquid per element, lagged by one step: mecha_links()). Adel has no element height: element_height (the
mean of its base's and its end's z, a blade element's along its blade) links it, by a translator formula
(adel_links()); CNW-Grass's facade took the mean height of the element's triangles (mean_triangle_heights).

The legacy's behaviours that this port corrects are named in CORRECTIONS (all applied by default, "recommended";
corrections=() reproduces the legacy bit for bit). Fixed without a name, since the legacy fails there: an element
whose previous Ci is 0 (the legacy divides by it), an air temperature of 0 °C (the relative change of Ts divided by
0), a canopy height of 0 (ln(0)), a cylinder of null diameter (0 / 0), the primitive scale without areas; and its
parameters are the plants' (the legacy's are module globals set by setattr, shared by every simulation), the stems'
characteristic dimension their diameter whatever their organ (the facade's 'pedoncule' for 'peduncle').
"""
import warnings
from collections import namedtuple
from dataclasses import dataclass, fields
from math import exp, log, sqrt

import numpy as np

from openalea.metafspm.coupling.component import (FunctionalComponent, forcing_variable, input_variable, parameter,
                                                   setting, state_variable)
from openalea.metafspm.coupling.geometry import CentreIn, FractionIn
from openalea.metafspm.coupling.translator import Translator
from openalea.metafspm.data_structure.configs import ScalesConfig as scales
from openalea.metafspm.solve.decorator import rate, state
from openalea.metafspm.testing.corrections import Corrections

ELEMENT = dict(scale=scales.SubOrgan)
PER_AREA = dict(scale=scales.SubOrgan, kind="intensive", weight="exchanging_area")
BY = "GasExchange"
BRAUNE = "Braune et al. (2009); Evers et al. (2010); Barillot et al. (2016); Gauthier et al. (2021)"
BERNACCHI = "Bernacchi et al. (2001)"
MULLER = "Müller et al. (2005)"
TUZET = "Tuzet et al. (2003)"
ALABASHA = "Alabasha et al. (2019)"
CAMPBELL = "Campbell and Norman (1998)"

# ── constants of the equations (the legacy's parameters module) ──────────────
VON_KARMAN = 0.40                   # von Kármán's constant (adim)
LATENT_HEAT = 2260E3                # latent heat of vaporisation of water (J kg-1)
AIR_HEAT_CAPACITY = 1256            # volumetric heat capacity of air (J m-3 K-1)
PSYCHROMETRIC_CONSTANT = 66E-3      # (kPa K-1)
GAS_CONSTANT = 8.3144               # (J mol-1 K-1)
KELVIN_DEGREE = 273.15              # 0 °C in K
ATMOSPHERIC_PRESSURE = 1.01325E5    # (Pa)
WATER_MOLAR_MASS = 18               # (g mol-1)
NITROGEN_MOLAR_MASS = 14            # (g mol-1)
SATURATION_A, SATURATION_B, SATURATION_C = 239, 17.4, 0.611   # Tetens' saturation vapour pressure (°C, adim, kPa)
SATURATION_EXPONENT = 2
RESISTANCE_EXPONENT = 2
CO2_DIFFUSIVITY_RATIO = 1.6         # gs for water vapour / gs for CO2
BOUNDARY_CO2_RATIO = 1.37           # gb for water vapour / gb for CO2 (1.6^(2/3))
ELECTRON_EXPONENT, ELECTRON_A, ELECTRON_B = 2, 4, 2           # J's non-rectangular hyperbola
RUBP_A, RUBP_B = 4, 8               # Aj = J (Ci − Γ*) / (4 Ci + 8 Γ*)
OXYGENATION_FACTOR = 0.5            # Vomax = Vc_max Ko Γ* / (0.5 Kc O2)
TPU_FACTOR = 3                      # Ap = (1 − Γ* / Ci)(3 TPU + Vo)

STOMATAL_MODELS = ("BWB", "Leuning", "Tuzet", "hydraulics")

NAMES = {
    "vapour_pressure_slope": "the slope of the saturation vapour pressure at the air temperature, Tetens' "
                             "derivative es(Ta) 17.4 · 239 / (Ta + 239)² with Ta in °C (the legacy adds 273.15 to Ta: "
                             "a slope about four times too small at the fixed point's first iteration)",
    "compensation_point_units": "the CO2 compensation point Γ of the Leuning, Tuzet and hydraulics models in µmol "
                                "mol-1, as Cs, its temperature in K (the legacy: Γ0 in mol mol-1 and Ta in °C minus "
                                "Tref' in K, Γ about 0.002 µmol mol-1 instead of about 27)",
    "primitive_aggregation": "at the primitive scale, every output the area-weighted mean of the primitives' (the "
                             "legacy: Ag only, the others the last primitive's)",
}
CORRECTIONS = Corrections(NAMES)

# The legacy's other parametrisation (parameters.SurfacicProteins False: the nitrogen dependences on the total surfacic
# nitrogen, structural included; Barillot et al. 2016, Gauthier et al. 2020), as a scenario
TOTAL_NITROGEN = {"Vc_max25_nitrogen_minimum": 0.17, "Jmax25_nitrogen_minimum": 0.17, "TPU25_nitrogen_minimum": 0.229,
                  "Rdark25_nitrogen_minimum": 0.118, "electron_yield_intercept": 0.2101}


# ── the surfacic contents, from a CN model's amounts ─────────────────────────

def surfacic_photosynthetic_nitrogen(proteins, green_area, conversion=1.06):
    """The surfacic non-structural nitrogen of the Farquhar model (g m-2) from the element's proteins (µmol N) and
    green area (m2): its surfacic proteins (g N m-2) times *conversion* (the legacy's Psurf_to_SLNnonstruct, its
    default; SurfacicProteins True)."""
    return proteins * 1E-6 * NITROGEN_MOLAR_MASS / green_area * conversion


def surfacic_total_nitrogen(nitrates, amino_acids, proteins, structural_nitrogen, green_area):
    """The surfacic nitrogen, structural included (g m-2), from the amounts of nitrates, amino acids and proteins
    (µmol N) and the structural nitrogen (mtg): the legacy's with SurfacicProteins False (with TOTAL_NITROGEN)."""
    return ((nitrates + amino_acids + proteins) * 1E-6 * NITROGEN_MOLAR_MASS + structural_nitrogen) / green_area


def surfacic_nonstructural_carbohydrates(sucrose, starch, fructan, green_area):
    """The surfacic non-structural carbohydrates (µmol C m-2) from the amounts of sucrose, starch and fructan (µmol C)
    and the green area (m2)."""
    return (sucrose + starch + fructan) / green_area


# ── the element's height, from the architecture ──────────────────────────────

def element_height(z1, z2):
    """An element's height (m): the mean of its base's and its end's heights (Adel's drawn ends, display_z1,
    display_z2), for a translator formula."""
    return (np.asarray(z1, dtype=float) + np.asarray(z2, dtype=float)) / 2


def mean_triangle_heights(triangles, n_elements):
    """Each element's height as CNW-Grass's facade took it: the mean over its triangles of their vertices' mean z
    (astk's get_height, then nanmean), NaN without triangles (a coupling.geometry.Triangles, n_elements elements)."""
    heights = triangles.vertices[..., 2].mean(axis=1)
    total = np.bincount(triangles.owner, weights=heights, minlength=n_elements)
    count = np.bincount(triangles.owner, minlength=n_elements)
    return np.divide(total, count, out=np.full(n_elements, np.nan), where=count > 0)


def adel_links(receiver="GasExchange", provider="ShootArchitecture") -> Translator:
    """The links from Adel's ShootArchitecture that sharing by name does not make: the elements' height
    (element_height of where Adel draws them, display_z1, display_z2: its topological ends z1, z2 stack the hidden
    parts the turtle draws on the visible ones). width, diameter and green_area are shared by their names."""
    return Translator().link(receiver, "height", provider, ("display_z1", "display_z2"), formula=element_height)


def emerged_links(receiver="GasExchange", provider="ShootArchitecture") -> Translator:
    """
    adel_links for a shoot crossing the soil surface (sown, on a crown below it: the Scene's soil_surface), the gas
    exchange of each element's part above it only, as Caribu's light is cut there (metafspm's data_structure.geometry, the
    domain "air"): its height, the centre of its triangles' parts above the soil (CentreIn: NaN wholly below, no gas
    exchange); its green area, the architecture's green_area times its green triangles' share above the soil
    (FractionIn, weighted by Adel's per-triangle "green": the senescence going down from the top; 0 without any),
    under a receiver-private name (green_area_in_air): the gas exchange reads it as its green_area, the
    architecture's green_area untouched. width, diameter: shared by their names.
    """
    return (Translator()
            .link(receiver, "height", provider, (CentreIn("air", axis="z"),), formula=_height)
            .link(receiver, "green_area", provider, ("green_area", FractionIn("air", weight="green", empty=0.)),
                  formula=np.multiply, names={"green_area": "green_area_in_air"}))


def _height(centre_height):
    """The elements' heights from the geometric source's (a module-level function: a checkpoint keeps the
    DataStructure's derivations)."""
    return np.asarray(centre_height, dtype=float)


def mecha_links(receiver="GasExchange", provider="MechaHydraulics") -> Translator:
    """The link from MECHA's hydraulics (MechaHydraulics, on an anatomy-mode DataStructure: its nodes the
    Compartments): each element's leaf_water_potential (MPa), the mean of MECHA's water_potential (hPa, 1e-4 MPa) over
    the element's Compartments weighted by their liquid_volume (the walls and protoplasts, not the air spaces and
    pores, whose water potential is the air's Kelvin equivalent). Read by the Tuzet and hydraulics stomatal models,
    MECHA's of the step before when the gas exchange runs first (a lagged feedback, outside MECHA's solve)."""
    return Translator().link(receiver, "leaf_water_potential", provider, {"water_potential": 1e-4},
                             source_scale=scales.Compartment, scale=scales.SubOrgan, aggregation="weighted_mean",
                             weight="liquid_volume")


# ── the equations, for one element (the legacy's model, its numerics kept) ───

def temperature_dependence(value_25, temperature, activation_enthalpy, reference_temperature, entropy=None,
                           deactivation_enthalpy=None):
    """A photosynthetic parameter at *temperature* (°C) from its value at 25 °C: its activation (Arrhenius) and,
    with an *entropy* and a *deactivation_enthalpy* (kJ mol-1 K-1, kJ mol-1), its deactivation."""
    temperature_kelvin = temperature + KELVIN_DEGREE
    activation = exp((activation_enthalpy * (temperature_kelvin - reference_temperature))
                     / (GAS_CONSTANT * 1E-3 * reference_temperature * temperature_kelvin))
    if entropy is None:
        deactivation = 1
    else:
        deactivation = (1 + exp((reference_temperature * entropy - deactivation_enthalpy)
                                / (reference_temperature * GAS_CONSTANT * 1E-3))) / (
            1 + exp((temperature_kelvin * entropy - deactivation_enthalpy) / (temperature_kelvin * GAS_CONSTANT * 1E-3)))
    return value_25 * activation * deactivation


def inhibition_by_carbohydrates(surfacic_NSC, p):
    """The relative decrease of Ag by the surfacic non-structural carbohydrates (µmol C m-2; Azcón-Bieto 1983)."""
    if surfacic_NSC <= p.NSC_inhibition_threshold:
        return 0
    return min(p.NSC_inhibition_maximum * (surfacic_NSC - p.NSC_inhibition_threshold)
               / (p.NSC_inhibition_constant + surfacic_NSC - p.NSC_inhibition_threshold), 1)


def photosynthesis(PAR, capacities, surfacic_NSC, surface_temperature, Ci, p, retroinhibition):
    """(Ag, An, Rd) (µmol m-2 s-1) of Farquhar's model at the surface temperature (°C) and the internal CO2 (µmol
    mol-1), the capacities at 25 °C (Vc_max25, Jmax25, TPU25, Rdark25, α) the element's."""
    Vc_max25, Jmax25, TPU25, Rdark25, alpha = capacities
    reference = p.reference_temperature
    Kc = temperature_dependence(p.KC25, surface_temperature, p.Kc_activation_enthalpy, reference)
    Ko = temperature_dependence(p.KO25, surface_temperature, p.Ko_activation_enthalpy, reference)
    Gamma = temperature_dependence(p.GAMMA25, surface_temperature, p.Gamma_activation_enthalpy, reference)

    Vc_max = temperature_dependence(Vc_max25, surface_temperature, p.Vc_max_activation_enthalpy, reference,
                                    p.Vc_max_entropy, p.Vc_max_deactivation_enthalpy)
    Ac = (Vc_max * (Ci - Gamma)) / (Ci + Kc * (1 + p.O2 / Ko))

    Jmax = temperature_dependence(Jmax25, surface_temperature, p.Jmax_activation_enthalpy, reference,
                                  p.Jmax_entropy, p.Jmax_deactivation_enthalpy)
    J = ((Jmax + alpha * PAR) - sqrt((Jmax + alpha * PAR) ** ELECTRON_EXPONENT
                                     - ELECTRON_A * p.THETA * alpha * PAR * Jmax)) / (ELECTRON_B * p.THETA)
    Aj = (J * (Ci - Gamma)) / (RUBP_A * Ci + RUBP_B * Gamma)

    if retroinhibition:
        Ag = min(Ac, Aj) * (1 - inhibition_by_carbohydrates(surfacic_NSC, p))
    else:
        TPU = temperature_dependence(TPU25, surface_temperature, p.TPU_activation_enthalpy, reference,
                                     p.TPU_entropy, p.TPU_deactivation_enthalpy)
        Vomax = (Vc_max * Ko * Gamma) / (OXYGENATION_FACTOR * Kc * p.O2)
        Vo = (Vomax * p.O2) / (p.O2 + Ko * (1 + Ci / Kc))
        Ap = (1 - Gamma / Ci) * (TPU_FACTOR * TPU + Vo)
        Ag = min(Ac, Aj, Ap)

    Rdark = temperature_dependence(Rdark25, surface_temperature, p.Rdark_activation_enthalpy, reference)
    Rd = Rdark * (p.light_respiration_fraction + (1 - p.light_respiration_fraction)
                  * p.light_respiration_base ** (PAR / p.light_respiration_scale))
    if Ag <= 0:             # Ci below Γ*, or Na below its minimum: no assimilation (Farquhar 1980; Caemmerer 2000)
        return 0., 0., Rd
    return Ag, Ag - Rd, Rd


def saturation_vapour_pressure(temperature):
    """Tetens' saturation vapour pressure (kPa) at *temperature* (°C)."""
    return SATURATION_C * exp((SATURATION_B * temperature) / (SATURATION_A + temperature))


def surface_CO2(ambient_CO2, An, p):
    """The CO2 at the element's surface (µmol mol-1; Prieto et al. 2012)."""
    return ambient_CO2 - An * (BOUNDARY_CO2_RATIO / p.boundary_layer_conductance)


def compensation_point(air_temperature, p, corrected):
    """Γ of the Leuning and Tuzet models (µmol mol-1 corrected, else the legacy's, see CORRECTIONS)."""
    if corrected:
        temperature, reference = air_temperature + KELVIN_DEGREE, p.compensation_point_reference * 1E6
    else:
        temperature, reference = air_temperature, p.compensation_point_reference
    difference = temperature - p.compensation_point_temperature
    return reference * (1 + p.compensation_point_linear * difference
                        + p.compensation_point_quadratic * (difference * difference))


def water_potential_factor(water_potential, p):
    """Tuzet's reduction of the conductance by the leaf water potential (MPa)."""
    return 1 / (1 + (water_potential / p.reference_water_potential) ** p.water_potential_shape)


def vapour_pressure_factor(air_temperature, relative_humidity, p):
    """Leuning's reduction of the conductance by the air's vapour pressure deficit."""
    saturation = saturation_vapour_pressure(air_temperature)
    deficit = saturation - relative_humidity * saturation
    return 1 / (1 + deficit / p.vapour_pressure_deficit_scale)


def stomatal_conductance(model, Ag, An, surfacic_nitrogen, ambient_CO2, relative_humidity, air_temperature,
                         water_potential, p, corrected_compensation):
    """gs, the stomatal conductance to water vapour (mol m-2 s-1), of *model* (STOMATAL_MODELS)."""
    Cs = surface_CO2(ambient_CO2, An, p)
    if model == "BWB":
        m = p.BWB_slope_factor * surfacic_nitrogen ** p.BWB_slope_exponent
        return p.minimum_stomatal_conductance + m * ((Ag * relative_humidity) / Cs)
    gamma = compensation_point(air_temperature, p, corrected_compensation)
    if model == "Leuning":
        factor = vapour_pressure_factor(air_temperature, relative_humidity, p)
        return p.minimum_stomatal_conductance + (CO2_DIFFUSIVITY_RATIO * p.stomatal_slope * Ag * factor) / (Cs - gamma)
    if model == "Tuzet":
        factor = water_potential_factor(water_potential, p)
        return p.minimum_stomatal_conductance + (CO2_DIFFUSIVITY_RATIO * p.stomatal_slope * Ag * factor) / (Cs - gamma)
    deficit_factor = vapour_pressure_factor(air_temperature, relative_humidity, p)
    potential_factor = water_potential_factor(water_potential, p)
    return p.minimum_stomatal_conductance + (CO2_DIFFUSIVITY_RATIO * p.stomatal_slope * Ag * deficit_factor
                                             * potential_factor) / (Cs - gamma)


def internal_CO2(ambient_CO2, An, conductance, p):
    """Ci (µmol mol-1), from An and gs."""
    return ambient_CO2 - An * ((CO2_DIFFUSIVITY_RATIO / conductance) + (BOUNDARY_CO2_RATIO
                                                                        / p.boundary_layer_conductance))


def organ_temperature(dimension, height, canopy_height, wind_speed, PAR, conductance, air_temperature,
                      surface_temperature, relative_humidity, lamina, p, reference_height, minimum_wind_speed,
                      corrected_slope):
    """(Ts in °C, E in mm s-1): the element's energy balance (Penman-Monteith) at its surface temperature of the
    iteration before (Ta at the first one)."""
    displacement = p.displacement_height_ratio * canopy_height
    roughness = p.roughness_length_ratio * canopy_height
    wind_speed = max(wind_speed, minimum_wind_speed)
    friction_velocity = (wind_speed * VON_KARMAN) / log((reference_height - displacement) / roughness)
    canopy_wind = (friction_velocity / VON_KARMAN) * log((canopy_height - displacement) / roughness)
    wind = canopy_wind * exp(p.wind_attenuation * (height / canopy_height - 1))

    if lamina:              # a horizontal plane in forced convection
        heat_resistance = p.lamina_boundary_layer_factor * sqrt(dimension / wind)
    elif dimension > 0:     # a vertical cylinder in forced convection
        heat_resistance = dimension / (p.cylinder_boundary_layer_factor * ((wind * dimension)
                                       / p.cylinder_boundary_layer_viscosity) ** p.cylinder_boundary_layer_exponent)
    else:                   # its limit for a null diameter
        heat_resistance = 0.
    turbulence_resistance = 1 / (VON_KARMAN ** RESISTANCE_EXPONENT * wind_speed) * (
        log((reference_height - displacement) / roughness)) ** RESISTANCE_EXPONENT

    net_radiation = (PAR * p.PAR_to_global) / p.PAR_photon_energy
    saturation_air = saturation_vapour_pressure(air_temperature)
    vapour_pressure = relative_humidity * saturation_air
    if surface_temperature == air_temperature:
        if corrected_slope:
            slope = ((SATURATION_B * SATURATION_A) / (air_temperature + SATURATION_A) ** SATURATION_EXPONENT) \
                * saturation_air
        else:
            air_kelvin = air_temperature + KELVIN_DEGREE
            slope = ((SATURATION_B * SATURATION_A) / (air_kelvin + SATURATION_A) ** SATURATION_EXPONENT) \
                * saturation_air
    else:
        saturation_surface = saturation_vapour_pressure(surface_temperature)
        surface_kelvin, air_kelvin = surface_temperature + KELVIN_DEGREE, air_temperature + KELVIN_DEGREE
        slope = (saturation_surface - saturation_air) / (surface_kelvin - air_kelvin)

    deficit = saturation_air - vapour_pressure
    water_resistance = p.boundary_layer_ratio * heat_resistance
    physical_conductance = (conductance * GAS_CONSTANT * (surface_temperature + KELVIN_DEGREE)) / ATMOSPHERIC_PRESSURE
    stomatal_resistance = 1 / physical_conductance
    evaporation = max(0., (slope * net_radiation + (AIR_HEAT_CAPACITY * deficit) / (heat_resistance
                                                                                    + turbulence_resistance))
                      / (LATENT_HEAT * (slope + PSYCHROMETRIC_CONSTANT * ((water_resistance + turbulence_resistance
                                                                          + stomatal_resistance)
                                                                         / (heat_resistance + turbulence_resistance)))))
    surface_temperature = air_temperature + ((heat_resistance + turbulence_resistance)
                                             * (net_radiation - LATENT_HEAT * evaporation)) / AIR_HEAT_CAPACITY
    return surface_temperature, evaporation


def _temperature_converged(surface_temperature, previous, tolerance):
    """Ts's relative change below *tolerance* (the legacy's); from 0 °C, its change in °C (the legacy divides by
    0)."""
    if previous == 0:
        return abs(surface_temperature - previous) < tolerance
    return abs((surface_temperature - previous) / previous) < tolerance


Weather = namedtuple("Weather", "air_temperature ambient_CO2 relative_humidity wind_speed")
Options = namedtuple("Options", "stomatal_model retroinhibition tolerance maximum_iterations reference_height "
                                "minimum_wind_speed corrections")


def element_exchange(PAR, surfacic_nitrogen, surfacic_NSC, capacities, dimension, height, canopy_height, lamina,
                     water_potential, Ci, surface_temperature, weather, p, options):
    """The fixed point of one element (or primitive), from its Ci and Ts: (Ag, An, Rd, gs, Ci, Ts, Tr, converged),
    Tr in mmol m-2 s-1, Ag the efficiency's for a non-lamina element."""
    corrected_compensation = "compensation_point_units" in options.corrections
    corrected_slope = "vapour_pressure_slope" in options.corrections
    count = 0
    while True:
        previous_Ci, previous_temperature = Ci, surface_temperature
        Ag, An, Rd = photosynthesis(PAR, capacities, surfacic_NSC, surface_temperature, Ci, p, options.retroinhibition)
        conductance = stomatal_conductance(options.stomatal_model, Ag, An, surfacic_nitrogen, weather.ambient_CO2,
                                           weather.relative_humidity, weather.air_temperature, water_potential, p,
                                           corrected_compensation)
        Ci = internal_CO2(weather.ambient_CO2, An, conductance, p)
        surface_temperature, evaporation = organ_temperature(
            dimension, height, canopy_height, weather.wind_speed, PAR, conductance, weather.air_temperature,
            surface_temperature, weather.relative_humidity, lamina, p, options.reference_height,
            options.minimum_wind_speed, corrected_slope)
        count += 1
        converged = (abs((Ci - previous_Ci) / previous_Ci) < options.tolerance
                     and _temperature_converged(surface_temperature, previous_temperature, options.tolerance))
        if converged or count >= options.maximum_iterations:
            break
    transpiration = (evaporation * 1E6) / WATER_MOLAR_MASS        # mm s-1 (kg m-2 s-1) to mmol m-2 s-1
    if not lamina:
        Ag = Ag * p.stem_efficiency
    return Ag, An, Rd, conductance, Ci, surface_temperature, transpiration, converged


def primitives_exchange(PARs, areas, surfacic_nitrogen, surfacic_NSC, capacities, dimension, height, canopy_height,
                        lamina, water_potential, Ci, weather, p, options):
    """An element's outputs over its primitives (one PAR and area each; [PARa], [1.] at the element's scale), each
    started from the previous one's Ci and Ts: (Ag, An, Rd, gs, Ci, Ts, Tr, converged), Ag their area-weighted mean,
    the others too with the correction primitive_aggregation (else the last primitive's)."""
    surface_temperature = weather.air_temperature
    results = []
    for PAR in PARs:
        result = element_exchange(PAR, surfacic_nitrogen, surfacic_NSC, capacities, dimension, height, canopy_height,
                                  lamina, water_potential, Ci, surface_temperature, weather, p, options)
        Ci, surface_temperature = result[4], result[5]
        results.append(result)
    total_area = sum(areas)
    weighted = [sum([result[k] * area for result, area in zip(results, areas)]) / total_area for k in range(7)] \
        if "primitive_aggregation" in options.corrections else list(results[-1][:7])
    weighted[0] = sum([result[0] * area for result, area in zip(results, areas)]) / total_area
    return tuple(weighted) + (all(result[7] for result in results),)


# ── the component ────────────────────────────────────────────────────────────

def _parameter(unit, description, default, references=BRAUNE, **options):
    """A per-plant parameter of the gas exchange (by GasExchange)."""
    return parameter(unit=unit, description=description, by=BY, default=default, references=references, **options)


@dataclass
class GasExchange(FunctionalComponent):
    """Each shoot element's gas exchange: its canopy's height (_canopy_height), its photosynthetic capacities from its
    nitrogen (_capacities), and the coupled photosynthesis, stomatal conductance and energy balance (_exchange)."""

    # --- Inputs ---
    PARa: float = input_variable(unit="µmol m-2 s-1", description="The PAR the element absorbs per unit area",
                                 by="CaribuLight", default=0., kind="intensive", **ELEMENT)
    green_area: float = input_variable(unit="m2", description="The element's green area (0: no gas exchange)",
                                       by="ShootArchitecture", default=0., **ELEMENT)
    height: float = input_variable(
        unit="m", description="The element's height above the soil", value_comment="NaN: no geometry, no gas "
        "exchange; Adel's: element_height(display_z1, display_z2) (adel_links()); CNW-Grass's facade: its "
        "triangles' mean height",
        by="ShootArchitecture", default=np.nan, kind="descriptor", **ELEMENT)
    width: float = input_variable(unit="m", description="The lamina element's width (its characteristic dimension "
                                  "for the heat transfer)", by="ShootArchitecture", default=0., kind="descriptor",
                                  **ELEMENT)
    diameter: float = input_variable(unit="m", description="The organ's diameter (a stem's or a sheath's "
                                     "characteristic dimension for the heat transfer)", by="ShootArchitecture",
                                     default=0., kind="descriptor", scale=scales.Organ)
    surfacic_nitrogen: float = input_variable(
        unit="g m-2", unit_comment="of nitrogen", description="The element's surfacic nitrogen Na (non-structural: "
        "surfacic_photosynthetic_nitrogen; with TOTAL_NITROGEN, structural included: surfacic_total_nitrogen)",
        value_comment="the legacy's NA_0 until a CN model gives it; 0: no gas exchange", by="CNMetabolism",
        default=2., kind="intensive", **ELEMENT)
    surfacic_NSC: float = input_variable(
        unit="µmol m-2", unit_comment="of C", description="The element's surfacic non-structural carbohydrates "
        "(sucrose, starch, fructan: surfacic_nonstructural_carbohydrates), which inhibit its photosynthesis",
        by="CNMetabolism", default=0., kind="intensive", **ELEMENT)
    leaf_water_potential: float = input_variable(
        unit="MPa", description="The element's leaf water potential (the Tuzet and hydraulics stomatal models)",
        value_comment="the legacy's default -0.1; MECHA's at the step before (a lagged feedback: mecha_links())",
        by="MECHA",
        default=-0.1, kind="intensive", optional=True, **ELEMENT)
    SAM_temperature: float = input_variable(
        unit="°C", description="The axis' apex temperature: the surface temperature of its elements that exchange "
        "nothing (hidden in the sheaths)", value_comment="NaN: the air's", by="Morphogenesis", default=np.nan,
        kind="intensive", optional=True, scale=scales.Axis)
    PARa_primitives: object = input_variable(
        unit="µmol m-2 s-1", description="The PAR each primitive of the element absorbs (primitive_scale)",
        by="CaribuLight", default=None, dtype="object", optional=True, **ELEMENT)
    primitive_areas: object = input_variable(
        unit="m2", description="The area of each primitive of the element (primitive_scale)", by="CaribuLight",
        default=None, dtype="object", optional=True, **ELEMENT)

    # --- Forcings ---
    air_temperature: float = forcing_variable(unit="°C", description="The air temperature")
    relative_humidity: float = forcing_variable(unit="adim", description="The air's relative humidity (a fraction)")
    ambient_CO2: float = forcing_variable(unit="µmol mol-1", description="The air's CO2 concentration")
    wind_speed: float = forcing_variable(unit="m s-1", description="The wind speed at the reference height (above "
                                         "the canopy, reference_height)")

    # --- State variables ---
    canopy_height: float = state_variable(unit="m", description="The plant's canopy height: its exchanging "
                                          "elements' highest", kind="descriptor", default=0., scale=scales.Plant)
    Vc_max25: float = state_variable(unit="µmol m-2 s-1", description="The element's maximal carboxylation rate at "
                                     "25 °C", default=0., **PER_AREA)
    Jmax25: float = state_variable(unit="µmol m-2 s-1", unit_comment="of electrons", description="The element's "
                                   "maximal electron transport rate at 25 °C", default=0., **PER_AREA)
    TPU25: float = state_variable(unit="µmol m-2 s-1", description="The element's triose phosphate utilisation rate "
                                  "at 25 °C", default=0., **PER_AREA)
    Rdark25: float = state_variable(unit="µmol m-2 s-1", description="The element's dark respiration at 25 °C",
                                    default=0., **PER_AREA)
    electron_yield: float = state_variable(unit="mol mol-1", unit_comment="of electrons per photon", description="The "
                                           "element's electron transport efficiency α", default=0., **PER_AREA)
    Ag: float = state_variable(unit="µmol m-2 s-1", unit_comment="of CO2", description="The element's gross "
                               "assimilation", default=0., **PER_AREA)
    An: float = state_variable(unit="µmol m-2 s-1", unit_comment="of CO2", description="The element's net "
                               "assimilation", default=0., **PER_AREA)
    Rd: float = state_variable(unit="µmol m-2 s-1", unit_comment="of CO2", description="The element's mitochondrial "
                               "respiration in the light", default=0., **PER_AREA)
    gs: float = state_variable(unit="mol m-2 s-1", unit_comment="of water vapour", description="The element's "
                               "stomatal conductance", default=0., **PER_AREA)
    Ci: float = state_variable(unit="µmol mol-1", description="The element's internal CO2 concentration",
                               value_comment="NaN or 0: the next step starts from the ambient CO2", default=np.nan,
                               **PER_AREA)
    Ts: float = state_variable(unit="°C", description="The element's surface temperature", value_comment="NaN: not "
                               "a shoot element", default=np.nan, **PER_AREA)
    Tr: float = state_variable(unit="mmol m-2 s-1", unit_comment="of water", description="The element's "
                               "transpiration rate", default=0., **PER_AREA)
    exchanging_area: float = state_variable(unit="m2", description="The element's area exchanging with the air: "
                                            "its green area (the one it reads: an assembly's, above the soil)",
                                            value_comment="0: an element exchanging nothing", default=0.,
                                            kind="extensive", **ELEMENT)
    transpiration: float = state_variable(unit="mmol s-1", unit_comment="of water", description="The element's "
                                          "transpiration, Tr · exchanging_area", default=0., kind="extensive",
                                          **ELEMENT)

    # --- Parameters: the nitrogen dependences (Gauthier et al. 2021, the legacy's SurfacicProteins) ---
    Vc_max25_nitrogen_slope: float = _parameter("µmol g-1 s-1", "Vc_max25's increase with Na", 84.965,
                                                unit_comment="of CO2 per g of N", min_value=0.)
    Jmax25_nitrogen_slope: float = _parameter("µmol g-1 s-1", "Jmax25's increase with Na", 117.6,
                                              unit_comment="of electrons per g of N", min_value=0.)
    TPU25_nitrogen_slope: float = _parameter("µmol g-1 s-1", "TPU25's increase with Na", 9.25,
                                             unit_comment="of CO2 per g of N", min_value=0.)
    Rdark25_nitrogen_slope: float = _parameter("µmol g-1 s-1", "Rdark25's increase with Na", 0.493,
                                               unit_comment="of CO2 per g of N", min_value=0.)
    electron_yield_nitrogen_slope: float = _parameter("m2 g-1", "α's increase with Na", 0.0413,
                                                      unit_comment="mol e- mol-1 photon per g N m-2", min_value=0.)
    Vc_max25_nitrogen_minimum: float = _parameter("g m-2", "Na below which Vc_max25 is 0", 0.,
                                                  value_comment="0.17 with TOTAL_NITROGEN", min_value=0.)
    Jmax25_nitrogen_minimum: float = _parameter("g m-2", "Na below which Jmax25 is 0", 0.,
                                                value_comment="0.17 with TOTAL_NITROGEN", min_value=0.)
    TPU25_nitrogen_minimum: float = _parameter("g m-2", "Na below which TPU25 is 0", 0.,
                                               value_comment="0.229 with TOTAL_NITROGEN", min_value=0.)
    Rdark25_nitrogen_minimum: float = _parameter("g m-2", "Na below which Rdark25 is 0", 0.,
                                                 value_comment="0.118 with TOTAL_NITROGEN", min_value=0.)
    electron_yield_intercept: float = _parameter("mol mol-1", "α at Na = 0 (β)", 0.2101 + 0.0083,
                                                 unit_comment="of electrons per photon",
                                                 value_comment="0.2101 with TOTAL_NITROGEN", min_value=0.)
    BWB_slope_factor: float = _parameter("adim", "δ1 of the BWB slope m = δ1 Na^δ2", 14.7, min_value=0.)
    BWB_slope_exponent: float = _parameter("adim", "δ2 of the BWB slope m = δ1 Na^δ2", -0.548)

    # --- Parameters: the temperature dependences (Braune et al. 2009; Bernacchi et al. 2001 for Kc, Ko, Rdark) ---
    Vc_max_activation_enthalpy: float = _parameter("kJ mol-1", "Vc_max's enthalpy of activation ΔHa", 89.7)
    Jmax_activation_enthalpy: float = _parameter("kJ mol-1", "Jmax's enthalpy of activation ΔHa", 48.9)
    TPU_activation_enthalpy: float = _parameter("kJ mol-1", "TPU's enthalpy of activation ΔHa", 47.)
    Kc_activation_enthalpy: float = _parameter("kJ mol-1", "Kc's enthalpy of activation ΔHa", 79.43,
                                               references=BERNACCHI)
    Ko_activation_enthalpy: float = _parameter("kJ mol-1", "Ko's enthalpy of activation ΔHa", 36.38,
                                               references=BERNACCHI)
    Gamma_activation_enthalpy: float = _parameter("kJ mol-1", "Γ*'s enthalpy of activation ΔHa", 35.)
    Rdark_activation_enthalpy: float = _parameter("kJ mol-1", "Rdark's enthalpy of activation ΔHa", 46.39,
                                                  references=BERNACCHI)
    Vc_max_deactivation_enthalpy: float = _parameter("kJ mol-1", "Vc_max's enthalpy of deactivation ΔHd", 149.3)
    Jmax_deactivation_enthalpy: float = _parameter("kJ mol-1", "Jmax's enthalpy of deactivation ΔHd", 152.3)
    TPU_deactivation_enthalpy: float = _parameter("kJ mol-1", "TPU's enthalpy of deactivation ΔHd", 152.3)
    Vc_max_entropy: float = _parameter("kJ mol-1 K-1", "Vc_max's entropy term ΔS", 0.486)
    Jmax_entropy: float = _parameter("kJ mol-1 K-1", "Jmax's entropy term ΔS", 0.495)
    TPU_entropy: float = _parameter("kJ mol-1 K-1", "TPU's entropy term ΔS", 0.495)
    reference_temperature: float = _parameter("K", "The reference temperature of the 25 °C values (Tref)", 298.15)

    # --- Parameters: photosynthesis ---
    O2: float = _parameter("µmol mol-1", "The intercellular O2 concentration", 21000, references=BERNACCHI,
                           min_value=0.)
    KC25: float = _parameter("µmol mol-1", "RuBisCO's Michaelis constant for CO2 at 25 °C (Kc)", 404,
                             references=BERNACCHI, min_value=0.)
    KO25: float = _parameter("µmol mol-1", "RuBisCO's Michaelis constant for O2 at 25 °C (Ko)", 278.4E3,
                             references=BERNACCHI, min_value=0.)
    GAMMA25: float = _parameter("µmol mol-1", "The CO2 compensation point without Rd at 25 °C (Γ*)", 39,
                                min_value=0.)
    THETA: float = _parameter("adim", "The curvature of J's light response (θ)", 0.72, min_value=0., max_value=1.)
    stem_efficiency: float = _parameter("adim", "The relative Ag of a non-lamina element (stems, sheaths)", 0.78,
                                        min_value=0., max_value=1.)
    light_respiration_fraction: float = _parameter("adim", "Rd / Rdark at high light (Rd_A)", 0.33,
                                                   references=MULLER, min_value=0., max_value=1.)
    light_respiration_base: float = _parameter("adim", "The base of Rd's decrease with the PAR (Rd_B)", 0.5,
                                               references=MULLER, min_value=0., max_value=1.)
    light_respiration_scale: float = _parameter("µmol m-2 s-1", "The PAR scale of Rd's decrease (Rd_C)", 15,
                                                references=MULLER, min_value=0.)
    NSC_inhibition_threshold: float = _parameter("µmol m-2", "The surfacic NSC above which they inhibit Ag (WSC_min)",
                                                 100000, references="Azcón-Bieto (1983)", unit_comment="of C",
                                                 min_value=0.)
    NSC_inhibition_maximum: float = _parameter("adim", "The maximal inhibition of Ag by the NSC", 1,
                                               references="Azcón-Bieto (1983)", min_value=0., max_value=1.)
    NSC_inhibition_constant: float = _parameter("µmol m-2", "The NSC's half-inhibition constant (K_Inhibition)",
                                                938000, references="Azcón-Bieto (1983)", unit_comment="of C",
                                                min_value=0.)

    # --- Parameters: stomatal conductance ---
    minimum_stomatal_conductance: float = _parameter("mol m-2 s-1", "gs in the dark (g0, GSMIN)", 0.05, min_value=0.)
    boundary_layer_conductance: float = _parameter("mol m-2 s-1", "The boundary layer's conductance to water vapour "
                                                   "(gb), for Cs and Ci", 3.5, references=MULLER, min_value=0.)
    stomatal_slope: float = _parameter("adim", "The slope m' of the Leuning, Tuzet and hydraulics models", 5.7,
                                       references=ALABASHA, min_value=0.)
    vapour_pressure_deficit_scale: float = _parameter("kPa", "D0 of the Leuning and hydraulics models", 30,
                                                      references=ALABASHA, min_value=0.)
    reference_water_potential: float = _parameter("MPa", "The leaf water potential halving gs (ψref; Tuzet, "
                                                  "hydraulics)", -0.65, references=ALABASHA, max_value=0.)
    water_potential_shape: float = _parameter("adim", "The shape n of gs's response to ψ (Tuzet, hydraulics)", 4,
                                              references=ALABASHA, min_value=0.)
    compensation_point_reference: float = _parameter("mol mol-1", "Γ0, the CO2 compensation point at Tref' of the "
                                                     "Leuning, Tuzet and hydraulics models", 28E-6, references=TUZET,
                                                     min_value=0.)
    compensation_point_linear: float = _parameter("K-1", "Γ1, Γ's linear temperature coefficient", 0.0509,
                                                  references=TUZET)
    compensation_point_quadratic: float = _parameter("K-2", "Γ2, Γ's quadratic temperature coefficient", 0.0010,
                                                     references=TUZET)
    compensation_point_temperature: float = _parameter("K", "Γ's reference temperature Tref'", 293.2,
                                                       references=TUZET)

    # --- Parameters: the energy balance ---
    wind_attenuation: float = _parameter("adim", "The wind's attenuation coefficient in the canopy (a)", 2.5,
                                         references=CAMPBELL, min_value=0.)
    displacement_height_ratio: float = _parameter("adim", "The zero plane displacement height / Zh", 0.7,
                                                  min_value=0., max_value=1.)
    roughness_length_ratio: float = _parameter("adim", "The roughness length / Zh", 0.1, min_value=0., max_value=1.)
    lamina_boundary_layer_factor: float = _parameter("", "rbh = factor · sqrt(w / u) of a lamina", 154,
                                                     unit_comment="s^0.5 m-1",
                                                     references="Monteith (1973); Finnigan and Raupach (1987)",
                                                     min_value=0.)
    cylinder_boundary_layer_factor: float = _parameter("m2 s-1", "rbh = w / (factor · (u w / ν)^c) of a cylinder",
                                                       1.2E-5, min_value=0.)
    cylinder_boundary_layer_viscosity: float = _parameter("m2 s-1", "ν, the air's kinematic viscosity, in a "
                                                          "cylinder's rbh", 1.5E-5, min_value=0.)
    cylinder_boundary_layer_exponent: float = _parameter("adim", "c, the exponent of a cylinder's rbh", 0.47)
    boundary_layer_ratio: float = _parameter("adim", "The boundary layer's resistance to water / to heat", 0.96,
                                             min_value=0.)
    PAR_to_global: float = _parameter("adim", "The absorbed global radiation / the absorbed PAR (in energy)", 1.48,
                                      min_value=0.)
    PAR_photon_energy: float = _parameter("µmol J-1", "The PAR's photons per joule", 4.55,
                                          references="Goudriaan and van Laar (1994)", min_value=0.)

    # --- Settings ---
    stomatal_model: str = setting(description="The stomatal conductance: 'BWB' (Ball, Woodrow and Berry, with Ag), "
                                  "'Leuning', 'Tuzet' (the leaf water potential) or 'hydraulics' (both)",
                                  default="BWB", choices=STOMATAL_MODELS)
    NSC_retroinhibition: bool = setting(description="Ag inhibited by the NSC (min(Ac, Aj)), else min(Ac, Aj, Ap)",
                                        default=True, choices=(True, False))
    primitive_scale: bool = setting(description="The photosynthesis of each primitive of the elements "
                                    "(PARa_primitives, primitive_areas), else of the element (PARa)", default=False,
                                    choices=(True, False))
    convergence_tolerance: float = setting(unit="adim", description="The relative change of Ci and Ts that ends the "
                                           "fixed point", default=0.01)
    maximum_iterations: int = setting(unit="adim", description="The fixed point's iterations at most", default=30)
    reference_height: float = setting(unit="m", description="The height of the wind's measurement above the soil "
                                      "(zr; for wheat, 2 m)", default=2)
    minimum_wind_speed: float = setting(unit="m s-1", description="The smallest wind speed used", default=0.1)
    lamina_labels: tuple = setting(description="The organs' labels of the laminae (a horizontal plane's heat "
                                   "transfer, Ag not reduced)", default=("Blade",))
    exchanging_labels: tuple = setting(description="The labels of the elements that exchange gas",
                                       default=("StemElement", "LeafElement", "LeafElement1"))
    hidden_labels: tuple = setting(description="The labels of the elements inside the sheaths' tube (no exchange, "
                                   "the apex temperature)", default=("HiddenElement",))
    corrections: object = setting(description="The legacy's behaviours corrected: 'recommended' (all), a tuple of "
                                  "names (CORRECTIONS), or () for the legacy's", default="recommended")

    def __post_init__(self):
        super().__post_init__()
        CORRECTIONS.resolve(self.corrections)

    @rate(include_inactive=True)
    def _canopy_height(self, height, green_area, exchanging_labels):
        """Each plant's canopy height: its exchanging elements' highest, with a green area (0 without any)."""
        data_structure = self.data_structure
        exchanging = np.isin(np.asarray(data_structure.labels()), list(exchanging_labels))
        area = np.asarray(green_area, dtype=float)
        values = np.where(exchanging & (area != 0) & np.isfinite(height), height, 0.)
        return np.maximum(data_structure.map_to(values, "node", "Plant", "max"), 0.)

    @rate(include_inactive=True)
    def _capacities(self, surfacic_nitrogen, Vc_max25_nitrogen_slope, Jmax25_nitrogen_slope, TPU25_nitrogen_slope,
                    Rdark25_nitrogen_slope, electron_yield_nitrogen_slope, Vc_max25_nitrogen_minimum,
                    Jmax25_nitrogen_minimum, TPU25_nitrogen_minimum, Rdark25_nitrogen_minimum,
                    electron_yield_intercept) -> dict:
        """The photosynthetic capacities at 25 °C, linear in the surfacic nitrogen."""
        nitrogen = np.asarray(surfacic_nitrogen, dtype=float)
        return {"Vc_max25": Vc_max25_nitrogen_slope * (nitrogen - Vc_max25_nitrogen_minimum),
                "Jmax25": Jmax25_nitrogen_slope * (nitrogen - Jmax25_nitrogen_minimum),
                "TPU25": TPU25_nitrogen_slope * (nitrogen - TPU25_nitrogen_minimum),
                "Rdark25": Rdark25_nitrogen_slope * (nitrogen - Rdark25_nitrogen_minimum),
                "electron_yield": electron_yield_nitrogen_slope * nitrogen + electron_yield_intercept}

    @state(include_inactive=True)
    def _exchange(self, PARa, green_area, height, width, diameter, surfacic_nitrogen,
                  surfacic_NSC, leaf_water_potential, SAM_temperature, PARa_primitives, primitive_areas, canopy_height,
                  Vc_max25, Jmax25, TPU25, Rdark25, electron_yield, Ci, BWB_slope_factor, BWB_slope_exponent,
                  Vc_max_activation_enthalpy, Jmax_activation_enthalpy, TPU_activation_enthalpy, Kc_activation_enthalpy,
                  Ko_activation_enthalpy, Gamma_activation_enthalpy, Rdark_activation_enthalpy,
                  Vc_max_deactivation_enthalpy, Jmax_deactivation_enthalpy, TPU_deactivation_enthalpy, Vc_max_entropy,
                  Jmax_entropy, TPU_entropy, reference_temperature, O2, KC25, KO25, GAMMA25, THETA, stem_efficiency,
                  light_respiration_fraction, light_respiration_base, light_respiration_scale, NSC_inhibition_threshold,
                  NSC_inhibition_maximum, NSC_inhibition_constant, minimum_stomatal_conductance,
                  boundary_layer_conductance, stomatal_slope, vapour_pressure_deficit_scale, reference_water_potential,
                  water_potential_shape, compensation_point_reference, compensation_point_linear,
                  compensation_point_quadratic, compensation_point_temperature, wind_attenuation,
                  displacement_height_ratio, roughness_length_ratio, lamina_boundary_layer_factor,
                  cylinder_boundary_layer_factor, cylinder_boundary_layer_viscosity, cylinder_boundary_layer_exponent,
                  boundary_layer_ratio, PAR_to_global, PAR_photon_energy) -> dict:
        """Each exchanging element's fixed point (element_exchange), at the weather of the step's end, over its green
        area; the others exchange nothing."""
        arguments = locals()
        data_structure = self.data_structure
        n = np.asarray(PARa).size
        labels = np.asarray(data_structure.labels())
        exchanging = np.isin(labels, list(self.exchanging_labels))
        hidden = np.isin(labels, list(self.hidden_labels))
        lamina = np.isin(np.asarray(data_structure.labels(location="Organ")), list(self.lamina_labels))
        per_node = {name: np.broadcast_to(np.asarray(arguments[name], dtype=float), (n,))
                    for name in ("PARa", "height", "width", "surfacic_nitrogen", "surfacic_NSC",
                                 "leaf_water_potential", "Vc_max25", "Jmax25", "TPU25", "Rdark25", "electron_yield",
                                 "Ci")}
        organ = np.asarray(data_structure.owner("Organ"))
        diameter = np.asarray(diameter, dtype=float)
        stem_diameter = np.where(organ >= 0, diameter[organ.clip(0)] if diameter.size else 0., 0.)
        dimension = np.where(lamina, per_node["width"], stem_diameter)
        plant_canopy = np.asarray(canopy_height, dtype=float)[np.asarray(data_structure.owner("Plant"))]
        apex = np.asarray(SAM_temperature, dtype=float)
        if apex.ndim:
            apex = apex[np.asarray(data_structure.owner("Axis"))]

        weather = Weather(*(float(self.forcing(name)) for name in Weather._fields))
        options = Options(self.stomatal_model, bool(self.NSC_retroinhibition), float(self.convergence_tolerance),
                          int(self.maximum_iterations), float(self.reference_height), float(self.minimum_wind_speed),
                          CORRECTIONS.resolve(self.corrections))
        apex = np.broadcast_to(np.where(np.isnan(apex), weather.air_temperature, apex), (n,))

        outputs = {name: np.zeros(n) for name in ("Ag", "An", "Rd", "gs", "Ci", "Tr")}
        outputs["Ts"] = np.where(exchanging | hidden, apex, np.nan)
        area = np.broadcast_to(np.asarray(green_area, dtype=float), (n,))
        computing = exchanging & (area > 0) & np.isfinite(per_node["height"]) & (plant_canopy > 0)
        computing &= per_node["surfacic_nitrogen"] != 0
        indices = np.flatnonzero(computing)
        if self.primitive_scale and indices.size and (PARa_primitives is None or primitive_areas is None):
            raise ValueError("GasExchange(primitive_scale=True): the elements' PARa_primitives and primitive_areas "
                             "are not given")
        parameter_names = ElementParameters._fields
        table = np.column_stack([np.broadcast_to(np.asarray(arguments[name], dtype=float), (n,))[indices]
                                 for name in parameter_names]) if indices.size else np.zeros((0, len(parameter_names)))
        kinds, kind_of = np.unique(table, axis=0, return_inverse=True)
        plant_parameters = [ElementParameters(*row.tolist()) for row in kinds]
        not_converged = 0
        for i, kind in zip(indices.tolist(), np.ravel(kind_of).tolist()):
            previous_Ci = float(per_node["Ci"][i])
            start_Ci = previous_Ci if previous_Ci > 0 else weather.ambient_CO2       # NaN or 0: the ambient CO2
            if self.primitive_scale:
                PARs = [float(value) for value in PARa_primitives[i]]
                areas = [float(value) for value in primitive_areas[i]]
            else:
                PARs, areas = [float(per_node["PARa"][i])], [1.0]
            capacities = tuple(float(per_node[name][i]) for name in ("Vc_max25", "Jmax25", "TPU25", "Rdark25",
                                                                      "electron_yield"))
            result = primitives_exchange(PARs, areas, float(per_node["surfacic_nitrogen"][i]),
                                         float(per_node["surfacic_NSC"][i]), capacities, float(dimension[i]),
                                         float(per_node["height"][i]), float(plant_canopy[i]), bool(lamina[i]),
                                         float(per_node["leaf_water_potential"][i]), start_Ci, weather,
                                         plant_parameters[kind], options)
            for name, value in zip(("Ag", "An", "Rd", "gs", "Ci", "Ts", "Tr"), result[:7]):
                outputs[name][i] = value
            not_converged += not result[7]
        if not_converged:
            warnings.warn(f"GasExchange: {not_converged} elements did not converge in {self.maximum_iterations} "
                          f"iterations (Ci, Ts within {self.convergence_tolerance})", RuntimeWarning)
        outputs["exchanging_area"] = np.where(computing, area, 0.)
        outputs["transpiration"] = outputs["Tr"] * outputs["exchanging_area"]
        return outputs


# The per-plant parameters of one element's equations, by their names in GasExchange (the nitrogen dependences are
# its capacities, computed before)
_CAPACITY_PARAMETERS = {"Vc_max25_nitrogen_slope", "Jmax25_nitrogen_slope", "TPU25_nitrogen_slope",
                        "Rdark25_nitrogen_slope", "electron_yield_nitrogen_slope", "Vc_max25_nitrogen_minimum",
                        "Jmax25_nitrogen_minimum", "TPU25_nitrogen_minimum", "Rdark25_nitrogen_minimum",
                        "electron_yield_intercept"}
ElementParameters = namedtuple("ElementParameters", [
    item.name for item in fields(GasExchange)
    if item.metadata.get("role") == "parameter" and item.name not in _CAPACITY_PARAMETERS])


def default_parameters(**changes) -> ElementParameters:
    """GasExchange's default parameters for element_exchange (outside a component), *changes* given."""
    values = {item.name: item.default for item in fields(GasExchange) if item.name in ElementParameters._fields}
    values.update(changes)
    return ElementParameters(**values)


def default_capacities() -> dict:
    """GasExchange's default parameters of the capacities' nitrogen dependences ({name: value})."""
    return {item.name: item.default for item in fields(GasExchange) if item.name in _CAPACITY_PARAMETERS}
