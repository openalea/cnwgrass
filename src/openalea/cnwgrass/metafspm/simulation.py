"""
CNW-Grass's gas exchange for users: a small wheat stand's photosynthesis and transpiration in a few lines, its
results as tables.

    from openalea.cnwgrass.metafspm.simulation import simulate, parameters

    result = simulate()                         a stand of 12 wheat plants (Adel's), lit by Caribu over a clear May
                                                day, hour by hour: each leaf's and stem's gas exchange
    result.summary                              per hour: the weather, the PAR the stand absorbs, its assimilation,
                                                its transpiration, its leaves' mean temperature and conductance
    result.elements()                           per leaf and stem element at the last hour: its plant, organ, height,
                                                green area, absorbed PAR and gas exchange
    result.plants()                             per plant: its green area, the CO2 it fixed and the water it lost
    result.draw()                               the stand in 3D, coloured by a variable (An by default)
    parameters()                                every parameter and setting: its value, unit and description
    simulate({"stomatal_model": "Leuning"})     another stomatal conductance model
    simulate(air_temperature=series)            a measured weather (each a number, a function of the hour, or a
                                                pandas Series indexed by hour)
    simulate(outputs="outputs")                 the run's files too: run.json, scenario.json, summary.csv,
                                                plants.csv, the plants' summaries (Plant/)

The stand: rows of wheat plants grown by Adel's port (ScenarioShoot: ShootArchitecture grown by ScenarioDevelopment,
CN-Wheat's leaf shapes, Adel's LeafShapes.cnwheat()), their bases at the soil surface, with development delays of up
to 60 °Cd, grown to *thermal_time* before the day; then, every hour, Caribu (CaribuLightModel, with the sun) lights it
and GasExchange computes each element's gas exchange from the light it absorbs (PARa, linked by name) and the
weather. The plants keep growing during the day, by the degree-days of the air temperature. Lengths in m, PAR in
µmol m-2 s-1.
"""
from dataclasses import dataclass, fields

import numpy as np
import pandas as pd

from openalea.adel.architecture import ShootArchitecture
from openalea.adel.development import ScenarioDevelopment, ScenarioShoot, random_delays, sown_stand
from openalea.adel.leaf_shapes import LeafShapes
from openalea.caribu.component import CaribuLightModel
from openalea.caribu.optics import PlantOptics
from openalea.caribu.simulation import clear_day_PAR
from openalea.metafspm.coupling.forcings import as_forcing, forcing_mean
from openalea.metafspm.scene.population import planting
from openalea.metafspm.scene.recorder import SceneRecorder
from openalea.metafspm.scene.scenarios import parameter_table
from openalea.metafspm.scene.scene import Scene

from .gas_exchange import GasExchange, adel_links

HOUR = 3600.
LATITUDE = 48.85                            # Paris, Caribu's default
EXCHANGING = ("StemElement", "LeafElement")

# The parameters a beginner changes first: the stomatal model, the photosynthetic capacities' nitrogen dependences,
# the light response's curvature, the stomata's opening, the stems' efficiency and the wind in the canopy
MAIN_PARAMETERS = ("stomatal_model", "Vc_max25_nitrogen_slope", "Jmax25_nitrogen_slope", "THETA", "BWB_slope_factor",
                   "minimum_stomatal_conductance", "stomatal_slope", "stem_efficiency", "wind_attenuation",
                   "corrections")

# The weather of the default day: a mild clear May day (°C, a fraction, µmol mol-1, m s-1)
WEATHER = dict(air_temperature=20., relative_humidity=0.6, ambient_CO2=410., wind_speed=2.)


class ExchangingWheat(ScenarioShoot):
    """The plant model: Adel's wheat shoot (ShootArchitecture grown by ScenarioDevelopment), its elements' optics
    (Caribu's PlantOptics, which the light reads) and their gas exchange (GasExchange: its PARa the light's, by name;
    its elements' height from Adel's z1, z2, adel_links). A subclass may compute the gas exchange with a subclass of
    GasExchange (gas_exchange_component), its links named after it (translator, adel_links(receiver=...))."""
    gas_exchange_component = GasExchange
    translator = adel_links()

    def __init__(self, data_structure, time_step):
        self.architecture = ShootArchitecture(data_structure=data_structure, time_step=time_step)
        self.development = ScenarioDevelopment(data_structure=data_structure, time_step=time_step)
        self.optics = PlantOptics(data_structure=data_structure, time_step=time_step)
        self.gas_exchange = self.gas_exchange_component(data_structure=data_structure, time_step=time_step)
        self.couple(self.development, self.architecture, self.optics, self.gas_exchange, translator=self.translator)
        if not data_structure.has("leaf_water_potential"):   # an optional input, held here to be given (simulate)
            data_structure.register("leaf_water_potential", GasExchange.leaf_water_potential, location="node")

    def grow(self):
        """The shoot grown to its stage, its development's start (minus each plant's delay), by one call before the
        scene's first step (ScenarioShoot.run_shoot: every phytomer the stage asks for; the scene's first step adds
        no degree-days, none passed since)."""
        self.run_shoot()
        self.optics.run_stage("state")                  # the optics before the first lighting

    def run(self):
        self.run_shoot()
        self.optics()                                   # the new elements' too
        self.gas_exchange()


def parameters(values: dict = None, main_only: bool = False) -> pd.DataFrame:
    """GasExchange's parameters and settings, one row each: its name, its value (the default, or the one in
    *values*), its unit, its description, whether it is one of the main ones (those first), and its bounds
    (min_value, max_value)."""
    table = parameter_table(GasExchange, values=values, main=MAIN_PARAMETERS).drop(columns="component")
    return table[table["main"]].reset_index(drop=True) if main_only else table


def _exchanging(ds) -> np.ndarray:
    """The elements that may exchange gas: Adel's stem and leaf elements with a green area."""
    return np.isin(np.asarray(ds.labels()), EXCHANGING) & (np.asarray(ds.get("green_area")) > 0)


class _Summary:
    """The summary row of the hour just computed; none before the first step."""

    def __init__(self, forcings: dict, ground: float):
        self.forcings, self.ground = forcings, ground

    def __call__(self, scene):
        if scene.iteration == 0:
            return None
        ds = scene.populations[0].data_structure
        light = scene.environment[0].light.data_structure
        weather = {name: forcing_mean(self.forcings[name], scene.time - scene.time_step, scene.time)  # the hour's
                   for name in ("PARi", *WEATHER)}
        exchanging_area = np.asarray(ds.get("exchanging_area"))
        exchanging = _exchanging(ds)
        area = exchanging_area[exchanging]

        def per_ground(name):
            return float(np.sum(np.asarray(ds.get(name)) * exchanging_area)) / self.ground

        def mean(name):
            return float(np.sum(np.asarray(ds.get(name))[exchanging] * area) / np.sum(area)) if area.sum() > 0 \
                else np.nan

        return dict(hour=scene.start[1] + scene.time / HOUR, **weather,
                    PAR_absorbed=float(np.sum(np.asarray(light.get("PARa")) * np.asarray(light.get("lit_area"))))
                    / self.ground,
                    gross_assimilation=per_ground("Ag"), net_assimilation=per_ground("An"),
                    transpiration=float(np.sum(ds.get("transpiration"))) / self.ground,
                    surface_temperature=mean("Ts"), stomatal_conductance=mean("gs"))


def _plants_row(scene):
    """Per plant, the hour just computed: its green area (m2), net assimilation (µmol CO2 s-1) and transpiration
    (mmol H2O s-1)."""
    if scene.iteration == 0:
        return None
    ds = scene.populations[0].data_structure
    green_area, exchanging_area = np.asarray(ds.get("green_area")), np.asarray(ds.get("exchanging_area"))

    def per_plant(values):
        return np.asarray(ds.map_to(np.asarray(values, dtype=float), "node", "Plant", "sum"))

    return pd.DataFrame(dict(hour=scene.start[1] + scene.time / HOUR, plant=np.arange(len(ds.entity_ids("Plant"))),
                             green_area=per_plant(green_area),
                             net_assimilation=per_plant(np.asarray(ds.get("An")) * exchanging_area),
                             transpiration=per_plant(ds.get("transpiration"))))


@dataclass
class Simulation:
    """A stand's day of gas exchange: its *summary* (a DataFrame, one row per hour), its plants' hours
    (*plant_hours*, one row per plant and hour) and its *scene* (metafspm's: the plants' DataStructure holds every
    element's variables)."""
    scene: object
    summary: pd.DataFrame
    plant_hours: pd.DataFrame
    ground_area: float

    @property
    def data_structure(self):
        return self.scene.populations[0].data_structure

    @property
    def architecture(self):
        return self.scene.populations[0].instance.architecture

    def elements(self) -> pd.DataFrame:
        """The leaf and stem elements with a green area at the last hour, one row each: its plant (0, 1, ...), organ
        (blade, sheath or internode), height (m), green area (m2), exchanging area (m2: its green area above the
        soil, 0 if it exchanges nothing), PARa (the PAR it absorbs, µmol m-2 s-1), Ag and An
        (its gross and net assimilation, µmol CO2 m-2 s-1), gs (its stomatal conductance, mol m-2 s-1), Ci (its
        internal CO2, µmol mol-1), Ts (its temperature, °C), Tr (its transpiration rate, mmol m-2 s-1) and
        transpiration (Tr times its exchanging area, mmol s-1)."""
        ds = self.data_structure
        exchanging = _exchanging(ds)
        organ = self.architecture.organ_names()[np.asarray(ds.owner("Organ"))]
        table = pd.DataFrame(dict(plant=np.asarray(ds.owner("Plant")), organ=organ,
                                  **{name: np.asarray(ds.get(name), dtype=float) for name in (
                                      "height", "green_area", "exchanging_area", "PARa", "Ag", "An", "gs", "Ci", "Ts",
                                      "Tr", "transpiration")}))
        return table[exchanging].reset_index(drop=True)

    def plants(self) -> pd.DataFrame:
        """Per plant: its green area at the last hour (m2), the CO2 it fixed (net, mmol) and the water it transpired
        (g) over the run."""
        hours = self.plant_hours.groupby("plant")
        return pd.DataFrame({"green_area": hours["green_area"].last(),
                             "net_assimilation": hours["net_assimilation"].sum() * HOUR / 1e3,
                             "water": hours["transpiration"].sum() * HOUR * 18e-3})

    def draw(self, ax=None, color_by: str = "An"):
        """The stand in 3D, its elements coloured by *color_by* (An by default; a variable of the elements, e.g.
        "PARa", "Ts", "Tr"); returns the matplotlib axes."""
        import matplotlib.pyplot as plt
        from openalea.metafspm.display import draw, value_colors

        if ax is None:
            ax = plt.figure(figsize=(8, 6)).add_subplot(projection="3d")
        values = np.nan_to_num(np.asarray(self.data_structure.get(color_by), dtype=float))
        limits = (float(values.min()), max(float(values.max()), float(values.min()) + 1e-12))
        draw(ax, self.architecture.triangles(), colors=value_colors(values, limits=limits))
        ax.figure.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*limits), cmap="viridis"), ax=ax, shrink=0.6,
                           pad=0.12, label=color_by)
        return ax


def _check(values: dict) -> dict:
    names = {item.name for item in fields(GasExchange)}
    declared = set(GasExchange.declared_parameters()) | set(GasExchange.declared_settings())
    for name in values:
        if name not in declared or name not in names:
            raise KeyError(f"simulate: '{name}' is no parameter or setting of GasExchange (parameters() lists them)")
    return values


def simulate(parameters: dict = None, hours: float = 24., start_hour: float = 0., day_of_year: float = 130.,
             thermal_time: float = 600., PARi=None, air_temperature=None, relative_humidity=None, ambient_CO2=None,
             wind_speed=None, surfacic_nitrogen=2., leaf_water_potential=None, plants: int = 12, rows: int = 3,
             row_spacing: float = 0.15, plant_spacing: float = 0.04, outputs: str = None,
             compute=None) -> Simulation:
    """The gas exchange of a stand of *plants* wheat plants in *rows* rows (*row_spacing*, *plant_spacing* m apart;
    repeated around it: an infinite canopy), grown to *thermal_time* °Cd since sowing, hour by hour for *hours* hours
    from *start_hour* (universal time) of day *day_of_year*, at Paris' latitude. *parameters* ({name: value}) changes
    GasExchange's parameters and settings: parameters() lists them all.

    The weather (each a number, a function of the hour of the day, or a pandas Series indexed by it): *PARi* (µmol
    m-2 s-1, the PAR above the stand; a clear day's by default, Caribu's clear_day_PAR), *air_temperature* (°C, 20 by
    default), *relative_humidity* (a fraction, 0.6), *ambient_CO2* (µmol mol-1, 410), *wind_speed* (m s-1 at 2 m,
    2). *surfacic_nitrogen* (g m-2; 2, CNW-Grass's default) is every element's nitrogen per unit area, given the same
    way (a CN model would compute it). *leaf_water_potential* (MPa) is every element's leaf water potential, which the
    stomatal models "Tuzet" and "hydraulics" read (a hydraulic model, MECHA, would compute it; by default
    GasExchange's, -0.1).

    The clock: a step reads its forcings' means over its hour (metafspm's forcings over the step), its sun at its end:
    the scene starts an hour before *start_hour*, so that its first step is the hour ending at *start_hour*, the
    summary's first hour (its weather the hour's means).

    *outputs*: a folder for the run's files (metafspm's SceneRecorder): run.json, scenario.json, summary.csv,
    plants.csv (the plants' hours) and Plant/ExchangingWheat/summaries.csv (every variable per plant, every hour).
    *compute*: the Scene's run configuration (metafspm's Compute: its threads, its workers the light's projection
    threads; None: Compute())."""
    scenario = _check(dict(parameters or {}))
    shoot = {"n_sect": 2, "leaves": LeafShapes.cnwheat(), "start": float(thermal_time)}
    stand = sown_stand(plants, n_rows=rows, row_spacing=row_spacing, plant_spacing=plant_spacing)
    xrange, yrange = int(np.ceil(plants / rows)) * plant_spacing, rows * row_spacing
    table = planting(stand, ExchangingWheat, {**shoot, **scenario},
                     [{"delay": float(delay)} for delay in random_delays(plants)], xrange=xrange, yrange=yrange)
    given = dict(PARi=PARi if PARi is not None else clear_day_PAR(day_of_year, LATITUDE),
                 air_temperature=air_temperature, relative_humidity=relative_humidity, ambient_CO2=ambient_CO2,
                 wind_speed=wind_speed, surfacic_nitrogen=surfacic_nitrogen,
                 leaf_water_potential=leaf_water_potential)
    forcings = {name: as_forcing(value if value is not None else WEATHER[name], time_unit=HOUR, origin="calendar")
                for name, value in given.items() if value is not None or name in WEATHER}
    ground = xrange * yrange
    recorder = SceneRecorder(outputs, summaries=outputs is not None,
                             tables={"summary": _Summary(forcings, ground), "plants": _plants_row})
    built = Scene(table, environment=[CaribuLightModel],
                  environment_scenarios=[{"sun": True, "latitude": LATITUDE, "interval": HOUR}], time_step=HOUR,
                  forcings=forcings, recorder=recorder, start=(day_of_year, start_hour - 1.), compute=compute)
    built.populations[0].instance.grow()
    built.simulate(int(round(hours)))
    return Simulation(scene=built, summary=recorder.table("summary"), plant_hours=recorder.table("plants"),
                      ground_area=ground)
