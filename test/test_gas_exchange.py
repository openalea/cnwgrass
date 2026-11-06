"""GasExchange (openalea.cnwgrass.metafspm.gas_exchange) against the legacy openalea.cnwgrass.gasexchange Simulation:
its equations element by element, bit for bit with the legacy's numerics (corrections=()), for every stomatal model
and option of the legacy; the component on Adel's plants (ShootArchitecture) in a metafspm Scene, bit for bit on
every element it computes, its parameters per plant; the corrections and the fixed failures."""
import os
from dataclasses import fields

import numpy as np
import pandas as pd
import pytest

from openalea.cnwgrass.gasexchange import converter, parameters as legacy_parameters, simulation
from openalea.cnwgrass.metafspm import gas_exchange as port
from openalea.metafspm.testing.comparison import assert_arrays_equal

OUTPUTS = ("Ag", "An", "Rd", "gs", "Ci", "Ts", "Tr")
WEATHERS = [(18.8, 360., 0.53, 2.2), (5., 400., 0.9, 0.05), (30., 420., 0.2, 5.)]
CANOPY_HEIGHT = 0.85
LEGACY = frozenset()
STANDALONE = os.path.join(os.path.dirname(__file__), os.pardir, "example", "Farquhar_standalone")


def python_floats(element):
    """A legacy element's inputs as Python floats: with numpy's, the legacy's x ** y is numpy's power, which may
    differ from the C library's pow by an ulp (the port computes in Python floats)."""
    return {name: [float(item) for item in value] if isinstance(value, list) else float(value)
            for name, value in element.items()}


def random_elements(n, seed=1, primitives=False):
    """Legacy element inputs {(plant, axis, metamer, organ, element): inputs}: blades and sheaths, in the dark or
    lit, with NSC or none, a previous Ci or none."""
    generator = np.random.default_rng(seed)
    elements = {}
    for index in range(n):
        blade = bool(generator.random() < 0.6)
        green_area = 10 ** generator.uniform(-5, -2.5)
        count = int(generator.integers(1, 5))
        elements[(1, "MS", index, "blade" if blade else "sheath", "LeafElement1" if blade else "StemElement")] = dict(
            width=generator.uniform(0.001, 0.02), height=generator.uniform(0.01, 0.8),
            PARa=generator.choice([0., generator.uniform(0, 2000)]), nitrates=generator.uniform(0, 5),
            amino_acids=generator.uniform(0, 20), proteins=generator.uniform(0, 800) * green_area / 3e-3,
            Nstruct=generator.uniform(0, 1e-3), green_area=green_area,
            sucrose=generator.uniform(0, 3e6) * green_area * generator.choice([0, 1]), starch=0., fructan=0.,
            PARa_prim=list(generator.uniform(0, 2000, count)) if primitives else [],
            area_prim=list(generator.uniform(1e-5, 1e-4, count)) if primitives else 0,
            Ci=generator.choice([np.nan, generator.uniform(100, 400)]), water_potential=generator.uniform(-2, -0.05))
    return {key: python_floats(element) for key, element in elements.items()}


def legacy_run(elements, weather, model="BWB", canopy_height=CANOPY_HEIGHT):
    """The legacy Simulation's outputs, its global parameters as they are."""
    legacy = simulation.Simulation(stomatal_model_name=model, hydraulics=model in ("Tuzet", "hydraulics"))
    legacy.initialize({"elements": elements, "axes": {(1, "MS"): {"SAM_temperature": weather[0],
                                                                  "height_canopy": canopy_height}}})
    legacy.run(*weather)
    return {key: {name: float(legacy.outputs[key][name]) for name in OUTPUTS} for key in elements}


def capacities(nitrogen, total=False):
    """The capacities at 25 °C of the default parameters (TOTAL_NITROGEN's with *total*), as _capacities."""
    values = dict(port.default_capacities(), **(port.TOTAL_NITROGEN if total else {}))
    return tuple(values[f"{name}_nitrogen_slope"] * (nitrogen - values[f"{name}_nitrogen_minimum"])
                 for name in ("Vc_max25", "Jmax25", "TPU25", "Rdark25")) \
        + (values["electron_yield_nitrogen_slope"] * nitrogen + values["electron_yield_intercept"],)


def port_run(elements, weather, model="BWB", retroinhibition=True, total=False, primitive=False,
             corrections=LEGACY, canopy_height=CANOPY_HEIGHT, parameters=None):
    """The ported equations' outputs for the legacy's elements (its surfacic contents by the port's functions)."""
    options = port.Options(model, retroinhibition, 0.01, 30, 2., 0.1, frozenset(corrections))
    weather = port.Weather(*weather)
    parameters = parameters or port.default_parameters()
    outputs = {}
    for key, element in elements.items():
        if total:
            nitrogen = port.surfacic_total_nitrogen(element["nitrates"], element["amino_acids"], element["proteins"],
                                                    element["Nstruct"], element["green_area"])
        else:
            nitrogen = port.surfacic_photosynthetic_nitrogen(element["proteins"], element["green_area"])
        NSC = port.surfacic_nonstructural_carbohydrates(element["sucrose"], element["starch"], element["fructan"],
                                                        element["green_area"])
        PARs, areas = (element["PARa_prim"], element["area_prim"]) if primitive else ([element["PARa"]], [1.0])
        Ci = element["Ci"] if element["Ci"] > 0 else weather.ambient_CO2
        result = port.primitives_exchange(PARs, areas, nitrogen, NSC, capacities(nitrogen, total), element["width"],
                                          element["height"], canopy_height, key[3] == "blade",
                                          element["water_potential"], Ci, weather, parameters, options)
        outputs[key] = dict(zip(OUTPUTS, result[:7]))
    return outputs


def assert_same(expected, actual, keys=None):
    keys = list(expected) if keys is None else keys
    for name in OUTPUTS:
        assert_arrays_equal([expected[key][name] for key in keys], [actual[key][name] for key in keys], name=name)


@pytest.fixture
def quiet(capsys):
    """The legacy prints the elements that do not converge."""
    yield
    capsys.readouterr()


@pytest.mark.parametrize("model", port.STOMATAL_MODELS)
@pytest.mark.parametrize("weather", WEATHERS)
def test_the_equations_reproduce_the_legacy(model, weather, quiet):
    elements = random_elements(400)
    assert_same(legacy_run(elements, weather, model), port_run(elements, weather, model))


@pytest.mark.parametrize("retroinhibition, total", [(False, False), (True, True), (False, True)])
def test_the_legacy_versions_reproduced(monkeypatch, retroinhibition, total, quiet):
    """NSC_Retroinhibition False (Ag = min(Ac, Aj, Ap)); SurfacicProteins False (the total nitrogen, TOTAL_NITROGEN)."""
    monkeypatch.setattr(legacy_parameters, "NSC_Retroinhibition", retroinhibition)
    if total:
        monkeypatch.setattr(legacy_parameters, "SurfacicProteins", False)
        monkeypatch.setattr(legacy_parameters, "PARAM_N", dict(
            legacy_parameters.PARAM_N, beta=0.2101,
            surfacic_nitrogen_min={"Vc_max25": 0.17, "Jmax25": 0.17, "TPU25": 0.229, "Rdark25": 0.118}))
    elements = random_elements(300, seed=2)
    for model in ("BWB", "hydraulics"):
        assert_same(legacy_run(elements, WEATHERS[0], model),
                    port_run(elements, WEATHERS[0], model, retroinhibition=retroinhibition, total=total))


def test_the_primitive_scale(monkeypatch, quiet):
    """The legacy's prim_scale: Ag the area-weighted mean of the primitives', the other outputs the last one's
    (reproduced with corrections=()); with primitive_aggregation all are area-weighted means."""
    monkeypatch.setattr(legacy_parameters, "prim_scale", True)
    elements = random_elements(200, seed=3, primitives=True)
    assert_same(legacy_run(elements, WEATHERS[0]), port_run(elements, WEATHERS[0], primitive=True))
    corrected = port_run(elements, WEATHERS[0], primitive=True, corrections={"primitive_aggregation"})
    legacy = port_run(elements, WEATHERS[0], primitive=True)
    assert all(corrected[key]["Ag"] == legacy[key]["Ag"] for key in elements)
    key = next(key for key, element in elements.items() if len(element["PARa_prim"]) > 2)
    single = [port_run({key: dict(elements[key], PARa_prim=elements[key]["PARa_prim"][:k + 1],
                                  area_prim=elements[key]["area_prim"][:k + 1])}, WEATHERS[0], primitive=True)[key]
              for k in range(len(elements[key]["PARa_prim"]))]
    areas = elements[key]["area_prim"]
    for name in ("An", "Rd", "Tr", "gs"):        # each primitive's own, started from the one before
        assert corrected[key][name] == pytest.approx(sum(row[name] * area for row, area in zip(single, areas))
                                                     / sum(areas), rel=1e-12)


def test_the_standalone_example(quiet):
    """The legacy's example (example/Farquhar_standalone), its HiddenElement computed as it has a height."""
    elements = converter.from_dataframe(pd.read_csv(os.path.join(STANDALONE, "elements_inputs.csv")),
                                        pd.read_csv(os.path.join(STANDALONE, "axes_inputs.csv")))["elements"]
    for element in elements.values():
        element.update(PARa_prim=[], area_prim=0)
    weather = (18.8, 360, 0.53, 2.2)
    for model in ("BWB", "hydraulics"):
        assert_same(legacy_run(elements, weather, model, canopy_height=0.7),
                    port_run(elements, weather, model, canopy_height=0.7))


@pytest.mark.parametrize("model", port.STOMATAL_MODELS)
def test_the_corrections(model, quiet):
    """Recommended: the slope of the vapour pressure moves the fixed point within its tolerance; the compensation
    point in µmol mol-1 lowers Cs − Γ by about 27 µmol mol-1, raising the Leuning and Tuzet conductances."""
    elements = random_elements(200, seed=4)
    legacy = port_run(elements, WEATHERS[0], model)
    corrected = port_run(elements, WEATHERS[0], model, corrections=port.CORRECTIONS.resolve("recommended"))
    lit = [key for key in elements if legacy[key]["Ag"] > 1]
    ratio = np.array([corrected[key]["gs"] / legacy[key]["gs"] for key in lit])
    if model == "BWB":
        assert np.abs(ratio - 1).max() < 0.05
    else:
        assert np.median(ratio) > 1.02
    assert all(np.isfinite(corrected[key][name]) for key in elements for name in OUTPUTS)


def test_the_legacy_failures_fixed():
    """An air at 0 °C with Ci converged at once (in the dark: the legacy divides Ts's change by 0), a cylinder without
    diameter (0 / 0): the port computes them (a previous Ci of 0: test_the_component_restarts_Ci)."""
    options = port.Options("BWB", True, 0.01, 30, 2., 0.1, port.CORRECTIONS.resolve("recommended"))
    parameters = port.default_parameters()
    for PAR, temperature, dimension in ((0., 0., 0.003), (800., 0., 0.003), (800., 15., 0.)):
        weather = port.Weather(temperature, 360., 0.6, 2.)
        result = port.element_exchange(PAR, 2., 0., capacities(2.), dimension, 0.4, 0.8, False, -0.1, 360.,
                                       temperature, weather, parameters, options)
        assert all(np.isfinite(result[:7])) and (result[0] > 0) == (PAR > 0) and result[6] > 0
    key, element = next(iter(random_elements(1).items()))
    with pytest.raises(ZeroDivisionError):
        legacy_run({key: dict(element, PARa=0., Ci=np.nan)}, (0., 360., 0.6, 2.))


# ── the component, on Adel's plants ──────────────────────────────────────────

adel = pytest.importorskip("openalea.adel.development")
from openalea.metafspm.coupling.choregrapher import Choregrapher  # noqa: E402
from openalea.metafspm.scene.population import Planting, planting  # noqa: E402
from openalea.metafspm.scene.scene import Scene  # noqa: E402
from openalea.metafspm.coupling.translator import Translator  # noqa: E402


class ExchangingShoot(adel.ScenarioShoot):
    """Adel's plant model with its gas exchange (its elements' height from their z1, z2: adel_links)."""
    translator = port.adel_links()

    def __init__(self, data_structure, time_step):
        super().__init__(data_structure, time_step)
        self.gas_exchange = port.GasExchange(data_structure=data_structure, time_step=time_step)
        self.couple(self.development, self.architecture, self.gas_exchange, translator=self.translator)

    def run(self):
        super().run()
        self.gas_exchange()


@pytest.fixture(scope="module")
def shoots():
    """Two plants grown 60 days, the second with its own THETA and stem_efficiency; the weather a mutable dict."""
    Choregrapher().reset()
    weather = dict(air_temperature=18.8, relative_humidity=0.53, ambient_CO2=360., wind_speed=2.2)
    table = planting([Planting(), Planting(x=0.2)], ExchangingShoot, {},
                     [{}, {"THETA": 0.8, "stem_efficiency": 0.6}])
    scene = Scene(table, time_step=86400., forcings={name: (lambda time, name=name: weather[name])
                                                     for name in weather})
    for _ in range(60):
        scene.run()
    return scene, weather


def _set_inputs(ds, seed):
    """Random inputs on every element (the heights through Adel's z1 = z2; the surfacic nitrogen from proteins per
    unit green area, as the legacy computes it), the previous Ci too."""
    generator = np.random.default_rng(seed)
    n = ds.n_nodes()
    height = generator.uniform(0.01, 0.6, n)
    proteins = generator.uniform(0, 2e5, n) * (generator.random(n) < 0.95)
    values = dict(PARa=generator.uniform(0, 2000, n) * (generator.random(n) < 0.8),
                  green_area=10 ** generator.uniform(-5, -2.5, n), width=generator.uniform(0.001, 0.02, n),
                  surfacic_nitrogen=port.surfacic_photosynthetic_nitrogen(proteins, 1.),
                  surfacic_NSC=generator.uniform(0, 1e6, n), z1=height, z2=height,
                  Ci=np.where(generator.random(n) < 0.5, np.nan, generator.uniform(100, 400, n)))
    for name, value in values.items():
        ds.set(name, value)
    ds.set("diameter", generator.uniform(0.001, 0.005, len(ds.entity_ids("Organ"))))
    return dict(values, proteins=proteins)


@pytest.mark.parametrize("model", ("BWB", "Tuzet"))
def test_the_component_reproduces_the_legacy(shoots, model, monkeypatch, quiet):
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    monkeypatch.setattr(component, "stomatal_model", model)
    monkeypatch.setattr(component, "corrections", ())
    weather.update(air_temperature=22., relative_humidity=0.4, ambient_CO2=410., wind_speed=1.5)
    inputs = _set_inputs(ds, seed=5)
    previous_Ci = inputs["Ci"].copy()
    component()

    labels, organs = ds.labels(), ds.labels(location="Organ")
    plant, organ = ds.owner("Plant"), ds.owner("Organ")
    exchanging = np.isin(labels, ["StemElement", "LeafElement"])
    computed = exchanging & (inputs["surfacic_nitrogen"] != 0)
    assert computed.sum() > 20 and (organs[computed] == "Blade").any() and (organs[computed] != "Blade").any()
    canopy = np.asarray(ds.get("canopy_height"))
    for p in range(canopy.size):
        assert canopy[p] == inputs["z1"][exchanging & (plant == p)].max()
    outputs = {name: np.asarray(ds.get(name)) for name in OUTPUTS}
    diameter = np.asarray(ds.get("diameter"))
    for p, changes in enumerate(({}, {"THETA": 0.8, "EFFICENCY_STEM": 0.6})):     # the second plant's parameters
        with monkeypatch.context() as patched:
            for name, value in changes.items():
                patched.setattr(legacy_parameters, name, value)
            nodes = np.flatnonzero(computed & (plant == p))
            elements = {}
            for i in nodes.tolist():
                lamina = organs[i] == "Blade"
                elements[(1, "MS", i, "blade" if lamina else "sheath", labels[i])] = python_floats(dict(
                    width=inputs["width"][i] if lamina else diameter[organ[i]], height=inputs["z1"][i],
                    PARa=inputs["PARa"][i], nitrates=0., amino_acids=0., Nstruct=0., green_area=1.,
                    proteins=inputs["proteins"][i], sucrose=inputs["surfacic_NSC"][i], starch=0., fructan=0.,
                    PARa_prim=[], area_prim=0, Ci=previous_Ci[i], water_potential=-0.1))
            legacy = legacy_run(elements, (22., 410., 0.4, 1.5), model, canopy_height=float(canopy[p]))
            for name in OUTPUTS:
                assert_arrays_equal([legacy[key][name] for key in elements], outputs[name][nodes], name=name)
    assert_arrays_equal(outputs["Tr"] * inputs["green_area"], ds.get("transpiration"), name="transpiration")
    others = ~computed
    assert (outputs["Ag"][others] == 0).all() and (outputs["Tr"][others] == 0).all()
    hidden = labels == "HiddenElement"
    assert (outputs["Ts"][hidden] == 22.).all()


def _halved(green_area):
    return 0.5 * np.asarray(green_area)


class HalfExchangingShoot(ExchangingShoot):
    """The gas exchange over half the elements' green area: a link of a receiver-private name (metafspm's
    Link(names=)), as an assembly gives it the green area above the soil (emerged_links)."""
    translator = Translator(port.adel_links().links).link("GasExchange", "green_area", "ShootArchitecture",
                                                          ("green_area",), formula=_halved,
                                                          names={"green_area": "half_green_area"})


def test_an_element_exchanges_over_the_green_area_it_reads(shoots, quiet):
    """The area an element exchanges over (exchanging_area) is the green area the gas exchange reads: given by a link
    of a receiver-private name (here half the architecture's), its fluxes per unit area as over the whole, its
    transpiration Tr times that area, the per-area variables weighted by it; the architecture's green_area keeps its
    meaning (the element's, untouched)."""
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    weather.update(air_temperature=20., relative_humidity=0.6, ambient_CO2=400., wind_speed=2.)
    inputs = _set_inputs(ds, seed=7)
    component()
    whole = {name: np.array(ds.get(name)) for name in OUTPUTS + ("exchanging_area", "transpiration")}
    computed = whole["exchanging_area"] > 0
    np.testing.assert_array_equal(whole["exchanging_area"][computed], inputs["green_area"][computed])
    assert next(f for f in fields(port.GasExchange) if f.name == "An").metadata["weight"] == "exchanging_area"
    Choregrapher().reset()
    halved = Scene(planting([Planting(), Planting(x=0.2)], HalfExchangingShoot, {},
                            [{}, {"THETA": 0.8, "stem_efficiency": 0.6}]),
                   time_step=86400., forcings={name: (lambda time, name=name: weather[name]) for name in weather})
    for _ in range(60):
        halved.run()
    half_population = halved.populations[0]
    half_ds, half_component = half_population.data_structure, half_population.instance.gas_exchange
    assert half_ds.n_nodes() == ds.n_nodes()
    _set_inputs(half_ds, seed=7)
    half_ds.set("Ci", inputs["Ci"])
    half_component()
    cut = {name: np.array(half_ds.get(name)) for name in OUTPUTS + ("exchanging_area", "transpiration")}
    np.testing.assert_array_equal(half_ds.get("green_area"), inputs["green_area"])      # the architecture's
    np.testing.assert_array_equal(half_component.read("green_area"), 0.5 * inputs["green_area"])
    np.testing.assert_array_equal(cut["exchanging_area"][computed], 0.5 * inputs["green_area"][computed])
    for name in OUTPUTS:                                             # per unit area: as over the whole green area
        np.testing.assert_allclose(cut[name][computed], whole[name][computed], rtol=1e-12, err_msg=name)
    np.testing.assert_allclose(cut["transpiration"], cut["Tr"] * cut["exchanging_area"], rtol=1e-15)


def test_the_component_declarations_and_the_organs(shoots, quiet):
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    weather.update(air_temperature=18.8, relative_humidity=0.53, ambient_CO2=360., wind_speed=2.2)
    _set_inputs(ds, seed=6)
    component()
    table = port.GasExchange.declarations().set_index("name")
    assert table.loc["transpiration", "kind"] == "extensive" and table.loc["Tr", "kind"] == "intensive"
    assert set(port.GasExchange.declared_forcings()) == {"air_temperature", "relative_humidity", "ambient_CO2",
                                                         "wind_speed"}
    assert table.loc["THETA", "min_value"] == 0. and table.loc["THETA", "references"]
    organ_transpiration = ds.map_to(ds.get("transpiration"), "node", "Organ", "sum")
    assert organ_transpiration.sum() == pytest.approx(np.sum(ds.get("transpiration")))
    assert (organ_transpiration > 0).sum() > 5


def test_the_component_restarts_Ci(shoots, quiet):
    """An element's Ci of 0 (the legacy's output for an element without exchange): the next step starts from the
    ambient CO2 (the legacy divides by it)."""
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    _set_inputs(ds, seed=8)
    ds.set("Ci", np.zeros(ds.n_nodes()))
    component()
    computed = np.isin(ds.labels(), ["StemElement", "LeafElement"]) & (np.asarray(ds.get("surfacic_nitrogen")) > 0)
    assert np.isfinite(ds.get("Ci")).all() and (np.asarray(ds.get("Ci"))[computed] > 0).all()


@pytest.mark.filterwarnings("ignore:GasExchange:RuntimeWarning")   # some strongly closed elements oscillate
def test_the_tuzet_feedback(shoots, monkeypatch, quiet):
    """Tuzet: a lower leaf water potential closes the stomata (the potential an input, e.g. MECHA's)."""
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    monkeypatch.setattr(component, "stomatal_model", "Tuzet")
    _set_inputs(ds, seed=7)
    ds.set("PARa", np.full(ds.n_nodes(), 1000.))
    if not ds.has("leaf_water_potential"):
        ds.register("leaf_water_potential", -0.1, location="node")
    conductances = []
    for potential in (-0.1, -0.65, -1.5):
        ds.set("leaf_water_potential", np.full(ds.n_nodes(), potential))
        component()
        conductances.append(np.asarray(ds.get("gs")).copy())
    lit = conductances[0] > 0.1
    assert lit.sum() > 5
    assert (conductances[0][lit] > conductances[1][lit]).all() and (conductances[1][lit] > conductances[2][lit]).all()


def test_the_component_at_the_primitive_scale(shoots, monkeypatch, quiet):
    """primitive_scale: each element's primitives (a light model's PARa_primitives and primitive_areas, object
    variables), as the legacy's prim_scale with corrections=()."""
    scene, weather = shoots
    population = scene.populations[0]
    ds, component = population.data_structure, population.instance.gas_exchange
    monkeypatch.setattr(component, "primitive_scale", True)
    monkeypatch.setattr(component, "corrections", ())
    monkeypatch.setattr(legacy_parameters, "prim_scale", True)
    weather.update(air_temperature=18.8, relative_humidity=0.53, ambient_CO2=360., wind_speed=2.2)
    inputs = _set_inputs(ds, seed=9)
    generator = np.random.default_rng(9)
    n = ds.n_nodes()
    PARs = [list(generator.uniform(0, 2000, int(generator.integers(1, 4)))) for _ in range(n)]
    areas = [list(generator.uniform(1e-5, 1e-4, len(values))) for values in PARs]
    for name, values in (("PARa_primitives", PARs), ("primitive_areas", areas)):
        if not ds.has(name):
            ds.register(name, None, location="node", dtype="object")
        ds.set(name, np.array(values + [None], dtype=object)[:-1])
    previous_Ci = inputs["Ci"].copy()
    component()
    labels, organs, plant = ds.labels(), ds.labels(location="Organ"), ds.owner("Plant")
    nodes = np.flatnonzero(np.isin(labels, ["StemElement", "LeafElement"]) & (inputs["surfacic_nitrogen"] != 0)
                           & (plant == 0))
    diameter = np.asarray(ds.get("diameter"))[ds.owner("Organ")]
    elements = {(1, "MS", i, "blade" if organs[i] == "Blade" else "sheath", labels[i]): python_floats(dict(
        width=inputs["width"][i] if organs[i] == "Blade" else diameter[i], height=inputs["z1"][i], PARa=0.,
        nitrates=0., amino_acids=0., Nstruct=0., green_area=1., proteins=inputs["proteins"][i],
        sucrose=inputs["surfacic_NSC"][i], starch=0., fructan=0., PARa_prim=PARs[i], area_prim=areas[i],
        Ci=previous_Ci[i])) for i in nodes.tolist()}
    legacy = legacy_run(elements, (18.8, 360., 0.53, 2.2), canopy_height=float(ds.get("canopy_height")[0]))
    for name in OUTPUTS:
        assert_arrays_equal([legacy[key][name] for key in elements], np.asarray(ds.get(name))[nodes], name=name)


class SownExchangingShoot(ExchangingShoot):
    """The gas exchange of a sown shoot over its elements' parts above the soil (emerged_links)."""
    translator = port.emerged_links()


def test_a_sown_shoot_exchanges_over_its_parts_above_the_soil(quiet):
    """emerged_links: each element's green area above the soil surface (the Scene's), from its green triangles'
    share there, read by the gas exchange as its green_area (Adel's untouched); its height the centre of its part
    above the soil, NaN wholly below (no gas exchange)."""
    Choregrapher().reset()
    weather = dict(air_temperature=18.8, relative_humidity=0.53, ambient_CO2=360., wind_speed=2.2)
    scene = Scene(planting([Planting(z=-0.05)], SownExchangingShoot, {}), time_step=86400., soil_surface=0.,
                  forcings={name: (lambda time, name=name: weather[name]) for name in weather})
    for _ in range(40):
        scene.run()
    ds, component = scene.populations[0].data_structure, scene.populations[0].instance.gas_exchange
    green_area, read = np.asarray(ds.get("green_area")), component.read("green_area")
    assert (read <= green_area * (1 + 1e-12)).all() and (read < green_area).any() and (read > 0).any()
    triangles = ds.triangles()
    lowest = np.full(ds.n_nodes(), np.inf)
    np.minimum.at(lowest, triangles.owner, triangles.vertices[..., 2].min(axis=1))
    highest = np.full(ds.n_nodes(), -np.inf)
    np.maximum.at(highest, triangles.owner, triangles.vertices[..., 2].max(axis=1))
    drawn = np.bincount(triangles.owner, minlength=ds.n_nodes()) > 0
    buried, emerged = drawn & (highest < 0), drawn & (lowest >= 0)
    assert buried.any() and (read[buried] == 0).all() and np.isnan(component.read("height")[buried]).all()
    np.testing.assert_allclose(read[emerged], green_area[emerged], rtol=1e-12)
    assert (component.read("height")[emerged] > 0).all()
