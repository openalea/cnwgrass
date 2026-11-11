"""The users' entry point (openalea.cnwgrass.metafspm.simulation): a small wheat stand grown by Adel, lit by Caribu,
its gas exchange hour by hour; its tables against the elements' variables, its parameters, settings and weather."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("openalea.adel.development")
pytest.importorskip("openalea.caribu.component")

from openalea.cnwgrass.metafspm.simulation import HOUR, parameters, simulate  # noqa: E402

SMALL = dict(plants=4, rows=2, thermal_time=300.)


@pytest.fixture(scope="module")
def noon():
    """Two hours, 11 and 12 h, of a small stand."""
    return simulate(hours=2, start_hour=11, **SMALL)


def test_the_summary_is_the_elements(noon):
    """The last hour's row: the elements' fluxes summed per m2 of ground, their temperature and conductance weighted
    by their exchanging area (the whole green area: the stand is above the soil); the hours from start_hour; a clear
    day's light absorbed."""
    summary, elements = noon.summary, noon.elements()
    assert summary["hour"].tolist() == [11., 12.]
    last = summary.iloc[-1]
    computing = elements["exchanging_area"] > 0
    assert computing.any()
    np.testing.assert_array_equal(elements["exchanging_area"][computing], elements["green_area"][computing])
    assert last["net_assimilation"] == pytest.approx((elements["An"] * elements["exchanging_area"]).sum()
                                                     / noon.ground_area, rel=1e-12)
    assert last["transpiration"] == pytest.approx(elements["transpiration"].sum() / noon.ground_area, rel=1e-12)
    assert last["surface_temperature"] == pytest.approx(np.average(elements["Ts"], weights=elements["exchanging_area"]))
    assert 0 < last["PAR_absorbed"] < last["PARi"] and last["net_assimilation"] > 0
    assert last["air_temperature"] == 20.
    blades = elements[elements["organ"] == "blade"]                 # a stem's Ag times stem_efficiency, its An not
    assert (blades["Ag"] >= blades["An"]).all()
    np.testing.assert_allclose(elements["transpiration"], elements["Tr"] * elements["exchanging_area"], rtol=1e-12)
    assert set(elements["organ"]) <= {"blade", "sheath", "internode"} and (elements["organ"] == "blade").any()


def test_the_plants_add_up_the_hours(noon):
    plants = noon.plants()
    hours = noon.plant_hours
    assert len(plants) == SMALL["plants"]
    assert plants["water"].sum() == pytest.approx(hours["transpiration"].sum() * HOUR * 18e-3, rel=1e-12)
    assert plants["net_assimilation"].sum() == pytest.approx(hours["net_assimilation"].sum() * HOUR / 1e3, rel=1e-12)
    last = hours[hours["hour"] == 12.]
    assert last["transpiration"].sum() == pytest.approx(noon.summary.iloc[-1]["transpiration"] * noon.ground_area)


def test_the_stand_grown_to_its_stage(noon):
    """The plants' thermal time: thermal_time minus their delays, and the degree-days of the hours before."""
    data_structure = noon.data_structure
    delays = np.asarray(data_structure.get("delay"))
    np.testing.assert_allclose(np.asarray(data_structure.get("thermal_time")), SMALL["thermal_time"] - delays + 20. / 24.,
                               rtol=1e-12)


def test_the_plants_are_sown_at_the_soil_surface(noon):
    """The plants' bases at the soil surface (their depth 0: not metafspm's stand default), so that every element is
    above it and exchanges over its whole green area."""
    data_structure = noon.data_structure
    np.testing.assert_array_equal(np.asarray(data_structure.get("z")), 0.)
    elements = noon.elements()
    assert (elements["height"] >= 0.).all()
    np.testing.assert_array_equal(elements["exchanging_area"], elements["green_area"])


def test_the_night():
    """In the dark: no assimilation, the stomata at their minimal conductance, the water still lost."""
    night = simulate(hours=1, start_hour=1, **SMALL)
    row = night.summary.iloc[0]
    assert row["PARi"] == 0 and row["gross_assimilation"] == 0 and row["net_assimilation"] == 0
    assert row["stomatal_conductance"] == pytest.approx(parameters().set_index("name").loc[
        "minimum_stomatal_conductance", "value"])
    assert row["transpiration"] > 0


def test_the_parameters_settings_and_weather(noon):
    main = parameters(main_only=True)
    assert main["name"].iloc[0] == "stomatal_model" and main["main"].all()
    assert parameters({"THETA": 0.8}).set_index("name").loc["THETA", "value"] == 0.8
    with pytest.raises(KeyError):
        simulate({"leaf_number": 9}, hours=1, **SMALL)
    leuning = simulate({"stomatal_model": "Leuning"}, hours=1, start_hour=12, **SMALL)
    assert leuning.scene.populations[0].instance.gas_exchange.stomatal_model == "Leuning"
    reference = noon.summary.iloc[-1]
    assert leuning.summary.iloc[0]["stomatal_conductance"] != pytest.approx(reference["stomatal_conductance"])
    weather = simulate(hours=1, start_hour=12, air_temperature=pd.Series([10., 30.], index=[0., 24.]),
                       surfacic_nitrogen=1., **SMALL)
    assert weather.summary.iloc[0]["air_temperature"] == pytest.approx(10. + 20. * 11.5 / 24.)   # 11 to 12 h's mean
    assert np.all(np.asarray(weather.data_structure.get("surfacic_nitrogen")) == 1.)
    assert weather.summary.iloc[0]["net_assimilation"] < reference["net_assimilation"]


def test_the_leaf_water_potential_closes_the_stomata():
    """Tuzet: a lower leaf water potential given to every element, lower conductances."""
    conductances = [simulate({"stomatal_model": "Tuzet"}, hours=1, start_hour=12, leaf_water_potential=potential,
                             **SMALL).summary.iloc[0]["stomatal_conductance"] for potential in (-0.1, -0.65)]
    assert conductances[1] < conductances[0]


def test_the_outputs(tmp_path):
    simulate(hours=1, start_hour=12, outputs=str(tmp_path), **SMALL)
    assert {"run.json", "scenario.json", "summary.csv", "plants.csv"} <= {path.name for path in tmp_path.iterdir()}
    assert len(pd.read_csv(tmp_path / "summary.csv")) == 1


def test_the_run_configuration_is_the_scenes():
    """compute= (metafspm's Compute) is the Scene's run configuration."""
    from openalea.metafspm.compute import Compute
    compute = Compute(threads=1, memory="keep")
    assert simulate(hours=1, start_hour=11, compute=compute, **SMALL).scene.compute is compute
