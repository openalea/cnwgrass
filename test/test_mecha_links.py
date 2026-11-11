"""
GasExchange's leaf water potential from MECHA (mecha_links): on an anatomy-mode DataStructure (GRANAP's sections,
MECHA's network over their Compartments), each element's leaf_water_potential (MPa) is the mean of MECHA's
water_potential (hPa) over its liquid nodes (weighted by liquid_volume, the air spaces and pores excluded), a name
of its own beside MECHA's: the two potentials' units differ.
"""
from dataclasses import fields

import numpy as np
import pytest

from openalea.cnwgrass.metafspm.gas_exchange import GasExchange, mecha_links

granap = pytest.importorskip("openalea.granap.metafspm.component")
mecha = pytest.importorskip("openalea.mecha.component")

from openalea.metafspm.coupling.component import StructuralComponent  # noqa: E402
from openalea.metafspm.coupling.composite_wrapper import CompositeModel  # noqa: E402
from openalea.metafspm.data_structure.configs import PropsConfig  # noqa: E402
from openalea.metafspm.scene.population import Planting, planting, population_data_structure  # noqa: E402

from openalea.granap.leaf_class import LeafAnatomy  # noqa: E402
from openalea.granap.metafspm.component import GranapAnatomy, OrganAnatomy, junction_rules  # noqa: E402
from openalea.granap.metafspm.library import AnatomyLibrary, preset_parameters  # noqa: E402
from openalea.granap.stem_class import StemAnatomy  # noqa: E402
from openalea.mecha.component import MechaHydraulics  # noqa: E402

SOIL = -1000.                                            # hPa
TRANSPIRATION = 0.05                                     # cm3 d-1, the leaf organ's


class Seedling(StructuralComponent):
    """A stem element bearing a leaf element and a root segment, 5 mm long each."""

    @classmethod
    def initiate_plant(cls, mtg, plant, parameters):
        s, L = mtg.scales, mtg.labels

        def organ(axis_label, kind, organ_label):
            axis = mtg.add_component(plant, **PropsConfig(scale=s.Axis, edge_type='/', label=getattr(L.Axis, axis_label)))
            unit = mtg.add_component(axis, **PropsConfig(scale=s.GrowthUnit, edge_type='/',
                                                       label=getattr(L.GrowthUnit, kind)))
            phytomer = mtg.add_component(unit, **PropsConfig(scale=s.Phytomer, edge_type='/',
                                                           label=getattr(L.Phytomer, kind)))
            return mtg.add_component(phytomer, **PropsConfig(scale=s.Organ, edge_type='/',
                                                           label=getattr(L.Organ, organ_label)))

        stem = mtg.add_component(organ("Shoot", "Shoot", "StemInternode"),
                               **PropsConfig(scale=s.SubOrgan, edge_type='/', label=L.SubOrgan.StemElement))
        leaf = mtg.add_component_with_topo(organ("Shoot", "Shoot", "Leaf"), stem,
                                         **PropsConfig(scale=s.SubOrgan, edge_type='+', label=L.SubOrgan.LeafElement))
        root = mtg.add_component_with_topo(organ("Root", "Root", "RootInternode"), stem,
                                         **PropsConfig(scale=s.SubOrgan, edge_type='+', label=L.SubOrgan.RootSegment))
        mtg.set_values("stem_radius", [stem], [1e-3])
        vertices = [stem, leaf, root]
        mtg.set_values("radius", vertices, [1e-3, 1e-3, 0.55e-3])
        mtg.set_values("length", vertices, [5e-3] * 3)
        mtg.set_values("distance_from_tip", vertices, [0.] * 3)


class Anatomy(GranapAnatomy):
    shoot = {"StemElement": OrganAnatomy(AnatomyLibrary(parameters=preset_parameters("for_monocot_stem"),
                                                        organ_class=StemAnatomy)),
             "LeafElement": OrganAnatomy(AnatomyLibrary(parameters=preset_parameters("for_monocot_leaf"),
                                                        organ_class=LeafAnatomy), size=None)}


class Initiated(CompositeModel):
    """The plants' anatomies, wired at the segments' junctions; no component (the test couples them)."""
    initiators, nodes, wiring = (Seedling, Anatomy), "Compartment", staticmethod(junction_rules)


class Coupled(CompositeModel):
    """MECHA's hydraulics and the gas exchange on one DataStructure, the leaf water potential linked."""
    translator = mecha_links()
    parts = {"hydraulics": MechaHydraulics, "gas_exchange": GasExchange}

    def __init__(self, data_structure, time_step):
        # an optional input, not registered until given: held at the elements, for the link to map onto them
        data_structure.register("leaf_water_potential", GasExchange.leaf_water_potential, location="SubOrgan")
        super().__init__(data_structure, time_step)


@pytest.fixture(scope="module")
def coupled():
    data_structure = population_data_structure(planting([Planting()], None), Initiated, time_step=3600.).data_structure
    Anatomy(data_structure=data_structure)()
    data_structure.register("barrier", np.ones(data_structure.n_nodes()), location="node")
    data_structure.register("soil_psi", np.full(data_structure.n_nodes(), SOIL), location="node")
    leaf = np.asarray(data_structure.labels(location="Organ", per_node=False)) == "Leaf"
    data_structure.register("leaf_transpiration", np.where(leaf, TRANSPIRATION, 0.), location="Organ")
    model = Coupled(data_structure, 3600.)
    model.hydraulics.collar_condition = 0.
    model.hydraulics()
    return data_structure, model


def test_the_gas_exchange_names_its_leaf_water_potential_apart_from_mechas():
    names = {f.name for f in fields(GasExchange)}
    assert "leaf_water_potential" in names and "water_potential" not in names
    assert GasExchange.__dataclass_fields__["leaf_water_potential"].metadata["unit"] == "MPa"
    assert "leaf_water_potential" not in {f.name for f in fields(MechaHydraulics)}
    (link,) = mecha_links().links
    assert (link.receiver, link.variable, link.provider) == ("GasExchange", "leaf_water_potential", "MechaHydraulics")
    assert link.sources == {"water_potential": 1e-4} and link.unit_check            # hPa -> MPa, checked


def test_each_element_gets_the_mean_of_its_liquid_in_mpa(coupled):
    data_structure, _ = coupled
    elements = data_structure.entity_ids("SubOrgan")
    owner = np.asarray(data_structure.owner("SubOrgan"))
    water, liquid = np.asarray(data_structure.get("water_potential")), np.asarray(data_structure.get("liquid_volume"))
    vapour = np.asarray(data_structure.get("vapour_node")) > 0
    expected = np.array([np.sum(water[owner == i] * liquid[owner == i]) / np.sum(liquid[owner == i])
                         for i in range(elements.size)]) * 1e-4
    assert data_structure.location("leaf_water_potential") == "SubOrgan"
    np.testing.assert_allclose(data_structure.get("leaf_water_potential"), expected, rtol=1e-12)
    leaf = data_structure.labels(location="SubOrgan", per_node=False) == "LeafElement"
    index = int(np.flatnonzero(leaf)[0])
    liquid_range = water[(owner == index) & ~vapour] * 1e-4
    assert liquid_range.min() <= data_structure.get("leaf_water_potential")[index] <= liquid_range.max() < 0
