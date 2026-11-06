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
    def initiate_plant(cls, g, plant, parameters):
        s, L = g.scales, g.labels

        def organ(axis_label, kind, organ_label):
            axis = g.add_component(plant, **PropsConfig(scale=s.Axis, edge_type='/', label=getattr(L.Axis, axis_label)))
            unit = g.add_component(axis, **PropsConfig(scale=s.GrowthUnit, edge_type='/',
                                                       label=getattr(L.GrowthUnit, kind)))
            phytomer = g.add_component(unit, **PropsConfig(scale=s.Phytomer, edge_type='/',
                                                           label=getattr(L.Phytomer, kind)))
            return g.add_component(phytomer, **PropsConfig(scale=s.Organ, edge_type='/',
                                                           label=getattr(L.Organ, organ_label)))

        stem = g.add_component(organ("Shoot", "Shoot", "StemInternode"),
                               **PropsConfig(scale=s.SubOrgan, edge_type='/', label=L.SubOrgan.StemElement))
        leaf = g.add_component_with_topo(organ("Shoot", "Shoot", "Leaf"), stem,
                                         **PropsConfig(scale=s.SubOrgan, edge_type='+', label=L.SubOrgan.LeafElement))
        root = g.add_component_with_topo(organ("Root", "Root", "RootInternode"), stem,
                                         **PropsConfig(scale=s.SubOrgan, edge_type='+', label=L.SubOrgan.RootSegment))
        g.property("stem_radius")[stem] = 1e-3
        for vertex, radius in ((stem, 1e-3), (leaf, 1e-3), (root, 0.55e-3)):
            g.property("radius")[vertex] = radius
            g.property("length")[vertex] = 5e-3
            g.property("distance_from_tip")[vertex] = 0.


class Anatomy(GranapAnatomy):
    shoot = {"StemElement": OrganAnatomy(AnatomyLibrary(parameters=preset_parameters("for_monocot_stem"),
                                                        organ_class=StemAnatomy)),
             "LeafElement": OrganAnatomy(AnatomyLibrary(parameters=preset_parameters("for_monocot_leaf"),
                                                        organ_class=LeafAnatomy), size=None)}


class Initiated:
    """The plants' anatomies, wired at the segments' junctions; no component (the test couples them)."""
    initiators, nodes, wiring = (Seedling, Anatomy), "Compartment", staticmethod(junction_rules)

    def __init__(self, data_structure, time_step):
        self.components = []

    def run(self):
        pass


class Coupled(CompositeModel):
    """MECHA's hydraulics and the gas exchange on one DataStructure, the leaf water potential linked."""
    translator = mecha_links()

    def __init__(self, data_structure, time_step):
        self.hydraulics = MechaHydraulics(data_structure=data_structure, time_step=time_step)
        self.gas_exchange = GasExchange(data_structure=data_structure, time_step=time_step)
        # an optional input, not registered until given: held at the elements, for the link to map onto them
        data_structure.register("leaf_water_potential", GasExchange.leaf_water_potential, location="SubOrgan")
        self.couple(self.hydraulics, self.gas_exchange)


@pytest.fixture(scope="module")
def coupled():
    ds = population_data_structure(planting([Planting()], None), Initiated, time_step=3600.).data_structure
    Anatomy(data_structure=ds)()
    ds.register("barrier", np.ones(ds.n_nodes()), location="node")
    ds.register("soil_psi", np.full(ds.n_nodes(), SOIL), location="node")
    leaf = np.asarray(ds.labels(location="Organ", per_node=False)) == "Leaf"
    ds.register("leaf_transpiration", np.where(leaf, TRANSPIRATION, 0.), location="Organ")
    model = Coupled(ds, 3600.)
    model.hydraulics.collar_condition = 0.
    model.hydraulics()
    return ds, model


def test_the_gas_exchange_names_its_leaf_water_potential_apart_from_mechas():
    names = {f.name for f in fields(GasExchange)}
    assert "leaf_water_potential" in names and "water_potential" not in names
    assert GasExchange.__dataclass_fields__["leaf_water_potential"].metadata["unit"] == "MPa"
    assert "leaf_water_potential" not in {f.name for f in fields(MechaHydraulics)}
    (link,) = mecha_links().links
    assert (link.receiver, link.variable, link.provider) == ("GasExchange", "leaf_water_potential", "MechaHydraulics")
    assert link.sources == {"water_potential": 1e-4} and link.unit_check            # hPa -> MPa, checked


def test_each_element_gets_the_mean_of_its_liquid_in_mpa(coupled):
    ds, _ = coupled
    elements = ds.entity_ids("SubOrgan")
    owner = np.asarray(ds.owner("SubOrgan"))
    water, liquid = np.asarray(ds.get("water_potential")), np.asarray(ds.get("liquid_volume"))
    vapour = np.asarray(ds.get("vapour_node")) > 0
    expected = np.array([np.sum(water[owner == i] * liquid[owner == i]) / np.sum(liquid[owner == i])
                         for i in range(elements.size)]) * 1e-4
    assert ds.location("leaf_water_potential") == "SubOrgan"
    np.testing.assert_allclose(ds.get("leaf_water_potential"), expected, rtol=1e-12)
    leaf = ds.labels(location="SubOrgan", per_node=False) == "LeafElement"
    index = int(np.flatnonzero(leaf)[0])
    liquid_range = water[(owner == index) & ~vapour] * 1e-4
    assert liquid_range.min() <= ds.get("leaf_water_potential")[index] <= liquid_range.max() < 0
