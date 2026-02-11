## FSPM Wheat
## Adel master du 2018-11-30


from openalea.adel.Stand import AgronomicStand
from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.AdelR import devCsv
from openalea.adel.plantgen_extensions import TillerEmission, TillerRegression, \
    AxePop, PlantGen, HaunStage
from openalea.adel.echap_leaf import echap_leaves

import os
import pandas as pd

## -- Create a full MTG

def create_init_MTG_with_tillers(nplants=1, sowing_density=250., plant_density=250.,
                           inter_row=0.15, nff=12, nsect=1, seed=1, leaves=echap_leaves(xy_model='Soissons_byleafclass')):
    stand = AgronomicStand(sowing_density=sowing_density,
                           plant_density=plant_density, inter_row=inter_row,
                           noise=0.04, density_curve_data=None)
    em = TillerEmission(
        primary_tiller_probabilities={'T1': 1., 'T2': 1, 'T3': 1,
                                      'T4': 1})
    reg = TillerRegression(ears_per_plant=3)
    axp = AxePop(MS_leaves_number_probabilities={str(nff): 1}, Emission=em, Regression=reg)
    plants = axp.plant_list(nplants=nplants)
    hs = HaunStage(mean_nff=nff)
    pgen = PlantGen(HSfit=hs)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort_values(['id_plt', 'id_cohort', 'N_phytomer'])
    devT = devCsv(axeT, dimT, phenT)
    adel = AdelDyn(nplants=nplants, nsect=nsect, devT=devT, stand=stand,
                     seed=seed, sample='sequence', leaves=leaves,  scene_unit='m')
    age = hs.TT(reg.hs_debreg(nff=nff)) # date a laquelle debut des regression. Problemes : 1) toutes les feuilles pas encore visibles, 2) il y a des feuilles senescentes
    g = adel.setup_canopy(age)
    return adel, g

adel, g = create_init_MTG_with_tillers(nff = 14)

adel.plot(g)

## -- Adapt the dimensions to inputs data

INPUTS_DIRPATH = ''
HIDDENZONE_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'hiddenzones_initial_state.csv')
ELEMENTS_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'elements_initial_state.csv')
hdz_inputs = pd.read_csv(HIDDENZONE_INPUTS_FILEPATH)
elt_inputs = pd.read_csv(ELEMENTS_INPUTS_FILEPATH)


elt_inputs.set_index(['plant', 'axis', 'metamer', 'organ', 'element'], inplace=True)
elt_dict = elt_inputs.to_dict('index')

org_id_list = [ i[:4] for i in elt_dict.keys() ]
width_dict = { (1,'MS',1,'blade'): 0.0030,(1,'MS',2,'blade'): 0.0033, (1,'MS',1,'sheath'): 0.003, (1,'MS',2,'sheath'): 0.003 }
age_dict = {(1,'MS',1,'blade'): 280,(1,'MS',2,'blade'): 180}

for vid_axe in g.components_at_scale(1,2):
    for vid_meta in g.components_at_scale(vid_axe,3):
        for vid_org in g.components_at_scale(vid_meta,4):
            org_id = (1, g.property('label')[vid_axe],  int(filter(str.isdigit,g.property('label')[vid_meta] )), g.property('label')[vid_org] )
            g.property('senesced_length')[vid_org] = 0.
            if org_id not in org_id_list:
                g.property('visible_length')[vid_org] = 0.
            elif org_id in width_dict.keys() :
                length_max = 0
                for k, v in elt_dict.iteritems():
                    if k[:4] == org_id:
                        length_max+= v['length']
                g.property('shape_mature_length')[vid_org] = length_max
                g.property('length')[vid_org] = length_max
                g.property('senesced_length')[vid_org] = 0.

                for vid_elt in g.components_at_scale(vid_org,5):
                    g.property('senesced_length')[vid_elt] = 0.
                    g.property('senesced_area')[vid_elt] = 0.

                if g.property('label')[vid_org] == 'blade' :
                    g.property('width')[vid_org] = width_dict[org_id]
                    g.property('shape_max_width')[vid_org] = width_dict[org_id]
                    g.property('visible_length')[vid_org] =  elt_dict[ tuple(list(org_id) + ['LeafElement1']) ]['length']
                    g.property('age')[vid_org] = age_dict[org_id]
                elif g.property('label')[vid_org] == 'sheath' :
                    g.property('diameter')[vid_org] = width_dict[org_id]
                    g.property('visible_length')[vid_org] = elt_dict[tuple(list(org_id) + ['StemElement'])]['length']

adel.update_geometry(g)

adel.plot(g)

# Save mtg
adel.save(g)
# save adel pars
adel.save_pars()

adel2 = AdelDyn(seed=1, scene_unit='m',leaves=echap_leaves(xy_model='Soissons_byleafclass'))
adel2.pars = adel2.read_pars()
g2 = adel2.load()
adel2.update_geometry(g2)
adel2.plot(g2)


