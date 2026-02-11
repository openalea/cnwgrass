from openalea.adel.dresser import blade_dimension, stem_dimension, dimension_table, AdelDressDyn
from openalea.adel.Stand import AgronomicStand
from openalea.adel.adel_dynamic import AdelWheatDyn
from openalea.adel.AdelR import devCsv
from openalea.adel.plantgen_extensions import TillerEmission, TillerRegression, \
    AxePop, PlantGen, HaunStage


def create_init_MTG_with_tillers(nplants=1, sowing_density=250., plant_density=250.,
                           inter_row=0.15, nff=12, nsect=1, seed=1, leaves=None):
    stand = AgronomicStand(sowing_density=sowing_density,
                           plant_density=plant_density, inter_row=inter_row,
                           noise=0.04, density_curve_data=None)
    em = TillerEmission(
        primary_tiller_probabilities={'T1': 1., 'T2': 0.5, 'T3': 0.5,
                                      'T4': 0.3})
    reg = TillerRegression(ears_per_plant=3)
    axp = AxePop(MS_leaves_number_probabilities={str(nff): 1}, Emission=em, Regression=reg)
    plants = axp.plant_list(nplants=nplants)
    hs = HaunStage(mean_nff=nff)
    pgen = PlantGen(HSfit=hs)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort_values(['id_plt', 'id_cohort', 'N_phytomer'])
    devT = devCsv(axeT, dimT, phenT)
    adel = AdelWheatDyn(nplants=nplants, nsect=nsect, devT=devT, stand=stand,
                     seed=seed, sample='sequence', leaves=leaves,  scene_unit='m')
    age = hs.TT(reg.hs_debreg(nff=nff)) # date a laquelle debut des regression. Problemes : 1) toutes les feuilles pas encore visibles, 2) il y a des feuilles senescentes
    g = adel.setup_canopy(age) # MG : je ne crois pas que id_cohort soit renvoyee par la fonction R runAdel
    return adel, g

def create_init_MTG_FSPMWHEAT(dirpath):
    """
    Creates the initial MTG used for FSPM-Wheat model

    : Parameters:
        - `dirpath` (:class:`str`) - The path to save the MTG

    : Returns:
        adel, g

    :Returns Type:
        :class:`Adel`, :class:`MTG`
    """

    # Organ dimensions
    # TODO : complete with zero up to desired  phtomer number
    blades = blade_dimension(length=[ 0.082, 0.092],
                             width = [0.0030,0.0033], #area=[ 0.00024, 0.00028],
                             ntop=[ 2, 1])

    stem = stem_dimension(ntop=[2, 1],
                          sheath=[ 0.0285, 0.029], d_sheath=[ 0.003, 0.003],
                          internode=[0, 0], d_internode=[ 0, 0])

    dimT = dimension_table(blades, stem)

    # Creates the stand
    stand = AgronomicStand(sowing_density=250, plant_density=250, inter_row=0.17, noise=0.03)

    # Adel instantiation
    adel = AdelDressDyn(dimT=dimT, dim_unit='m', nplants= 1, scene_unit='m', seed=1234)
    g = adel.canopy()
    adel.plot(g)

    # Save MTG
    adel.save(g, dir=dirpath)
    return adel, g

adel, g = create_init_MTG_FSPMWHEAT(r'./adelwheat')
# adel, g = create_init_MTG_with_tillers(nff = 14)
# adel.save(g, dir= r'./adelwheat' )
