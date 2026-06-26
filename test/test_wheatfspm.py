import os
import pandas as pd

from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves
from openalea.integration import caribu_facade
from openalea.integration import cnmetabolism_facade
from openalea.integration import morphogenesis_facade
from openalea.integration import gasexchange_facade
from openalea.integration import growth_facade
from openalea.integration import senescence_facade
from openalea.integration import hydraulics_facade

# -- SIMULATION PARAMETERS --
# Length of the simulation (in hours)
SIMULATION_LENGTH = 5

# define the time step in hours for each simulator
CARIBU_TIMESTEP = 4
SENESCENCE_TIMESTEP = 1
GASEXCHANGE_TIMESTEP = 1
MORPHOGENESIS_TIMESTEP = 1
GROWTH_TIMESTEP = 1
CNMETABOLISM_TIMESTEP = 1
hydraulics_TIMESTEP = 1

# Define default plant density (culm m-2)
PLANT_DENSITY = {1: 250.}

# precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 2

# number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# Name of the CSV files which will contain the outputs of the model
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

# -- INPUTS CONFIGURATION --

# Path of the directory which contains the inputs of the model
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUTS_DIRPATH = os.path.join(SCRIPT_DIR, 'inputs')

# Name of the CSV files which describes the initial state of the system
AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'

# Read the inputs from CSV files and create inputs dataframes
inputs_dataframes = {}
for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                        ORGANS_INITIAL_STATE_FILENAME,
                        HIDDENZONES_INITIAL_STATE_FILENAME,
                        ELEMENTS_INITIAL_STATE_FILENAME,
                        SOILS_INITIAL_STATE_FILENAME):
    inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
    inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)

# -- OUTPUTS CONFIGURATION --

# Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
UPDATE_SHARED_DF = False

# create empty dataframes to shared data between the models
shared_axes_inputs_outputs_df = pd.DataFrame()
shared_organs_inputs_outputs_df = pd.DataFrame()
shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
shared_elements_inputs_outputs_df = pd.DataFrame()
shared_soils_inputs_outputs_df = pd.DataFrame()

# define lists of dataframes to store the inputs and the outputs of the models at each step.
axes_all_data_list = []
organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
hiddenzones_all_data_list = []
elements_all_data_list = []
soils_all_data_list = []

all_simulation_steps = []  # to store the steps of the simulation

# -- POSTPROCESSING CONFIGURATION --

# Name of the CSV files which will contain the postprocessing of the model
AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'


# -- ADEL and MTG CONFIGURATION --
adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
g = adel_wheat.load(directory=INPUTS_DIRPATH)


def test_cnmetabolism():
    # -- CNMETABOLISM --
    # Initial states
    cnmetabolism_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        [i for i in cnmetabolism_facade.cnmetabolism_converter.AXES_VARIABLES if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    cnmetabolism_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in cnmetabolism_facade.cnmetabolism_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    cnmetabolism_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in cnmetabolism_facade.cnmetabolism_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    cnmetabolism_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in cnmetabolism_facade.cnmetabolism_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    cnmetabolism_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in cnmetabolism_facade.cnmetabolism_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    cnmetabolism_facade_ = cnmetabolism_facade.CNMetabolismFacade(g,
                                                   CNMETABOLISM_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   {},
                                                   cnmetabolism_axes_initial_state,
                                                   cnmetabolism_organs_initial_state,
                                                   cnmetabolism_hiddenzones_initial_state,
                                                   cnmetabolism_elements_initial_state,
                                                   cnmetabolism_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df,
                                                   tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5},
                                                   update_shared_df=UPDATE_SHARED_DF)

    # Run facade
    cnmetabolism_facade_.run(12, 10)


def test_morphogenesis():
    # Initial states
    morphogenesis_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
    morphogenesis_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    morphogenesis_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

    phytoT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'phytoT.csv'))

    # Facade initialisation
    morphogenesis_facade_ = morphogenesis_facade.MorphogenesisFacade(g,
                                                            MORPHOGENESIS_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            morphogenesis_axes_initial_state,
                                                            morphogenesis_hiddenzones_initial_state,
                                                            morphogenesis_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT_df,
                                                            update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    morphogenesis_facade_.run(12, 10)
    adel_wheat.update_geometry(g)


def test_caribu():
    caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                shared_elements_inputs_outputs_df,
                                                adel_wheat,
                                                update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    caribu_facade_.run(True, energy=1000, DOY=1, hourTU=12, latitude=48.85, sun_sky_option='sky',
                       heterogeneous_canopy=False, plant_density=PLANT_DENSITY[1])

def test_gasexchange():
    # adhoc modification for test
    g.property('diameter').update({19: 0.003, 34: 0.003})

    # Initial states
    gasexchange_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    gasexchange_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

    # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
    update_parameters_gasexchange = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

    # Facade initialisation
    gasexchange_facade_ = gasexchange_facade.GasExchangeFacade(g,
                                                                     gasexchange_elements_initial_state,
                                                                     gasexchange_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df,
                                                                     update_parameters=update_parameters_gasexchange,
                                                                     update_shared_df=UPDATE_SHARED_DF)

    # Run facade
    gasexchange_facade_.run(12, 400, 0.8, 1.5)

def test_growth():
    # Initial states
    growth_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
    growth_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    growth_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]
    growth_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]

    # Facade initialisation
    growth_facade_ = growth_facade.GrowthFacade(g,
                                                               GROWTH_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               growth_hiddenzones_initial_state,
                                                               growth_elements_initial_state,
                                                               growth_root_initial_state,
                                                               growth_axes_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_hiddenzones_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               update_shared_df=UPDATE_SHARED_DF)
        # Run facade
    growth_facade_.run()

def test_senescence():
    # Initial states
    senescence_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[
        inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        senescence_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
        [i for i in senescence_facade.converter.SENESCENCE_ROOTS_INPUTS if
         i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    senescence_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        senescence_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
        [i for i in senescence_facade.converter.SENESCENCE_ELEMENTS_INPUTS if
         i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    senescence_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        senescence_facade.converter.AXES_TOPOLOGY_COLUMNS +
        [i for i in senescence_facade.converter.SENESCENCE_AXES_INPUTS if
         i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    senescence_facade_ = senescence_facade.SENESCENCEFacade(g,
                                                               SENESCENCE_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescence_roots_initial_state,
                                                               senescence_axes_initial_state,
                                                               senescence_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    senescence_facade_.run()

def test_hydraulics():
    # Initial states
    hydraulics_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        [i for i in hydraulics_facade.hydraulics_converter.AXES_VARIABLES if
         i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    hydraulics_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in hydraulics_facade.hydraulics_converter.ORGANS_VARIABLES if
         i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    hydraulics_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in hydraulics_facade.hydraulics_converter.HIDDENZONE_VARIABLES if
         i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    hydraulics_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in hydraulics_facade.hydraulics_converter.ELEMENTS_VARIABLES if
         i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    hydraulics_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in hydraulics_facade.hydraulics_converter.SOILS_VARIABLES if
         i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    hydraulics_facade_ = hydraulics_facade.hydraulicsFacade(g,
                                                                  hydraulics_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                  {},
                                                                  hydraulics_axes_initial_state,
                                                                  hydraulics_hiddenzones_initial_state,
                                                                  hydraulics_elements_initial_state,
                                                                  hydraulics_organs_initial_state,
                                                                  hydraulics_soils_initial_state,
                                                                  shared_axes_inputs_outputs_df,
                                                                  shared_hiddenzones_inputs_outputs_df,
                                                                  shared_elements_inputs_outputs_df,
                                                                  shared_organs_inputs_outputs_df,
                                                                  shared_soils_inputs_outputs_df,
                                                                  update_shared_df=UPDATE_SHARED_DF)

    # Run facade
    hydraulics_facade_.run()
