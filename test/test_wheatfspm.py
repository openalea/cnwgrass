import os
import pandas as pd

from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves
from openalea.fspmwheat import caribu_facade
from openalea.fspmwheat import cnwheat_facade
from openalea.fspmwheat import elongwheat_facade
from openalea.fspmwheat import farquharwheat_facade
from openalea.fspmwheat import fspmwheat_facade
from openalea.fspmwheat import growthwheat_facade
from openalea.fspmwheat import senescwheat_facade

"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test:

        * the run of a simulation with/without interpolation of the forcings,
        * the logging,
        * the postprocessing,
        * and the graphs generation.

"""

# -- SIMULATION PARAMETERS --
# Length of the simulation (in hours)
SIMULATION_LENGTH = 5

# define the time step in hours for each simulator
CARIBU_TIMESTEP = 4
SENESCWHEAT_TIMESTEP = 1
FARQUHARWHEAT_TIMESTEP = 1
ELONGWHEAT_TIMESTEP = 1
GROWTHWHEAT_TIMESTEP = 1
CNWHEAT_TIMESTEP = 1

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


def test_cnwheat():
    # -- CNWHEAT --
    # Initial states
    cnwheat_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   {},
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df,
                                                   update_shared_df=UPDATE_SHARED_DF)

    # Run facade
    cnwheat_facade_.run(12, 10, {'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5})


def test_elongwheat():
    # Initial states
    elongwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS + [i for i in
                                                                   elongwheat_facade.simulation.HIDDENZONE_INPUTS
                                                                   if i in
                                                                   inputs_dataframes[
                                                                       HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.ELEMENT_INPUTS
                                                                if
                                                                i in
                                                                inputs_dataframes[
                                                                    ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.AXIS_INPUTS if
                                                             i in inputs_dataframes[
                                                                 AXES_INITIAL_STATE_FILENAME].columns]].copy()

    phytoT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'phytoT.csv'))

    # Facade initialisation
    elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                            ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            elongwheat_axes_initial_state,
                                                            elongwheat_hiddenzones_initial_state,
                                                            elongwheat_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT_df,
                                                            None,
                                                            update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    elongwheat_facade_.run(12, 10)
    adel_wheat.update_geometry(g)


def test_caribu():
    caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                shared_elements_inputs_outputs_df,
                                                adel_wheat,
                                                update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    caribu_facade_.run(True, energy=1000, DOY=1, hourTU=12, latitude=48.85, sun_sky_option='sky',
                       heterogeneous_canopy=False, plant_density=PLANT_DENSITY[1])

def test_farquharwheat():
    # adhoc modification for test
    g.property('diameter').update({19: 0.003, 34: 0.003})
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
    update_parameters_farquharwheat = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df,
                                                                     update_parameters_farquharwheat,
                                                                     update_shared_df=UPDATE_SHARED_DF)

    # Run facade
    farquharwheat_facade_.run(12, 400, 0.8, 1.5)

def test_growthwheat():
    # Initial states
    growthwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ELEMENT_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ROOT_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(g,
                                                               GROWTHWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               growthwheat_hiddenzones_initial_state,
                                                               growthwheat_elements_initial_state,
                                                               growthwheat_root_initial_state,
                                                               growthwheat_axes_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_hiddenzones_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               None,
                                                               update_shared_df=UPDATE_SHARED_DF)
        # Run facade
    growthwheat_facade_.run()

def test_senescwheat():
    # adhoc modification for test
    g.property('delta_teq').update({2: 3600})
    g.property('delta_teq_roots').update({2: 3600})
    g.property('sum_TT').update({2: 0})

    # Initial states
    senescwheat_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[
        inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS if
         i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS if
         i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS if
         i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                               SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescwheat_roots_initial_state,
                                                               senescwheat_axes_initial_state,
                                                               senescwheat_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               None,
                                                               update_shared_df=UPDATE_SHARED_DF)
    # Run facade
    senescwheat_facade_.run()

def test_fspmwheat():
    fspmwheat_facade.FSPMWheatFacade(g)