import os
import pandas as pd

from openalea.cnwgrass.integration.runner import run as runner

"""
    main
    ~~~~

    This example illustrates the coupling of CNW-Grass with an external soil model, which provides root N uptake as an input.
    The initial conditions and weather data are the same as "Vegetative_stages".
    No coupling with hydraulics.

"""

simulation_length = 2500  # hours
METEO_FILENAME = 'meteo_Ljutovac2002.csv'

RERmax_vegetative_stages_example = {'morphogenesis': {'RERmax': {5: 3.35e-06, 6: 2.1e-06, 7: 2.e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}}}

# get nitrate uptake data
NITRATES_UPTAKE_FORCINGS_FILENAME = 'nitrates_uptake_forcings.csv'
nitrates_uptake_data_filepath = os.path.join('inputs', NITRATES_UPTAKE_FORCINGS_FILENAME)
nitrates_uptake_data_df = pd.read_csv(nitrates_uptake_data_filepath)
nitrates_uptake_data_grouped = nitrates_uptake_data_df.groupby( ['t', 'plant', 'axis', 'organ'])

def force_nitrates_uptake(t, population, g):
    """Force the nitrates uptake data of the population at `t` from input grouped dataframes

    :param int t: the current time step (s)
    :param cnmetabolism.model.Population population: the population of organs and elements provided by the cnmetabolism model
    :param openalea.mtg.mtg.MTG g: the MTG

    """
    mtg_plants_iterator = g.components_iter(g.root)
    for plant in population.plants:
        cnmetabolism_plant_index = plant.index
        while True:
            mtg_plant_vid = next(mtg_plants_iterator)
            if int(g.index(mtg_plant_vid)) == cnmetabolism_plant_index:
                break
        mtg_axes_iterator = g.components_iter(mtg_plant_vid)
        for axis in plant.axes:
            # Update Nitrates uptake in dataframe
            group = nitrates_uptake_data_grouped.get_group((t, plant.index, axis.label, 'roots'))
            nitrates_uptake_data_to_use = group.loc[
                group.first_valid_index(), group.columns.intersection(['Uptake_Nitrates'])].dropna().to_dict()
            axis.roots.__dict__.update(nitrates_uptake_data_to_use)

            # Update Nitrates uptake in MTG
            cnmetabolism_axis_label = axis.label
            while True:
                mtg_axis_vid = next(mtg_axes_iterator)
                if g.label(mtg_axis_vid) == cnmetabolism_axis_label:
                    break
            mtg_roots_properties = g.get_vertex_property(mtg_axis_vid)['roots']
            mtg_roots_properties.update(nitrates_uptake_data_to_use)



if __name__ == '__main__':
    runner(simulation_length=simulation_length, forced_start_time=0,
                    run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False,
                    tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}, heterogeneous_canopy=True,
                    external_soil_model=True, step_callback={'nitrate_uptake': force_nitrates_uptake},
                    update_parameters_all_models=RERmax_vegetative_stages_example,
                    METEO_FILENAME=METEO_FILENAME, MANAGEMENT_FILENAME='management.csv')
