from openalea.fspmwheat.fspmwheat_runner import run as fspmwheat_runner


"""
    main
    ~~~~

    This example simulates a field experiment on winter wheat (cv Soissons) described in Ljutovac (2002).
    The simulation starts at leaf 4 emergence on December 1998 and finishes on April 1999 at the beginning of internode elongation.
    This example integrates the co-regulation of leaf growth by the trophic and hydraulic status.
    Simulation results are described in Acker et al (2026).https://doi.org/10.1093/jxb/erag248.  

"""

simulation_length = 2500  # hours
METEO_FILENAME = 'meteo_Ljutovac2002.csv'

# Drought and rehydration scenario
drought_trigger = None#{'green_area': 0.002388}                   # plant green area at which the drought treatment starts (m2)
rehydration_scenario = None#{'stop_drought_SRWC': 20.,      # SRWC at which the drought event stops (%)
                         # 'SRWC_target': 80.,            # Target SRWC for rehydration
                         # 'rehydration_duration': 15.}   # duration of the rehydration period (days)



if __name__ == '__main__':
    fspmwheat_runner(simulation_length=simulation_length, forced_start_time=0,
                    run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False,
                    hydraulics=True, stomatal_model_name='hydraulics', drought_trigger=drought_trigger, rehydration_scenario=rehydration_scenario,
                    tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}, heterogeneous_canopy=True,
                    METEO_FILENAME=METEO_FILENAME, MANAGEMENT_FILENAME='management.csv')
