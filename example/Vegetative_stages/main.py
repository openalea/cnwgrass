from openalea.cnwgrass.integration.runner import run as runner


"""
    main
    ~~~~

    This example simulates a field experiment on winter wheat (cv Soissons) described in Ljutovac (2002).
    The simulation starts at leaf 4 emergence on December 1998 and finishes on April 1999 at the beginning of internode elongation.
    No coupling with hydraulics.
    Simulation results are described in Gauthier et al (2020).https://doi.org/10.1093/jxb/eraa276.

"""

simulation_length = 2500  # hours
METEO_FILENAME = 'meteo_Ljutovac2002.csv'


if __name__ == '__main__':
    runner(simulation_length=simulation_length, forced_start_time=0,
                    run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False,
                    tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}, heterogeneous_canopy=True,
                    METEO_FILENAME=METEO_FILENAME, MANAGEMENT_FILENAME='management.csv')
