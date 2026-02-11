import multiprocessing as mp
import os
import time

import pandas as pd
from example.Scenarios_monoculms import rearrange_graphs
from example.Scenarios_monoculms import rearrange_postprocessing
from example.Scenarios_monoculms import run_fspmwheat

scenarios_df = pd.read_csv(os.path.join('inputs', 'scenarios_list.csv'), index_col='Scenario')
scenarios_df['Scenario'] = scenarios_df.index
scenarios = scenarios_df.Scenario

if __name__ == '__main__':
    tstart = time.time()
    num_processes = mp.cpu_count()
    p = mp.Pool(num_processes)

    mp_solutions = p.map(run_fspmwheat.run_fspmwheat, list(scenarios))
    p.terminate()
    p.join()

    if 'Generate_Graphs' in scenarios_df.columns and any(scenarios_df.Generate_Graphs):
        rearrange_graphs.rearrange_graphs(scenarios=list(scenarios))
    if 'Run_Postprocessing' in scenarios_df.columns and any(scenarios_df.Run_Postprocessing):
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=1999, scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=3499, scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['organs_postprocessing'], t=3499, scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['performance_indices'], scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['leaf_traits'], scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['canopy_dynamics_daily'], scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['Conc_phloem'], scenarios=list(scenarios))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['RUE'], scenarios=list(scenarios))

    tend = time.time()
    tmp = (tend - tstart) / 60.
    print("multiprocessing: %8.3f minutes" % tmp)
