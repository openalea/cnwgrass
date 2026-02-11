import os
import pandas as pd

# Get the list of scenarios
subdirectories_outputs = os.listdir('outputs')
scenarios = []
for i in subdirectories_outputs:
    if i[:9] == 'Scenario_':
        scenarios.append(int(i[9:]))
scenarios_inputs = scenarios

INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']


def rearrange_postprocessing(postprocessing_tables, t=None, scenarios=None, merge_with_scenarios_list=True, scenarios_list_columns=None,
                             outputs_dir_path='outputs'):
    """
    For each graph type, create a common directory with the graphs of all the scenarios.

    :param list postprocessing_tables: List of names of postprocessing tables that we want to compare for each scenario.
    :param int t: Time step of the model for which we want to compare the postprocessing of each scenario.
    :param list scenarios: List of scenarios (numbers) for which the graphs will be rearranged.
    :param bool merge_with_scenarios_list: if True, the output datafram will be merge with a subset of the scenarios_list table (only the columns in scenarios_list_columns)
    :param list scenarios_list_columns: List of columns name of the scenarios_list table that are of interest.
    :param str outputs_dir_path: the path with the outputs of the scenarios
    """

    if scenarios is None:
        scenarios = scenarios_inputs

    # Create directory with general postprocessing
    general_res_dir = os.path.join(outputs_dir_path, 'General_postprocessing')
    if not os.path.exists(general_res_dir):
        os.mkdir(general_res_dir)

    # read the description of the scenarios if necessary
    if merge_with_scenarios_list:
        if scenarios_list_columns is None:
            scenarios_list_columns = ['Scenario', 'Scenario_label', 'Plant_Density',
                                      'constant_Conc_Nitrates', 'fertilization_U_3500', 'fertilization_quantity', 'fertilization_interval',
                                      'PAR', 'Inputs_PlantSoil_Dirpath']
        scenarios_df = pd.read_csv(os.path.join('inputs', 'scenarios_list.csv'), index_col='Scenario')
        scenarios_df['Scenario'] = scenarios_df.index

        scenarios_list_columns = [c for c in scenarios_list_columns if c in scenarios_df.columns] + [c for c in scenarios_df.columns if c.startswith('P_')]

    #: For each postprocessing table, extract the prostprocessings at t of the scenarios
    for pp_table in postprocessing_tables:
        res = []
        res_index_column = []
        res_file_name = None
        scenarios_with_res = []

        #: For each scenario, extract the prostprocessings at t
        for scenario in scenarios:

            scenario_name = 'Scenario_%.4d' % scenario
            scenario_postprocessing_dirpath = os.path.join(outputs_dir_path, scenario_name, 'postprocessing')

            pp_table_name = pp_table
            if pp_table == 'Conc_phloem':
                pp_table_name = 'organs_postprocessing'
            scenario_postprocessing_table_dirpath = os.path.join(scenario_postprocessing_dirpath, pp_table_name + '.csv')

            if not os.path.exists(scenario_postprocessing_table_dirpath):
                continue
            pp_res = pd.read_csv(scenario_postprocessing_table_dirpath)
            if pp_table == 'Conc_phloem':
                pp_res = pp_res[pp_res.organ == 'phloem'][['t', 'Conc_Amino_Acids', 'Conc_Sucrose']].copy()
                pp_res.reset_index(drop=True, inplace=True)
            if (t is not None) and ('t' in pp_res.columns):
                if t not in pp_res.t:
                    continue
                scenarios_with_res.append(scenario)
                res.append(pp_res[pp_res.t == t])
                if not res_file_name:
                    res_file_name = pp_table + '_' + str(int(t))
            else:
                scenarios_with_res.append(scenario)
                res.append(pp_res)
                if not res_file_name:
                    res_file_name = pp_table

            # Find the index of the postprocessing table
            if not res_index_column:
                res_index_column = ['Scenario']
                res_index_column.extend([c for c in INDEX_COLUMNS if c in pp_res.columns])

        # convert list of postprocessing into datafram
        res_df = pd.concat(res, keys=scenarios_with_res, sort=False)
        res_df.reset_index(0, inplace=True)
        res_df.rename({'level_0': 'Scenario'}, axis=1, inplace=True)
        res_df = res_df.reindex(res_index_column + res_df.columns.difference(res_index_column).tolist(), axis=1, copy=False)

        # Add info about the scenario
        if merge_with_scenarios_list:
            scenarios_description = scenarios_df[scenarios_list_columns]
            scenarios_description.reset_index(drop=True, inplace=True)
            res_df.reset_index(drop=True, inplace=True)
            res_df = scenarios_description.merge(res_df, how='outer', on='Scenario')

        # write the datafram
        res_dir = os.path.join(general_res_dir, res_file_name + '.csv')
        res_df.to_csv(res_dir, na_rep='NA', index=False)


if __name__ == '__main__':
    rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=1999)
    rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=3499)
    rearrange_postprocessing(postprocessing_tables=['organs_postprocessing'], t=3499)
    rearrange_postprocessing(postprocessing_tables=['hiddenzones_postprocessing'], t=3499)
    rearrange_postprocessing(postprocessing_tables=['performance_indices'])
    rearrange_postprocessing(postprocessing_tables=['leaf_traits'])
    rearrange_postprocessing(postprocessing_tables=['canopy_dynamics_daily'])
    rearrange_postprocessing(postprocessing_tables=['Conc_phloem'])
    rearrange_postprocessing(postprocessing_tables=['RUE'])
