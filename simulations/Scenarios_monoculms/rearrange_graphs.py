import os
import shutil

# Graphs to be rearranged
graphs_titles_inputs = ['N_content_roots_axis', 'N_content_shoot_axis', 'N_content_blade', 'shoot_roots_ratio_axis', 'green_area_blade', 'Ag_blade',
                        'lamina_Wmax', 'lamina_W_Lg', 'leaf_Lmax', 'internode_Lmax', 'sheath_Lmax', 'phyllochron', 'Conc_Sucrose_phloem', 'Conc_Amino_acids_phloem', 'S_Amino_Acids_roots', 'LAI', 'RER_comparison',
                        'Conc_Nitrates_roots', 'Conc_Sucrose_roots', 'sum_dry_mass_roots_axis', 'sum_dry_mass_shoot_axis', 'sum_dry_mass_axis', 'Summary', 'SSLW','LSSW',
                        'RUE', 'SLN_blade','sum_N_g_axis',
                        'SLN_nonstruct_blade', 'Conc_Nitrates_Soil',
                        'C_usages_cumulated', 'Photosynthetic_yield_blade', 'Photosynthetic_yield_plante', 'Photosynthetic_yield_plante2', 'Ts_blade', 'internode_Lmax', 'NNI_axis', 'N_dilution',
                        'S_Amino_Acids_blade', 'HATS_LATS_roots', 'senesced_mstruct_roots', 'synthetized_mstruct_roots', 'mstruct_axis', 'mstruct_roots',
                        'Cont_WSC_DM_axis', 'Cont_WSC_DM_hz','Cont_WSC_DM_blade','Cont_WSC_DM_sheath', 'Photosynthetic_rates_F6',
                        'Transpiration_blade']

# Get the list of scenarios
titi = os.listdir('outputs')
scenarios = []
for i in titi:
    if i[:9] == 'Scenario_':
        scenarios.append(int(i[9:]))
scenarios_inputs = scenarios


def rearrange_graphs(graphs_titles=None, scenarios=None, outputs_dir_path=None):
    """
    For each graph type, create a common directory with the graphs of all the scenarios.

    :param list graphs_titles: List of graph title to be rearranged.
    :param list scenarios: List of scenarios (numbers) for which the graphs will be rearranged.
    :param str outputs_dir_path: the path with the outputs of the scenarios
    """

    if scenarios is None:
        scenarios = scenarios_inputs
    if graphs_titles is None:
        graphs_titles = graphs_titles_inputs
    if not outputs_dir_path:
        outputs_dir_path = 'outputs'

    # Create directory with graphs by type
    graph_dir = os.path.join(outputs_dir_path, 'Graphs_by_type')
    if not os.path.exists(graph_dir):
        os.mkdir(graph_dir)

    #: For each graph type, create a common directory
    for graph in graphs_titles:

        graph_dir = os.path.join(outputs_dir_path, 'Graphs_by_type', graph)
        if not os.path.exists(graph_dir):
            os.mkdir(graph_dir)

        #: For each scenario, copy and rename the graphs in the common directory
        for scenario in scenarios:

            scenario_name = 'Scenario_%.4d' % scenario

            scenario_dir = os.path.join(outputs_dir_path, scenario_name, 'graphs')
            graph_src = os.path.join(scenario_dir, graph + '.PNG')
            graph_dest = os.path.join(graph_dir, graph + '.PNG')
            graph_renamed = os.path.join(graph_dir, scenario_name + '.PNG')  # os.path.join(graph_dir, scenario_name+'_'+scenario_label+'.PNG')

            if os.path.exists(graph_renamed):
                os.remove(graph_renamed)
            if os.path.exists(graph_src):
                shutil.copyfile(graph_src, graph_dest)
                os.rename(graph_dest, graph_renamed)


if __name__ == '__main__':
    rearrange_graphs()
