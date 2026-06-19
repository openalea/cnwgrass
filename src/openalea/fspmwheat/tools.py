# -*- coding: latin-1 -*-

import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
import matplotlib.ticker as mtick
import statsmodels.api as sm

from openalea.adel.mtg import to_plantgl
from openalea.plantgl.all import Viewer, Vector3

"""
    fspmwheat.tools
    ~~~~~~~~~~~~~~~

    This module provides convenient tools needed by the facades.

"""

OUTPUTS_INDEXES = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']  #: All the possible indexes

def combine_dataframes_inplace(model_dataframe, shared_column_indexes, shared_dataframe_to_update):
    """Combine `model_dataframe` and `shared_dataframe_to_update` in-place:

           * re-index `model_dataframe` and `shared_dataframe_to_update` by `shared_column_indexes`,
           * use method pd.DataFrame.combine_first(),
           * reset to the right types in `shared_dataframe_to_update`,
           * reorder the columns: first columns in `shared_column_indexes`, then others columns alphabetically,
           * and reset the index in `shared_dataframe_to_update`.

    :param pandas.DataFrame model_dataframe: dataframe to use for updating `shared_dataframe_to_update`.
    :param list shared_column_indexes: The indexes to re-index `model_dataframe` and `shared_dataframe_to_update` before combining them.
    :param pandas.DataFrame shared_dataframe_to_update: The dataframe to update.

    note:: `shared_dataframe_to_update` is updated in-place. Thus, `shared_dataframe_to_update` keeps the same object's memory address.

    """

    # re-index the dataframes to have common indexes
    if len(shared_dataframe_to_update) == 0:
        shared_dataframe_to_update_reindexed = shared_dataframe_to_update
    else:
        shared_dataframe_to_update.sort_values(shared_column_indexes, inplace=True)
        shared_dataframe_to_update_reindexed = pd.DataFrame(shared_dataframe_to_update.values.tolist(),
                                                            index=sorted(shared_dataframe_to_update.groupby(shared_column_indexes).groups.keys()),
                                                            columns=shared_dataframe_to_update.columns)

    model_dataframe.sort_values(shared_column_indexes, inplace=True)
    model_dataframe_reindexed = pd.DataFrame(model_dataframe.values.tolist(),
                                             index=sorted(model_dataframe.groupby(shared_column_indexes).groups.keys()),
                                             columns=model_dataframe.columns)

    # combine model and shared re-indexed dataframes
    if model_dataframe_reindexed.empty and shared_dataframe_to_update.empty:
        new_shared_dataframe = model_dataframe_reindexed.copy()
        for new_header in shared_dataframe_to_update_reindexed.columns.difference(model_dataframe_reindexed.columns):
            new_shared_dataframe[new_header] = ""
    else:
        new_shared_dataframe = model_dataframe_reindexed.combine_first(shared_dataframe_to_update_reindexed)

    # reset to the right types in the combined dataframe
    dtypes = model_dataframe_reindexed.dtypes.combine_first(shared_dataframe_to_update_reindexed.dtypes)
    for column_name, data_type in dtypes.items():
        if pd.api.types.is_integer_dtype(data_type) and new_shared_dataframe[column_name].isnull().values.any():  # Used to keep bool values
            data_type = float  # will return an error if data_type is integer
        new_shared_dataframe[column_name] = new_shared_dataframe[column_name].astype(data_type)

    # reorder the columns
    new_shared_dataframe = new_shared_dataframe.reindex(shared_column_indexes + sorted(new_shared_dataframe.columns.difference(shared_column_indexes)), axis=1)

    # update the shared dataframe in-place
    shared_dataframe_to_update.drop(shared_dataframe_to_update.index, axis=0, inplace=True)
    shared_dataframe_to_update.drop(shared_dataframe_to_update.columns, axis=1, inplace=True)
    shared_dataframe_to_update['dataframe_to_update_index'] = new_shared_dataframe.index
    shared_dataframe_to_update.set_index('dataframe_to_update_index', inplace=True)
    for column in new_shared_dataframe.columns:
        shared_dataframe_to_update[column] = new_shared_dataframe[column]
    shared_dataframe_to_update.reset_index(0, drop=True, inplace=True)
    

def plot_outputs(outputs, x_name, y_name, x_label='', y_label='', x_lim=None, title=None, meteo_data=None, filters={}, plot_filepath=None,
                        colors=[], linestyles=[], explicit_label=True, kwargs={}):
    """Plot `outputs`, with x=`x_name` and y=`y_name`.

    The general algorithm is:

        * find the scale of `outputs` and keep only the needed columns,
        * apply `filters` to `outputs` and make groups according to the scale,
        * plot each group as a new line,
        * save or display the plot.

    :param pandas.DataFrame outputs: The outputs of CN-Wheat.
    :param str x_name: x-axis of the plot.
    :param str y_name: y-axis of the plot.
    :param str x_label: The x label of the plot. Default is ''.
    :param str or unicode y_label: The y label of the plot. Default is ''.
    :param float x_lim: the x-axis limit.
    :param str title: the title of the plot. If None (default), create a title which is the concatenation of `y_name` and each scales which cardinality is one.
    :param pandas.DataFrame meteo_data: the meteo dataframe having the mapping between t (hours) and calendar dates
    :param dict filters: A dictionary whose keys are the columns of `outputs` for which we want to apply a specific filter.
          These columns can be one or more element of :const:`OUTPUTS_INDEXES`.
          The value associated to each key is a criteria that the rows of `outputs`
          must satisfy to be plotted. The values can be either one value or a list of values.
          If no value is given for any column, then all rows are plotted (default).
    :param list colors: The colors for lines. If empty, let matplotlib default line colors.
    :param list linestyles: The styles for lines. If empty, let matplotlib default line styles.
    :param str plot_filepath: The file path to save the plot. If `None`, do not save the plot but display it.
    :param bool explicit_label: True: makes the line label from concatenation of each scale id (default).
                              - False: makes the line label from concatenation of scales containing several distinct elements.
    :param dict kwargs: key arguments to be passed to matplolib

    Examples::

            import pandas as pd
            cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # in this example, 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose'.
            plot(cnwheat_output_df, x_name = 't', y_name = 'Conc_Sucrose', x_label='Time (Hour)', y_label=u'[Sucrose] (µmol g$^{-1}$ mstruct)', title='{} = f({})'.format('Conc_Sucrose', 't'), filters={'plant': 1, 'axis': 'MS', 'organ': 'Lamina', 'element': 1})

    """

    # finds the scale of `outputs`
    group_keys = [key for key in OUTPUTS_INDEXES if key in outputs and key != x_name and key != y_name]

    # make a group_keys with first letter of each key in upper case
    group_keys_upper = [group_key[0].upper() + group_key[1:] for group_key in group_keys]

    # create a mapping to associate each key to its index in group_keys
    group_keys_mapping = dict([(key, index) for (index, key) in enumerate(group_keys)])

    # keep only the needed columns (to make the grouping faster)
    outputs = outputs[group_keys + [x_name, y_name]]

    # apply filters to outputs
    for key, value in filters.items():
        if key in outputs:
            # convert to list if needed
            try:
                _ = iter(value)
            except TypeError:
                values = [value]
            else:
                values = value
                # handle strings too
                if isinstance(values, str):
                    values = [values]
            # select data from outputs
            outputs = outputs[outputs[key].isin(values)]

    # do not plot if there is nothing to plot
    if outputs[y_name].isnull().all():
        return

    # compute the cardinality of each group keys and create the title if needed
    subtitle_groups = []
    labels_groups = []
    for i in range(len(group_keys)):
        group_key = group_keys[i]
        group_cardinality = outputs[group_key].nunique()
        if group_cardinality == 1:
            group_value = outputs[group_key][outputs.first_valid_index()]
            subtitle_groups.append('{}: {}'.format(group_keys_upper[i], group_value))
        else:
            labels_groups.append(group_key)
    if title is None:  # we need to create the title
        title = y_name + '\n' + ' - '.join(subtitle_groups)

    # makes groups according to the scale
    outputs_grouped = outputs.groupby(group_keys)

    # plots each group as a new line
    fig, ax = plt.subplots()

    matplot_colors_cycler = cycle(colors)
    matplot_linestyles_cycler = cycle(linestyles)

    for outputs_group_name, outputs_group in outputs_grouped:
        line_label_list = []
        if explicit_label:
            # concatenate the keys of the group name
            line_label_list.extend(['{}: {}'.format(group_keys_upper[group_keys_mapping[output_group_name]], outputs_group_name) for output_group_name in outputs_group_name])
        else:
            # construct a label with only the essential keys of the group name ; the essential keys are those for which cardinality is non-zero
            for label_group in labels_groups:
                label_group_index = group_keys_mapping[label_group]
                line_label_list.append('{}: {}'.format(group_keys_upper[label_group_index], outputs_group_name[label_group_index]))

        kwargs['label'] = ' - '.join(line_label_list)

        # apply user colors
        try:
            color = next(matplot_colors_cycler)
        except StopIteration:
            pass
        else:
            kwargs['color'] = color

        # apply user lines style
        try:
            linestyle = next(matplot_linestyles_cycler)
        except StopIteration:
            pass
        else:
            kwargs['linestyle'] = linestyle

        # plot the line
        ax.plot(outputs_group[x_name], outputs_group[y_name], **kwargs)

    if y_name not in ('water_potential', 'osmotic_water_potential'):
        ax.set_ylim(bottom=0.)

    if x_lim is not None:
        ax.set_xlim(left=0, right=x_lim)
    else:
        ax.set_xlim(left=0)

    if meteo_data is not None:
        meteo_data['Date'] = pd.to_datetime(meteo_data['Date'], format='%d/%m/%Y')
        ax2 = ax.twiny()
        ax2.set_xticks(ax.get_xticks())
        ax2.set_xticklabels(meteo_data.loc[ax.get_xticks()]['Date'].dt.strftime('%d/%m'))
        ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
        ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
        ax2.spines['bottom'].set_position(('outward', 35))

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if kwargs['label']:
        ax.legend(prop={'size': 6}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
    ax.set_title(title)
    plt.tight_layout()

    if plot_filepath is None:
        # display the plot
        plt.show()
    else:
        # save the plot
        plt.savefig(plot_filepath, dpi=200, format='PNG', bbox_inches='tight')
        plt.close()


def additional_graphs(axes_outputs, hz_outputs, elements_outputs,
                      axes_postprocessing, hz_postprocessing, elements_postprocessing, organs_postprocessing,
                      plant_density, RER_max_param, GRAPHS_DIRPATH, data_obs):

    """
    :param pandas.DataFrame axes_outputs:
    :param pandas hz_outputs:
    :param pandas.DataFrame elements_outputs:
    :param pandas.DataFrame axes_postprocessing:
    :param pandas.DataFrame hz_postprocessing:
    :param pandas.DataFrame elements_postprocessing:
    :param pandas.DataFrame organs_postprocessing:
    :param int plant_density:
    :param dict RER_max_param:
    :param str GRAPHS_DIRPATH:
    :param pandas.DataFrame data_obs:
    """

    colors = ['blue', 'darkorange', 'green', 'red', 'darkviolet', 'gold', 'magenta', 'brown', 'darkcyan', 'grey', 'lime']
    colors = colors + colors

    # 0) Phyllochron
    df_MS = axes_outputs[axes_outputs['axis'] == 'MS']

    grouped_df = hz_postprocessing[hz_postprocessing['axis'] == 'MS'].groupby(['plant', 'metamer'])[['t', 'leaf_is_emerged']]
    grouped_iter = iter(grouped_df)
    next(grouped_iter, None)  # ignoring the first leaf as we can't calculate its phyllochron
    leaf_emergence = {}
    for group_name, data in grouped_iter:
        plant, metamer = group_name[0], group_name[1]
        if True not in data['leaf_is_emerged'].unique():
            continue
        leaf_emergence_t = data[data['leaf_is_emerged'] == True].iloc[0]['t']
        leaf_emergence[(plant, metamer)] = leaf_emergence_t

    phyllochron = {'plant': [], 'metamer': [], 'phyllochron': []}
    for key, leaf_emergence_t in sorted(leaf_emergence.items()):
        plant, metamer = key[0], key[1]
        if (plant, metamer - 1) not in leaf_emergence.keys():
            continue
        phyllochron['plant'].append(plant)
        phyllochron['metamer'].append(metamer)
        prev_leaf_emergence_t = leaf_emergence[(plant, metamer - 1)]
        if df_MS[(df_MS['t'] == leaf_emergence_t) | (df_MS['t'] == prev_leaf_emergence_t)].sum_TT.count() == 2:
            phyllo_DD = df_MS[(df_MS['t'] == leaf_emergence_t)].sum_TT.values[0] - df_MS[(df_MS['t'] == prev_leaf_emergence_t)].sum_TT.values[0]
        else:
            phyllo_DD = np.nan
        phyllochron['phyllochron'].append(phyllo_DD)

    if len(phyllochron['metamer']) > 0:
        fig, ax = plt.subplots()
        plt.xlim((int(min(phyllochron['metamer']) - 1), int(max(phyllochron['metamer']) + 1)))
        plt.ylim(ymin=0, ymax=150)
        ax.plot(phyllochron['metamer'], phyllochron['phyllochron'], color='b', marker='o')
        for i, j in zip(phyllochron['metamer'], phyllochron['phyllochron']):
            ax.annotate(str(int(round(j, 0))), xy=(i, j + 2), ha='center')
        ax.set_xlabel('Leaf number')
        ax.set_ylabel('Phyllochron (Degree Day)')
        ax.set_title('phyllochron')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'phyllochron' + '.PNG'))
        plt.close()

    # 1) Comparison Dimensions with Ljutovac 2002
    df_hz_filtered = hz_outputs[(hz_outputs['axis'] == 'MS') & (hz_outputs['plant'] == 1) & ~np.isnan(hz_outputs.leaf_Lmax)].copy()
    df_IN = df_hz_filtered[~ np.isnan(df_hz_filtered.internode_Lmax)]
    last_value_idx = df_hz_filtered.groupby(['metamer'])['t'].transform(max) == df_hz_filtered['t']
    df_hz_filtered_end = df_hz_filtered[last_value_idx].copy()
    df_hz_filtered_end['lamina_Wmax'] = df_hz_filtered_end.leaf_Wmax
    df_hz_filtered_end['lamina_W_Lg'] = df_hz_filtered_end.leaf_Wmax / df_hz_filtered_end.lamina_Lmax
    bchmk = data_obs.loc[data_obs.metamer >= min(df_hz_filtered_end.metamer)]
    bchmk['lamina_W_Lg'] = bchmk.lamina_Wmax / bchmk.lamina_Lmax
    last_value_idx = df_IN.groupby(['metamer'])['t'].transform(max) == df_IN['t']
    df_IN_last = df_IN[last_value_idx].copy()
    res = df_hz_filtered_end[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'lamina_W_Lg', 'SSLW', 'LSSW']].merge(df_IN_last[['metamer', 'internode_Lmax']], left_on='metamer',
                                                                                                                          right_on='metamer', how='outer').copy()

    var_list = ['leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'internode_Lmax']
    for var in list(var_list):
        fig, ax = plt.subplots()
        plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
        plt.ylim(ymin=0, ymax=np.nanmax(list(res[var] * 100 * 1.05) + list(bchmk[var] * 1.05)))

        tmp = res[['metamer', var]].drop_duplicates().sort_values('metamer').reset_index(drop=True)

        line1 = ax.plot(tmp.metamer, tmp[var] * 100, color='c', marker='o')
        line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')

        ax.set_ylabel(var + ' (cm)')
        ax.set_title(var)
        ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
        plt.close()

    var = 'lamina_W_Lg'
    fig, ax = plt.subplots()
    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
    plt.ylim(ymin=0, ymax=np.nanmax(list(res[var] * 1.05) + list(bchmk[var] * 1.05)))
    tmp = res[['metamer', var]].drop_duplicates().sort_values('metamer').reset_index(drop=True)
    line1 = ax.plot(tmp.metamer, tmp[var], color='c', marker='o')
    line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')
    ax.set_ylabel(var)
    ax.set_title(var)
    ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
    plt.close()

    # 1bis) Comparison Structural Masses vs. adaptation from Bertheloot 2008

    # SSLW Laminae
    bchmk = pd.DataFrame.from_dict({1: 15, 2: 23, 3: 25, 4: 18, 5: 22, 6: 25, 7: 20, 8: 23, 9: 26, 10: 28, 11: 31}, orient='index').rename(columns={0: 'SSLW'})
    bchmk.index.name = 'metamer'
    bchmk = bchmk.reset_index()
    bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

    fig, ax = plt.subplots()
    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
    plt.ylim(ymin=0, ymax=50)

    tmp = res[['metamer', 'SSLW']].drop_duplicates().sort_values('metamer').reset_index(drop=True)

    line1 = ax.plot(tmp.metamer, tmp.SSLW, color='c', marker='o')
    line2 = ax.plot(bchmk.metamer, bchmk.SSLW, color='orange', marker='o')

    ax.set_ylabel('Structural Specific Lamina Weight (g.m-2)')
    ax.set_title('Structural Specific Lamina Weight')
    ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=3)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SSLW.PNG'))
    plt.close()

    # LWS Sheaths
    bchmk = pd.DataFrame.from_dict({1: 0.08, 2: 0.09, 3: 0.11, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55, 11: 0.65}, orient='index').rename(columns={0: 'LSSW'})
    bchmk.index.name = 'metamer'
    bchmk = bchmk.reset_index()
    bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

    fig, ax = plt.subplots()
    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
    plt.ylim(ymin=0, ymax=0.8)

    tmp = res[['metamer', 'LSSW']].drop_duplicates().sort_values('metamer').reset_index(drop=True)

    line1 = ax.plot(tmp.metamer, tmp.LSSW, color='c', marker='o')
    line2 = ax.plot(bchmk.metamer, bchmk.LSSW, color='orange', marker='o')

    ax.set_ylabel('Lineic Structural Sheath Weight (g.m-1)')
    ax.set_title('Lineic Structural Sheath Weight')
    ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=2)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'LSSW.PNG'))
    plt.close()

    # 2) LAI
    elements_postprocessing['green_area_rep'] = elements_postprocessing.green_area * elements_postprocessing.nb_replications
    grouped_df = elements_postprocessing[(elements_postprocessing.axis == 'MS') & (elements_postprocessing.element == 'LeafElement1')].groupby(['t', 'plant'])
    LAI_dict = {'t': [], 'plant': [], 'LAI': []}
    for name, data in grouped_df:
        t, plant = name[0], name[1]
        LAI_dict['t'].append(t)
        LAI_dict['plant'].append(plant)
        LAI_dict['LAI'].append(data['green_area_rep'].sum() * plant_density)

    plot_outputs(pd.DataFrame(LAI_dict), 't', 'LAI', x_label='Time (Hour)', y_label='LAI', plot_filepath=os.path.join(GRAPHS_DIRPATH, 'LAI.PNG'), explicit_label=False)

    # 3) RER during the exponentiel-like phase

    # - RER parameters
    rer_param = dict((k, v) for k, v in RER_max_param)

    # - Simulated RER

    # import simulation outputs
    data_RER = hz_postprocessing[(hz_postprocessing.axis == 'MS') & (hz_postprocessing.metamer >= 1)].copy()
    data_RER.sort_values(['t', 'metamer'], inplace=True)

    # - Time previous leaf emergence
    tmp = data_RER[data_RER.leaf_is_emerged]
    leaf_em = tmp.groupby('metamer', as_index=False)['t'].min()
    leaf_em['t_em'] = leaf_em.t
    prev_leaf_em = leaf_em
    prev_leaf_em.metamer = leaf_em.metamer + 1

    data_RER2 = pd.merge(data_RER, prev_leaf_em[['metamer', 't_em']], on='metamer')
    data_RER2 = data_RER2[data_RER2.t <= data_RER2.t_em]

    # - SumTimeEq
    df_MS['SumTimeEq'] = np.cumsum(df_MS.delta_teq)
    data_RER3 = pd.merge(data_RER2, df_MS[['t', 'SumTimeEq']], on='t')

    # - logL
    data_RER3['logL'] = np.log(data_RER3.leaf_L)

    # - Estimate RER
    RER_sim = {}
    for leaf in data_RER3.metamer.drop_duplicates():
        Y = data_RER3.logL[data_RER3.metamer == leaf]
        X = data_RER3.SumTimeEq[data_RER3.metamer == leaf]
        X = sm.add_constant(X)
        mod = sm.OLS(Y, X)
        fit_RER = mod.fit()
        RER_sim[leaf] = fit_RER.params['SumTimeEq']
    if len(RER_sim) != 0:
        # - Graph
        fig, ax1 = plt.subplots()
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        x, y = zip(*sorted(RER_sim.items()))
        ax1.plot(x, y, label=r'Simulated RER', linestyle='-', color='g')
        ax1.errorbar(data_obs.metamer, data_obs.RER, yerr=data_obs.RER_confint, marker='o', color='g', linestyle='', label="Observed RER", markersize=2)
        ax1.plot(list(rer_param.keys()), list(rer_param.values()), marker='*', color='k', linestyle='', label="Model parameters")

        # Formatting
        ax1.set_ylabel(u'Relative Elongation Rate at 12°C (s$^{-1}$)')
        ax1.legend(prop={'size': 12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=3, mode="expand", borderaxespad=0.)
        ax1.legend(loc='upper left')
        ax1.set_xlabel('Phytomer rank')
        ax1.set_ylim(bottom=0., top=6e-6)
        ax1.set_xlim(left=2)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RER_comparison.PNG'), format='PNG', bbox_inches='tight', dpi=200)
        plt.close()

    # 4) Total C production vs. Root C allcoation
    df_roots = organs_postprocessing[organs_postprocessing['organ'] == 'roots'].copy()
    df_roots['day'] = df_roots['t'] // 24 + 1
    df_roots['Unloading_Sucrose_tot'] = df_roots['Unloading_Sucrose'] * df_roots['mstruct']
    Unloading_Sucrose_tot = df_roots.groupby(['day'])['Unloading_Sucrose_tot'].agg('sum')
    days = df_roots['day'].unique()

    axes_postprocessing['day'] = axes_postprocessing['t'] // 24 + 1
    Total_Photosynthesis = axes_postprocessing.groupby(['day'])['Tillers_Photosynthesis'].agg('sum')

    elements_postprocessing['day'] = elements_postprocessing['t'] // 24 + 1
    elements_postprocessing['sum_respi_tillers'] = elements_postprocessing['sum_respi'] * elements_postprocessing['nb_replications']
    Shoot_respiration = elements_postprocessing.groupby(['day'])['sum_respi_tillers'].agg('sum')
    Net_Photosynthesis = Total_Photosynthesis - Shoot_respiration

    share_net_roots_live = Unloading_Sucrose_tot.to_numpy() / Net_Photosynthesis.replace(0, np.nan).to_numpy() * 100

    fig, ax = plt.subplots()
    line1 = ax.plot(days, Net_Photosynthesis, label=u'Net_Photosynthesis')
    line2 = ax.plot(days, Unloading_Sucrose_tot, label=u'C unloading to roots')

    ax2 = ax.twinx()
    line3 = ax2.plot(days, share_net_roots_live, label=u'Net C Shoot production sent to roots (%)', color='red')

    lines = line1 + line2 + line3
    labs = [line.get_label() for line in lines]
    ax.legend(lines, labs, loc='center left', prop={'size': 10}, framealpha=0.5, bbox_to_anchor=(1, 0.815), borderaxespad=0.)

    ax.set_xlabel('Days')
    ax2.set_ylim([0, 200])
    ax.set_ylabel(u'C (µmol C.day$^{-1}$ )')
    ax2.set_ylabel(u'Ratio (%)')
    ax.set_title('C allocation to roots')
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'C_allocation.PNG'), dpi=200, format='PNG', bbox_inches='tight')

    # 5) C usages relative to Net Photosynthesis
    if not elements_postprocessing.empty:
        df_phloem = organs_postprocessing[organs_postprocessing['organ'] == 'phloem'].copy()
        df_phloem['day'] = df_phloem['t'] // 24 + 1

        AMINO_ACIDS_C_RATIO = 4.15  #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
        AMINO_ACIDS_N_RATIO = 1.25  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

        # Photosynthesis
        elements_postprocessing['Photosynthesis_tillers'] = elements_postprocessing['Photosynthesis'].fillna(0) * elements_postprocessing['nb_replications'].fillna(1.)
        Tillers_Photosynthesis_Ag = elements_postprocessing.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
        C_usages = pd.DataFrame({'t': Tillers_Photosynthesis_Ag['t']})
        C_usages['C_produced'] = np.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)

        # Respiration
        C_usages['Respi_roots'] = np.cumsum(axes_postprocessing.C_respired_roots)
        C_usages['Respi_shoot'] = np.cumsum(axes_postprocessing.C_respired_shoot)

        # Exudation
        C_usages['exudation'] = np.cumsum(axes_postprocessing.C_exuded.fillna(0))

        # Structural growth
        C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        C_usages['Structure_roots'] = np.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

        hz_postprocessing['C_consumption_mstruct'] = hz_postprocessing.sucrose_consumption_mstruct.fillna(0) + hz_postprocessing.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        hz_postprocessing['C_consumption_mstruct_tillers'] = hz_postprocessing['C_consumption_mstruct'] * hz_postprocessing['nb_replications']
        C_consumption_mstruct_shoot = hz_postprocessing.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
        C_usages['Structure_shoot'] = np.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True)).apply(float)

        # Non structural C
        df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS.iloc[0]
        C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

        elements_postprocessing['C_NS'] = elements_postprocessing.sucrose.fillna(0) + elements_postprocessing.fructan.fillna(0) + elements_postprocessing.starch.fillna(0) + (
                elements_postprocessing.amino_acids.fillna(0) + elements_postprocessing.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        elements_postprocessing['C_NS_tillers'] = elements_postprocessing['C_NS'] * elements_postprocessing['nb_replications'].fillna(1.)
        C_elt = elements_postprocessing.groupby(['t']).agg({'C_NS_tillers': 'sum'})

        hz_postprocessing['C_NS'] = hz_postprocessing.sucrose.fillna(0) + hz_postprocessing.fructan.fillna(0) + (hz_postprocessing.amino_acids.fillna(0) + hz_postprocessing.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        hz_postprocessing['C_NS_tillers'] = hz_postprocessing['C_NS'] * hz_postprocessing['nb_replications'].fillna(1.)
        C_hz = hz_postprocessing.groupby(['t']).agg({'C_NS_tillers': 'sum'})

        df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO

        C_NS_autre = df_roots.set_index('t').C_NS.add(C_elt.C_NS_tillers, fill_value=0).add(C_hz.C_NS_tillers, fill_value=0)
        C_NS_autre_init = C_NS_autre - C_NS_autre[0]
        C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

        # Total
        C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) / \
                               C_usages.C_produced.replace(0, np.nan).to_numpy()

        # ----- Graph
        fig, ax = plt.subplots()
        ax.plot(C_usages.t, C_usages.Structure_shoot / C_usages.C_produced * 100,
                label=u'Structural mass - Shoot', color='g')
        ax.plot(C_usages.t, C_usages.Structure_roots / C_usages.C_produced * 100,
                label=u'Structural mass - Roots', color='r')
        ax.plot(C_usages.t, (C_usages.NS_phloem + C_usages.NS_other) / C_usages.C_produced * 100, label=u'Non-structural C', color='darkorange')
        ax.plot(C_usages.t, (C_usages.Respi_roots + C_usages.Respi_shoot) / C_usages.C_produced.replace(0, np.nan).to_numpy() * 100, label=u'C loss by respiration', color='b')
        ax.plot(C_usages.t, C_usages.exudation / C_usages.C_produced.replace(0, np.nan).to_numpy() * 100, label=u'C loss by exudation', color='c')

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel('Time (h)')
        ax.set_ylabel(u'Carbon usages : Photosynthesis (%)')
        ax.set_ylim(bottom=0, top=100.)

        fig.suptitle(u'Total cumulated usages are ' + str(round(C_usages.C_budget.tail(1) * 100, 0)) + u' % of Photosynthesis')

        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'C_usages_cumulated.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

    # 6) RUE
    if not elements_postprocessing.empty:
        elements_postprocessing['PARa_MJ'] = elements_postprocessing['PARa'] * elements_postprocessing['green_area'] * elements_postprocessing['nb_replications'] * 3600 / 4.6 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
        elements_postprocessing['RGa_MJ'] = elements_postprocessing['PARa'] * elements_postprocessing['green_area'] * elements_postprocessing['nb_replications'] * 3600 / 2.02 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
        PARa = elements_postprocessing.groupby(['day'])['PARa_MJ'].agg('sum')
        PARa_cum = np.cumsum(PARa)
        days = elements_postprocessing['day'].unique()

        sum_dry_mass_shoot = axes_postprocessing.groupby(['day'])['sum_dry_mass_shoot'].agg('max')
        sum_dry_mass = axes_postprocessing.groupby(['day'])['sum_dry_mass'].agg('max')

        # RUE_shoot = np.polyfit(PARa_cum, sum_dry_mass_shoot, 1)[0]
        sum_dry_mass_shoot_numeric = pd.to_numeric(sum_dry_mass_shoot, errors='coerce')
        sum_dry_mass_plant_numeric = pd.to_numeric(sum_dry_mass, errors='coerce')
        para_cum_series = pd.Series(PARa_cum, index=PARa.index)
        df_temp_shoot = pd.DataFrame({'X': para_cum_series, 'Y': sum_dry_mass_shoot_numeric}).dropna()
        df_temp_plant = pd.DataFrame({'X': para_cum_series, 'Y': sum_dry_mass_plant_numeric}).dropna()
        RUE_shoot = np.polyfit(df_temp_shoot['X'], df_temp_shoot['Y'], 1)[0]
        RUE_plant = np.polyfit(df_temp_plant['X'], df_temp_plant['Y'], 1)[0]

        fig, ax = plt.subplots()
        ax.plot(PARa_cum, sum_dry_mass_shoot.dropna(), label='Shoot dry mass (g)')
        ax.plot(PARa_cum, sum_dry_mass.dropna(), label='Plant dry mass (g)')
        ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
        ax.set_xlabel('Cumulative absorbed PAR (MJ)')
        ax.set_ylabel('Dry mass (g)')
        ax.set_title('RUE')
        plt.text(max(PARa_cum) * 0.02, max(sum_dry_mass) * 0.95, 'RUE shoot : {0:.2f} , RUE plant : {1:.2f}'.format(round(RUE_shoot, 2), round(RUE_plant, 2)))
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE.PNG'), dpi=200, format='PNG', bbox_inches='tight')

        fig, ax = plt.subplots()
        ax.plot(days, sum_dry_mass_shoot.dropna(), label='Shoot dry mass (g)')
        ax.plot(days, sum_dry_mass.dropna(), label='Plant dry mass (g)')
        ax.plot(days, PARa_cum, label='Cumulative absorbed PAR (MJ)')
        ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
        ax.set_xlabel('Days')
        ax.set_title('RUE investigations')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE2.PNG'), dpi=200, format='PNG', bbox_inches='tight')

    # 7) Sum thermal time
    fig, ax = plt.subplots()
    ax.plot(df_MS['t'], df_MS['sum_TT'])
    ax.set_xlabel('Hours')
    ax.set_ylabel('Thermal Time')
    ax.set_title('Thermal Time')
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SumTT.PNG'), dpi=200, format='PNG', bbox_inches='tight')

    # 7) Residual N : ratio_N_mstruct_max
    df_elt_MS = elements_outputs.loc[elements_outputs.axis == 'MS']
    df_elt_MS = df_elt_MS.loc[df_elt_MS.mstruct != 0]
    df_elt_MS['N_content_total'] = df_elt_MS['N_content_total'] * 100
    x_name = 't'
    x_label = 'Time (Hour)'
    graph_variables_ph_elements = {'N_content_total': u'N content in green + senesced tissues (% mstruct)'}
    for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
        for variable_name, variable_label in graph_variables_ph_elements.items():

            graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
            plot_outputs(elements_outputs,
                      x_name=x_name,
                      y_name=variable_name,
                      x_label=x_label,
                      y_label=variable_label,
                      colors=[colors[i - 1] for i in elements_outputs.metamer.unique().tolist()],
                      filters={'organ': org_ph},
                      plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                      explicit_label=False)


def color_MTG_Nitrogen(g, df, t, SCREENSHOT_DIRPATH):
    
    def color_map(N):
        if 0 <= N <= 0.5:  # TODO: organe senescent (prendre prop)
            vid_colors = [150, 100, 0]
        elif 0.5 < N < 5:  # Fvertes
            vid_colors = [int(255 - N*51), int(255 - N * 20), 50]
        else:
            vid_colors = [0, 155, 0]
        return vid_colors

    def calculate_Total_Organic_Nitrogen(amino_acids, proteins, Nstruct):
        """Total amount of organic N (amino acids + proteins + Nstruct).

        :param float amino_acids: Amount of amino acids (µmol N)
        :param float proteins: Amount of proteins (µmol N)
        :param float Nstruct: Structural N mass (g)

        :return: Total amount of organic N (mg)
        :rtype: float
        """
        return (amino_acids + proteins) * 14E-3 + Nstruct * 1E3

    colors = {}

    groups_df = df.groupby(['plant', 'axis', 'metamer', 'organ', 'element'])
    for vid in g.components_at_scale(g.root, scale=5):
        pid = int(g.index(g.complex_at_scale(vid, scale=1)))
        axid = g.property('label')[g.complex_at_scale(vid, scale=2)]
        mid = int(g.index(g.complex_at_scale(vid, scale=3)))
        org = g.property('label')[g.complex_at_scale(vid, scale=4)]
        elid = g.property('label')[vid]
        id_map = (pid, axid, mid, org, elid)
        if id_map in groups_df.groups.keys():
            N = (g.property('proteins')[vid] * 14E-3) / groups_df.get_group(id_map)['mstruct'].iloc[0]
            # N = (calculate_Total_Organic_Nitrogen(g.property('amino_acids')[vid], g.property('proteins')[vid], g.property('Nstruct')[vid])) / g.property('mstruct')[vid]
            colors[vid] = color_map(N)
        else:
            g.property('geometry')[vid] = None

    # plantgl
    s = to_plantgl(g, colors=colors)[0]
    Viewer.add(s)
    Viewer.camera.setPosition(Vector3(83.883, 12.3239, 93.4706))
    Viewer.camera.lookAt(Vector3(0., 0, 50))
    Viewer.saveSnapshot(os.path.join(SCREENSHOT_DIRPATH, 'Day_{}.png'.format(t/24+1)))


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None,
                              precision=4, overwrite_desired_data=False):
    """Compare

            difference = actual_data_df - desired_data_df

       to

            tolerance = 10**-precision * (1 + abs(desired_data_df))

        where

            desired_data_df = pd.read_csv(os.path.join(data_dirpath, desired_data_filename))

        If difference > tolerance, then raise an AssertionError.

    :param str data_dirpath: The path of the directory where to find the data to compare.
    :param pandas.DataFrame actual_data_df: The computed data.
    :param str desired_data_filename: The file name of the expected data.
    :param str actual_data_filename: If not None, save the computed data to `actual_data_filename`, in directory `data_dirpath`. Default is None.
    :param int precision: The precision to use for the comparison. Default is `4`.
    :param bool overwrite_desired_data: If True the comparison between actual and desired data is not run. Instead, the desired data will be overwritten using actual data. To be used with caution.
    """

    relative_tolerance = 10 ** -precision
    absolute_tolerance = relative_tolerance

    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    if actual_data_filename is not None:
        # save actual outputs to CSV file
        actual_data_filepath = os.path.join(data_dirpath, actual_data_filename)
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(precision))

    if overwrite_desired_data:
        warnings.warn('!!! Unit test is running with overwrite_desired_data !!!')
        desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
        actual_data_df.to_csv(desired_data_filepath, na_rep='NA', index=False)

    else:
        # keep only numerical data (np.testing can compare only numerical data)
        for column in ('axis', 'organ', 'element', 'is_growing'):
            if column in desired_data_df.columns:
                del desired_data_df[column]
                del actual_data_df[column]

        # convert the actual outputs to floats
        actual_data_df = actual_data_df.astype(np.float)

        # compare actual data to desired data
        np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, relative_tolerance,
                                   absolute_tolerance)