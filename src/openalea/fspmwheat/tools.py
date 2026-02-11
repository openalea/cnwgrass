# -*- coding: latin-1 -*-

import os

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from alinea.adel.mtg import to_plantgl
from openalea.plantgl.all import Viewer, Vector3

"""
    fspmwheat.tools
    ~~~~~~~~~~~~~~~

    This module provides convenient tools needed by the facades.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""


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

    .. note:: `shared_dataframe_to_update` is updated in-place. Thus, `shared_dataframe_to_update` keeps the same object's memory address.

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
        if np.issubdtype(np.int64, data_type) and new_shared_dataframe[column_name].isnull().values.any():  # Used to keep bool values
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
    
    
def plot_linear_regression(x_array, y_array, x_label='x', y_label='y', plot_filepath=None):
    """Perform a linear regression of `x_array` vs `y_array`
    and create a plot showing the fit against the original data.
    If `plot_filepath` is not None, save the plot to a PNG file. Otherwise display the plot.

    This is derived from http://learningzone.rspsoc.org.uk/index.php/Learning-Materials/Python-Scripting/6.4-Fitting-linear-equations,
    which is under license CC BY-NC-SA 3.0 (https://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US).

    :param numpy.ndarray x_array: The first set of measurements.
    :param numpy.ndarray y_array: The second set of measurements.
    :param str x_label: The label of the abscissa axis. Default is 'x'.
    :param str y_label: The label of the ordinates axis. Default is 'y'.
    :param str plot_filepath: The file path to save the plot in. If `None`, do not save the plot.

    :Examples:

    >>> import pandas as pd
    >>> modelmaker_output_df = pd.read_csv('modelmaker_output.csv') # 'modelmaker_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # 'cnwheat_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> plot_linear_regression(modelmaker_output_df.Sucrose_Phloem,
                               cnwheat_output_df.Sucrose_Phloem,
                               x_label='modelmaker_{}'.format('Sucrose_Phloem'),
                               y_label='cnwheat_{}'.format('Sucrose_Phloem'),
                               plot_filepath='compare.png')

    """
    # Perform fit
    (aCoeff, bCoeff, rVal, _, _) = stats.linregress(x_array, y_array)

    # Use fits to predict y output for a range of diameters
    x_samples_array = np.linspace(min(x_array), max(x_array), 1000)
    y_predict_array = aCoeff * x_samples_array + bCoeff

    # Create a string, showing the form of the equation (with fitted coefficients) and r squared value.
    # Coefficients are rounded to two decimal places.
    equation = 'y = {} x + {} (R$^2$ = {})'.format(round(aCoeff, 2), round(bCoeff, 2), round(rVal**2, 2))

    plt.figure()

    # Plot fit against original data
    plt.plot(x_array, y_array, '.')
    plt.plot(x_samples_array, y_predict_array)
    plt.title('{} vs {}'.format(x_label, y_label))

    x_label = 'x = {}'.format(x_label)
    plt.xlabel(x_label)
    y_label = 'y = {}'.format(y_label)
    plt.ylabel(y_label)

    plt.legend(['x vs y', equation])

    # Save plot
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()
        

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
