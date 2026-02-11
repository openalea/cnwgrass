# -*- coding: latin-1 -*-

import os
import sys
from itertools import cycle
import warnings
import logging
import logging.config
import json

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

"""
    cnwheat.tools
    ~~~~~~~~~~~~~

    This module provides tools to help for the validation of the outputs: 
    
        * plot of multiple variables on the same graph, 
        * set up of loggers,
        * quantitative comparison test,
        * and progress-bar to follow the evolution of long simulations.  

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
    .. seealso:: Barillot et al. 2016.
"""

OUTPUTS_INDEXES = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']  #: All the possible indexes of CN-Wheat outputs


class DataWarning(UserWarning):
    """Raised when there is no data to plot for a variable."""
    def __init__(self, variable, keys):
        self.message = 'No data to plot for variable {} at {}.'.format(variable, keys)

    def __str__(self):
        return repr(self.message)


# show all DataWarning (not only the first one which occurred)
warnings.simplefilter('always', DataWarning)


def plot_cnwheat_ouputs(outputs, x_name, y_name, x_label='', y_label='', x_lim=None, title=None, filters={}, plot_filepath=None, colors=[], linestyles=[], explicit_label=True, kwargs={}):
    """Plot `outputs`, with x=`x_name` and y=`y_name`.

    The general algorithm is:

        * find the scale of `outputs` and keep only the needed columns,
        * apply `filters` to `outputs` and make groups according to the scale,
        * plot each group as a new line,
        * save or display the plot.

    :param pandas.DataFrame outputs: The outputs of CN-Wheat.
    :param str x_name: x axis of the plot.
    :param str y_name: y axis of the plot.
    :param str x_label: The x label of the plot. Default is ''.
    :param str or unicode y_label: The y label of the plot. Default is ''.
    :param float x_lim: the x-axis limit.
    :param str title: the title of the plot. If None (default), create a title which is the concatenation of `y_name` and each scales which cardinality is one.
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

    :Examples:

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # in this example, 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose'.
    >>> plot(cnwheat_output_df, x_name = 't', y_name = 'Conc_Sucrose', x_label='Time (Hour)', y_label=u'[Sucrose] (µmol g$^{-1}$ mstruct)', title='{} = f({})'.format('Conc_Sucrose', 't'), filters={'plant': 1, 'axis': 'MS', 'organ': 'Lamina', 'element': 1})

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
            # construct a label with only the essential keys of the group name ; the essential keys are those for which cardinality is non zero
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

    ax.set_ylim(bottom=0.)

    if x_lim is not None:
        ax.set_xlim(left=0, right=x_lim)
    else:
        ax.set_xlim(left=0)

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


def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_model=False, log_compartments=False, log_derivatives=False, 
                  remove_old_logs=False):
    """Setup logging configuration.

    :param str config_filepath: The file path of the logging configuration.
    :param int level: The global level of the logging. Use either
          `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or
          `logging.CRITICAL`.
    :param bool log_model: if `True`, log the messages from :mod:`cnwheat.model`. `False` otherwise.
    :param bool log_compartments: if `True`, log the values of the compartments. `False` otherwise.
    :param bool log_derivatives: if `True`, log the values of the derivatives. `False` otherwise.
    :param bool remove_old_logs: if `True`, remove all files in the logs directory documented in `config_filepath`.
    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        if remove_old_logs:
            logs_dir = os.path.dirname(os.path.abspath(config['handlers']['file_info']['filename']))
            for logs_file in os.listdir(logs_dir):
                os.remove(os.path.join(logs_dir, logs_file))
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    cnwheat_model_logger = logging.getLogger('cnwheat.model')
    cnwheat_model_logger.disabled = not log_model  # set to False to log messages from cnwheat.model
    logging.getLogger('cnwheat.compartments').disabled = not log_compartments  # set to False to log the compartments
    logging.getLogger('cnwheat.derivatives').disabled = not log_derivatives  # set to False to log the derivatives


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None, precision=4, overwrite_desired_data=False):
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
    
    relative_tolerance = 10**-precision
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
        np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, relative_tolerance, absolute_tolerance)


class ProgressBarError(Exception): pass


class ProgressBar(object):
    """
    Display a console progress bar.
    """

    def __init__(self, bar_length=20, title='', block_character='#', uncomplete_character='-'):
        if bar_length <= 0:
            raise ProgressBarError('bar_length <= 0')
        self.bar_length = bar_length  #: the number of blocks in the progress bar. MUST BE GREATER THAN ZERO !
        self.t_max = 1  #: the maximum t that the progress bar can display. MUST BE GREATER THAN ZERO !
        self.block_interval = 1  #: the time interval of each block. MUST BE GREATER THAN ZERO !
        self.last_upper_t = 0  #: the last upper t displayed by the progress bar
        self.progress_mapping = {}  #: a mapping to optimize the refresh rate
        self.title = title  #: the title to write on the left side of the progress bar
        self.block_character = block_character  #: the character to represent a block
        self.uncomplete_character = uncomplete_character  #: the character to represent the uncompleted part of the progress bar

    def set_t_max(self, t_max):
        """"Set :attr:`t_max` and update other attributes accordingly.
        """
        if t_max <= 0:
            raise ProgressBarError('t_max <= 0')
        self.t_max = t_max
        self.block_interval = self.t_max / self.bar_length
        self.last_upper_t = 0
        self.progress_mapping.clear()

    def update(self, t):
        """Update the progress bar if needed.
        """
        t = min(t, self.t_max)
        if t < self.last_upper_t:
            return
        else:
            self.last_upper_t = t
        t_inf = t // self.block_interval * self.block_interval
        if t_inf not in self.progress_mapping:
            progress = t / self.t_max
            block = int(round(self.bar_length * progress))
            text = "\r{0}: [{1}] {2:>5d}% ".format(self.title, self.block_character * block + self.uncomplete_character * (self.bar_length - block), int(progress*100))
            self.progress_mapping[t_inf] = text
            sys.stdout.write(self.progress_mapping[t_inf])
            sys.stdout.flush()
