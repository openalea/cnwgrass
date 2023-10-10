"""
    turgorgrowth.postprocessing
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.postprocessing` defines post-processing to apply
    on Turgor-Wheat outputs, and provides a front-end to automatize the generation of graphs
    for validation of the outputs.

    Please use front-ends :func:`postprocessing` and :func:`generate_graphs`.

    :license: CeCILL-C, see LICENSE for details.
"""

from __future__ import division  # use "//" to do integer division
import os

import pandas as pd

from turgorgrowth import simulation as turgorgrowth_simulation, parameters as turgorgrowth_parameters
from turgorgrowth import tools as turgorgrowth_tools

#: the time index
T_INDEX = turgorgrowth_simulation.Simulation.T_INDEX

#: the index to locate the plants in the modeled system
PLANTS_INDEXES = turgorgrowth_simulation.Simulation.PLANTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PLANTS_INDEXES`
PLANTS_T_INDEXES = turgorgrowth_simulation.Simulation.PLANTS_T_INDEXES
#: plants post-processing variables
PLANTS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PLANTS_T_INDEXES`, :attr:`PLANTS_RUN_VARIABLES <turgorgrowth.simulation.Simulation.PLANTS_RUN_VARIABLES>` and :attr:`PLANTS_POSTPROCESSING_VARIABLES`
PLANTS_RUN_POSTPROCESSING_VARIABLES = set(PLANTS_T_INDEXES + turgorgrowth_simulation.Simulation.PLANTS_RUN_VARIABLES + PLANTS_POSTPROCESSING_VARIABLES)

#: the indexes to locate the soil in the modeled system
SOIL_INDEXES = turgorgrowth_simulation.Simulation.SOIL_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`SOIL_INDEXES`
SOIL_T_INDEXES = turgorgrowth_simulation.Simulation.SOIL_T_INDEXES
#: soil post-processing variables
SOIL_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`SOIL_T_INDEXES`, :attr:`SOIL_RUN_VARIABLES <turgorgrowth.simulation.Simulation.SOIL_RUN_VARIABLES>` and :attr:`SOIL_POSTPROCESSING_VARIABLES`
SOIL_RUN_POSTPROCESSING_VARIABLES = SOIL_T_INDEXES + turgorgrowth_simulation.Simulation.SOIL_RUN_VARIABLES + SOIL_POSTPROCESSING_VARIABLES

#: the indexes to locate the axes in the modeled system
AXES_INDEXES = turgorgrowth_simulation.Simulation.AXES_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`AXES_INDEXES`
AXES_T_INDEXES = turgorgrowth_simulation.Simulation.AXES_T_INDEXES
#: axes post-processing variables
AXES_POSTPROCESSING_VARIABLES = []

#: concatenation of :attr:`AXES_T_INDEXES`, :attr:`AXES_RUN_VARIABLES <cnwheat.simulation.Simulation.AXES_RUN_VARIABLES>` and :attr:`AXES_POSTPROCESSING_VARIABLES`
AXES_RUN_POSTPROCESSING_VARIABLES = set(AXES_T_INDEXES + turgorgrowth_simulation.Simulation.AXES_RUN_VARIABLES + AXES_POSTPROCESSING_VARIABLES)

#: the indexes to locate the phytomers in the modeled system
PHYTOMERS_INDEXES = turgorgrowth_simulation.Simulation.PHYTOMERS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
PHYTOMERS_T_INDEXES = turgorgrowth_simulation.Simulation.PHYTOMERS_T_INDEXES
#: phytomers post-processing variables
PHYTOMERS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PHYTOMERS_T_INDEXES`, :attr:`PHYTOMERS_RUN_VARIABLES <cnwheat.simulation.Simulation.PHYTOMERS_RUN_VARIABLES>` and :attr:`PHYTOMERS_POSTPROCESSING_VARIABLES`
PHYTOMERS_RUN_POSTPROCESSING_VARIABLES = set(PHYTOMERS_T_INDEXES + turgorgrowth_simulation.Simulation.PHYTOMERS_RUN_VARIABLES + PHYTOMERS_POSTPROCESSING_VARIABLES)

#: the indexes to locate the organs in the modeled system
ORGANS_INDEXES = turgorgrowth_simulation.Simulation.ORGANS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ORGANS_INDEXES`
ORGANS_T_INDEXES = turgorgrowth_simulation.Simulation.ORGANS_T_INDEXES
#: organs post-processing variables
ORGANS_POSTPROCESSING_VARIABLES = []
ORGANS_RUN_VARIABLES_ADDITIONAL = []
#: concatenation of :attr:`ORGANS_T_INDEXES`, :attr:`ORGANS_RUN_VARIABLES <cnwheat.simulation.Simulation.ORGANS_RUN_VARIABLES>` and :attr:`ORGANS_POSTPROCESSING_VARIABLES`
ORGANS_RUN_POSTPROCESSING_VARIABLES = set(ORGANS_T_INDEXES + turgorgrowth_simulation.Simulation.ORGANS_RUN_VARIABLES + ORGANS_POSTPROCESSING_VARIABLES + ORGANS_RUN_VARIABLES_ADDITIONAL)

#: the indexes to locate the hidden zones in the modeled system
HIDDENZONE_INDEXES = turgorgrowth_simulation.Simulation.HIDDENZONE_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
HIDDENZONE_T_INDEXES = turgorgrowth_simulation.Simulation.HIDDENZONE_T_INDEXES
#: hidden zones post-processing variables
HIDDENZONE_POSTPROCESSING_VARIABLES = []
HIDDENZONE_RUN_VARIABLES_ADDITIONAL = []
#: concatenation of :attr:`HIDDENZONE_T_INDEXES`, :attr:`HIDDENZONE_RUN_VARIABLES <turgorgrowth.simulation.Simulation.HIDDENZONE_RUN_VARIABLES>` and :attr:`HIDDENZONE_POSTPROCESSING_VARIABLES`
HIDDENZONE_RUN_POSTPROCESSING_VARIABLES = HIDDENZONE_T_INDEXES + turgorgrowth_simulation.Simulation.HIDDENZONE_RUN_VARIABLES + HIDDENZONE_RUN_VARIABLES_ADDITIONAL + HIDDENZONE_POSTPROCESSING_VARIABLES

#: the indexes to locate the elements in the modeled system
ELEMENTS_INDEXES = turgorgrowth_simulation.Simulation.ELEMENTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
ELEMENTS_T_INDEXES = turgorgrowth_simulation.Simulation.ELEMENTS_T_INDEXES
#: elements post-processing variables
ELEMENTS_POSTPROCESSING_VARIABLES = []
ELEMENTS_RUN_VARIABLES_ADDITIONAL = []
#: concatenation of :attr:`ELEMENTS_T_INDEXES`, :attr:`ELEMENTS_RUN_VARIABLES <turgorgrowth.simulation.Simulation.ELEMENTS_RUN_VARIABLES>` and :attr:`ELEMENTS_POSTPROCESSING_VARIABLES`
ELEMENTS_RUN_POSTPROCESSING_VARIABLES = ELEMENTS_T_INDEXES + turgorgrowth_simulation.Simulation.ELEMENTS_RUN_VARIABLES + ELEMENTS_RUN_VARIABLES_ADDITIONAL + ELEMENTS_POSTPROCESSING_VARIABLES


# ---------------------------------------------------
#           POST-PROCESSING FUNCTIONS               -
#           DO NOT USE THEM DIRECTLY                -
# ---------------------------------------------------

class Roots:
    """
    Post-processing to apply on Xylem outputs.
    """
    pass

class Xylem:
    """
    Post-processing to apply on Xylem outputs.
    """
    pass

class HiddenZone:
    """
    Post-processing to apply on HiddenZone outputs.
    """
    pass


class Element:
    """
    Post-processing to apply on Element outputs.
    """
    pass

class Organ:
    """
    Post-processing to apply on Organ outputs.
    """
    pass

# ------------------------------------------------------------------------------------------------
#                                   POST-PROCESSING FRONT-END                                    -
#           PLEASE USE THIS FUNCTION TO APPLY POST-PROCESSING ON THE OUTPUT OF TURGOR-GROWTH     -
# ------------------------------------------------------------------------------------------------

def postprocessing(plants_df=None, axes_df=None, metamers_df=None, hiddenzones_df=None, organs_df=None, elements_df=None, soil_df=None, delta_t=1):
    """
    Compute post-processing from Turgor-Growth outputs, and format the post-processing to :class:`dataframes <pandas.DataFrame>`.

    For each post-processing output dataframe:

        * compute post-processing from Turgor-Growth outputs,
        * concatenate Turgor-Growth outputs and post-processing and place the results in a jointed dataframe,
        * reorder the columns of the dataframes according to
          :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`, :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`,  :attr:`ORGAN_RUN_POSTPROCESSING_VARIABLES`
        * and convert the indexes of plants and metamers to integers (if relevant).

    :param pandas.DataFrame hiddenzones_df: Turgor-Growth outputs at hidden zone scale (see :attr:`simulation.Simulation.HIDDENZONE_RUN_VARIABLES`)
    :param pandas.DataFrame elements_df: Turgor-Growth outputs at element scale (see :attr:`simulation.Simulation.ELEMENTS_RUN_VARIABLES`)
    :param pandas.DataFrame organ_df: Turgor-Growth outputs at xylem scale (see :attr:`simulation.Simulation. ORGAN_RUN_VARIABLES`)
    :param pandas.DataFrame soil_df: Turgor-Growth outputs at soil scale (see :attr:`simulation.Simulation.SOIL_RUN_VARIABLES`)
    :param float delta_t: Delta t between 2 outputs (in seconds).

    :return: :class:`dataframes <pandas.DataFrame>` of post-processing for each scale:
            * hidden zone (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
            * element (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
            * xylem (see :attr:`XYLEM_RUN_POSTPROCESSING_VARIABLES`)
            * and soil (see :attr:`SOIL_RUN_POSTPROCESSING_VARIABLES`)


    :rtype tuple [pandas.DataFrame]
    """

    returned_dataframes = []

    # # plants
    # if plants_df is not None:
    #     pp_plants_df = pd.concat([plants_df, pd.DataFrame(columns=PLANTS_POSTPROCESSING_VARIABLES)], sort=False)
    #     pp_plants_df = pp_plants_df.reindex(PLANTS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
    #     pp_plants_df['plant'] = pp_plants_df['plant'].astype(int)
    #     returned_dataframes.append(pp_plants_df)
    # else:
    #     returned_dataframes.append(pd.DataFrame({'A': []}))

    # axes
    if axes_df is not None:
        axes_df = axes_df[axes_df['axis'] == 'MS'].copy()  # TODO : Temporary !
        pp_axes_df = pd.concat([axes_df, pd.DataFrame(columns=AXES_POSTPROCESSING_VARIABLES)], sort=False)
        pp_axes_df = pp_axes_df.reindex(AXES_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_axes_df['plant'] = pp_axes_df['plant'].astype(int)
        returned_dataframes.append(pp_axes_df)
        # xylem
        # roots

    # # metamers
    # if metamers_df is not None:
    #     pp_metamers_df = pd.concat([metamers_df, pd.DataFrame(columns=PHYTOMERS_POSTPROCESSING_VARIABLES)], sort=False)
    #     pp_metamers_df = pp_metamers_df.reindex(PHYTOMERS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
    #     pp_metamers_df[['plant', 'metamer']] = pp_metamers_df[['plant', 'metamer']].astype(int)
    #     returned_dataframes.append(pp_metamers_df)
    # else:
    #     returned_dataframes.append(pd.DataFrame({'A': []}))

    # hidden zones
    if hiddenzones_df is not None:
        pp_hiddenzones_df = pd.concat([hiddenzones_df, pd.DataFrame(columns=HIDDENZONE_POSTPROCESSING_VARIABLES)], sort=True)
        pp_hiddenzones_df = pp_hiddenzones_df.reindex(columns=HIDDENZONE_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_hiddenzones_df[['plant', 'metamer']] = pp_hiddenzones_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_hiddenzones_df)

    # elements
    if elements_df is not None:
        pp_elements_df = pd.concat([elements_df, pd.DataFrame(columns=ELEMENTS_POSTPROCESSING_VARIABLES)], sort=True)
        pp_elements_df = pp_elements_df.reindex(columns=ELEMENTS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        ##pp_elements_df[['plant', 'metamer']] = pp_elements_df[['plant', 'metamer']].astype(int)
        ##returned_dataframes.append(pp_elements_df)

        grouped = elements_df.groupby('organ')
        for organ_type, parameters_class in (('blade', turgorgrowth_parameters.LAMINA_ELEMENT_PARAMETERS), ('internode', turgorgrowth_parameters.INTERNODE_ELEMENT_PARAMETERS), ('sheath', turgorgrowth_parameters.SHEATH_ELEMENT_PARAMETERS)):
            if organ_type not in grouped.groups:
                continue
            group = grouped.get_group(organ_type)
            if len(group) == 0:  # TODO: faire mm tri que dans simulation de cnwheat (surface nulle)
                continue
            curr_organ_elements_df = elements_df.loc[group.index]
            pp_curr_organ_elements_df = pp_elements_df.loc[group.index]
        pp_elements_df = pp_elements_df.reindex(columns=ELEMENTS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_elements_df[['plant', 'metamer']] = pp_elements_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_elements_df)

    # organs
    if organs_df is not None and axes_df is not None:
        axes_df = axes_df[axes_df['axis'] == 'MS'].copy()  # TODO : Temporary !
        pp_organs_df = pd.concat([organs_df, pd.DataFrame(columns=ORGANS_POSTPROCESSING_VARIABLES)], sort=False)
        # roots
        ##roots_df = organs_df.loc[organs_df.organ == 'roots']
        # xylem
        xylem_df = organs_df.loc[organs_df.organ == 'xylem']
        pp_organs_df = pp_organs_df.reindex(columns=ORGANS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_organs_df['plant'] = pp_organs_df['plant'].astype(int)
        returned_dataframes.append(pp_organs_df)

    # soil
    if soil_df is not None:
        pp_soil_df = pd.concat([soil_df, pd.DataFrame(columns=SOIL_POSTPROCESSING_VARIABLES)])
        pp_soil_df = pp_soil_df.reindex(columns=SOIL_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_soil_df[['plant']] = pp_soil_df[['plant']].astype(int)
        returned_dataframes.append(pp_soil_df)

    return tuple(returned_dataframes)


# ---------------------------------------------------------------------
#                       GRAPHS GENERATION FRONT-END                   -
#           PLEASE USE THIS FUNCTION FOR THE GENERATION OF GRAPHS     -
# ---------------------------------------------------------------------

def generate_graphs(axes_df=None, hiddenzones_df=None, organs_df=None, elements_df=None, soil_df=None, graphs_dirpath='.'):

    """
    Generate graphs to validate the outputs of Turgor-Growth, and save them in directory `graphs_dirpath`.

    :param pandas.DataFrame hiddenzones_df: Turgor-Growth outputs at hidden zone scale (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame elements_df: Turgor-Growth outputs at element scale (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame organ_df: Turgor-Growth outputs at organ scale (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame soil_df: Turgor-Growth outputs at soil scale (see :attr:`SOIL_RUN_POSTPROCESSING_VARIABLES`)
    :param str graphs_dirpath: the path of the directory to save the generated graphs
    """

    x_name = 't'
    x_label = 'Time (Hour)'
    colors = ['blue', 'darkorange', 'green', 'red', 'darkviolet', 'gold', 'magenta', 'brown', 'darkcyan', 'grey', 'lime']
    colors = colors + colors + colors + colors + colors

    # 1) Photosynthetic organs
    if elements_df is not None:
        graph_variables_ph_elements = {'length': u'Length (m)', 'water_flux_from_hz' : u'Water flux from hz to element (g H2O)',
                                       'osmotic_water_potential': u'Osmotic water potential (MPa)', 'thickness': u'Thickness (m)', 'total_water_potential': u'Total water potential (MPa)',
                                       'turgor_water_potential': u'Turgor water potential (MPa)', 'water_content': u'Water content (g)', 'water_influx': u'Water flow from Xylem (g)',  'width': u'Width (m)',
                                       'resistance': u'Resistance (MPa s g$^{-1}$)', 'volume': u'Volume m3)', 'sucrose': u'Sucrose', 'proteins': u'Proteins', 'amino_acids': u'Amino acids'}
    
        for org_ph in (['blade'], ['sheath'], ['internode']):
            for variable_name, variable_label in graph_variables_ph_elements.items():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                turgorgrowth_tools.plot_turgorgrowth_ouputs(elements_df,
                                                            x_name=x_name,
                                                            y_name=variable_name,
                                                            x_label=x_label,
                                                            y_label=variable_label,
                                                            colors=[colors[i - 1] for i in elements_df.metamer.unique().tolist()],
                                                            filters={'organ': org_ph},
                                                            plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                            explicit_label=False)

    # 2) Hidden zones
    if hiddenzones_df is not None:
        graph_variables_hiddenzones = { 'leaf_pseudo_age': u'Leaf pseudo age (°Cd)', 'phi_length': u'Extensibility parameter  for length (Mpa-1)', 'phi_width': u'Extensibility parameter  for width (Mpa-1)', 'phi_thickness': u'Extensibility parameter  for thickness (Mpa-1)',
                                        'length_leaf_emerged': u'length of emerged part of the growing leaf (m)',
                                        'leaf_L': 'Total leaf length (m)', 'length': u'Length of hz (m)',
                                       'osmotic_water_potential': u'Osmotic water potential (MPa)', 'width': u'Width (m)', 'total_water_potential': u'Total water potential (MPa)',
                                       'turgor_water_potential': u'Turgor water potential (MPa)', 'water_content': u'Water content (g)', 'water_influx': u'Water flow Xylem (g)', 'water_outflow': u'Water flow HZ (g)',
                                       'resistance': u'Resistance (MPa s g$^{-1}$)', 'thickness': u'Thickness (m)', 'volume': u'Volume m3)', 'sucrose': u'Sucrose', 'proteins': u'Proteins', 'amino_acids': u'Amino acids'}
    
        for variable_name, variable_label in graph_variables_hiddenzones.items():
            graph_name = variable_name + '_hz' + '.PNG'
            turgorgrowth_tools.plot_turgorgrowth_ouputs(hiddenzones_df,
                                                        x_name=x_name,
                                                        y_name=variable_name,
                                                        x_label=x_label,
                                                        y_label=variable_label,
                                                        colors=[colors[i - 1] for i in hiddenzones_df.metamer.unique().tolist()],
                                                        filters={'plant': 1, 'axis': 'MS'},
                                                        plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                        explicit_label=False)

    # 3) Roots and xylem
    if organs_df is not None:
        graph_variables_organs = {'total_water_potential': u'Total water potential (MPa)',
                                  'soil_water_potential': u'Soil water potential (Mpa)', 'SRWC': u'SRWC (%)'}

        for org in (['roots'], ['xylem']):
            for variable_name, variable_label in graph_variables_organs.items():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                turgorgrowth_tools.plot_turgorgrowth_ouputs(organs_df,
                                                        x_name=x_name,
                                                        y_name=variable_name,
                                                        x_label=x_label,
                                                        y_label=variable_label,
                                                        colors=['blue'],
                                                        filters={'organ': org},
                                                        plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                        explicit_label=False)

    # 4) Axes
    graph_variables_axes = {'Total_Transpiration': u'Total transpiration (g H2O)', 'Growth':u'Growth (g H2O)', 'total_water_influx':u' Water flux from xylem to HZ and photosynthetic organs (g H2O)'}

    for variable_name, variable_label in graph_variables_axes.items():
        graph_name = variable_name + '_axis' + '.PNG'
        turgorgrowth_tools.plot_turgorgrowth_ouputs(axes_df,
                                                    x_name=x_name,
                                                    y_name=variable_name,
                                                    x_label=x_label,
                                                    y_label=variable_label,
                                                    colors=['blue'],
                                                    filters={'plant': 1, 'axis': 'MS'},
                                                    plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                    explicit_label=False)

    # 5) Soil
    if soil_df is not None:
        graph_variables_soil = {'soil_water_potential': u'Total water potential of soil (MPa)'}
        for variable_name, variable_label in graph_variables_soil.items():
            graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
            turgorgrowth_tools.plot_turgorgrowth_ouputs(soil_df,
                                                            x_name=x_name,
                                                            y_name=variable_name,
                                                            x_label=x_label,
                                                            y_label=variable_label,
                                                            colors=['blue'],
                                                            filters={'plant': 1, 'axis': 'MS'},
                                                            plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                            explicit_label=False)
