# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os

import statsmodels.api as sm

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.image as mpimg

from openalea.elongwheat import parameters as elongwheat_parameters
from example.example_fspmwheat.Scenarios_monoculms import tools

# ----- Import scenarios list and description
scenarios_df = pd.read_csv(os.path.join('inputs', 'scenarios_list.csv'), index_col='Scenario')
scenarios_df['Scenario'] = scenarios_df.index
if 'Scenario_label' not in scenarios_df.keys():
    scenarios_df['Scenario_label'] = ''
else:
    scenarios_df['Scenario_label'] = scenarios_df['Scenario_label'].fillna('')
scenarios = scenarios_df.Scenario


# ------- Def functions to make graphs
def graph_RER(scenario, scenario_graphs_dirpath=None, scenario_outputs_dirpath=None):
    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_outputs_dirpath:
        scenario_outputs_dirpath = os.path.join('outputs', scenario_name, 'outputs')

    # ------ RERmax from parameters : from parameters file or from scenarios_list.csv
    scenario_dict = tools.buildDic(scenarios_df.to_dict('index')[scenario])

    rer_param = {}
    if scenario_dict.get('elongwheat'):
        rer_param = scenario_dict['elongwheat'].get('RERmax', {})
    if rer_param == {}:
        rer_param = dict((k, v) for k, v in elongwheat_parameters.RERmax.items())

    # ------ Simulated RER

    # import simulation outputs
    data_RER = pd.read_csv(os.path.join(scenario_outputs_dirpath, 'hiddenzones_outputs.csv'))
    data_RER = data_RER[(data_RER.axis == 'MS') & (data_RER.metamer >= 4)].copy()
    data_RER.sort_values(['t', 'metamer'], inplace=True)
    data_teq = pd.read_csv(os.path.join(scenario_outputs_dirpath, 'axes_outputs.csv'))
    data_teq = data_teq[data_teq.axis == 'MS'].copy()

    # Time previous leaf emergence
    tmp = data_RER[data_RER.leaf_is_emerged]
    leaf_em = tmp.groupby('metamer', as_index=False)['t'].min()
    leaf_em['t_em'] = leaf_em.t
    prev_leaf_em = leaf_em
    prev_leaf_em.metamer = leaf_em.metamer + 1

    data_RER2 = pd.merge(data_RER, prev_leaf_em[['metamer', 't_em']], on='metamer')
    data_RER2 = data_RER2[data_RER2.t <= data_RER2.t_em]

    # SumTimeEq
    data_teq['SumTimeEq'] = np.cumsum(data_teq.delta_teq)
    data_RER3 = pd.merge(data_RER2, data_teq[['t', 'SumTimeEq']], on='t')

    # logL
    data_RER3['logL'] = np.log(data_RER3.leaf_L)

    # Estimate RER
    RER_sim = {}
    for leaf in data_RER3.metamer.drop_duplicates():
        Y = data_RER3.logL[data_RER3.metamer == leaf]
        X = data_RER3.SumTimeEq[data_RER3.metamer == leaf]
        X = sm.add_constant(X)
        mod = sm.OLS(Y, X)
        fit_RER = mod.fit()
        RER_sim[int(leaf)] = fit_RER.params['SumTimeEq']

    # ------ Graph
    fig, (ax1) = plt.subplots(1, figsize=(4, 3))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    x, y = zip(*sorted(RER_sim.items()))
    ax1.plot(x, y, label=r'Simulated RER', linestyle='-', color='g')

    # obs
    # ax1.errorbar(rer_obs.leaf, rer_obs.RER, yerr=rer_obs.RER_confint, marker='o', color='g', linestyle='', label="Observed RER", markersize=2)

    # Parameters
    ax1.plot( list(rer_param.keys()), list(rer_param.values()), marker='*', color='k', linestyle='', label="Model parameters")

    # Formatting
    ax1.set_ylabel(u'Relative Elongation Rate at 12°C (s$^{-1}$)')
    ax1.legend(prop={'size': 12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=3, mode="expand", borderaxespad=0.)
    ax1.legend(loc='upper left')
    ax1.set_xlabel('Phytomer rank')
    ax1.set_xlim(left=4)
    ax1.set_ylim(bottom=0., top=6e-6)
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'RER_comparison.PNG'), format='PNG', bbox_inches='tight', dpi=600)
    plt.close()


def graph_C_usages(scenario, scenario_graphs_dirpath=None, scenario_postprocessing_dirpath=None):
    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_postprocessing_dirpath:
        scenario_postprocessing_dirpath = os.path.join('outputs', scenario_name, 'postprocessing')

    # --- Import simulations prostprocessings
    df_elt = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'elements_postprocessing.csv'))
    df_org = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'organs_postprocessing.csv'))
    df_axe = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'axes_postprocessing.csv'))
    df_hz = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'hiddenzones_postprocessing.csv'))
    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_phloem = df_org[df_org['organ'] == 'phloem'].copy()

    # --- C usages relatif to Net Photosynthesis
    AMINO_ACIDS_C_RATIO = 4.15  #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.25  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

    # Photosynthesis
    df_elt['Photosynthesis_tillers'] = df_elt['Photosynthesis'].fillna(0) * df_elt['nb_replications'].fillna(1.)
    Tillers_Photosynthesis_Ag = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
    C_usages = pd.DataFrame({'t': Tillers_Photosynthesis_Ag['t']})
    C_usages['C_produced'] = np.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)

    # Respiration
    C_usages['Respi_roots'] = np.cumsum(df_axe.C_respired_roots)
    C_usages['Respi_shoot'] = np.cumsum(df_axe.C_respired_shoot)

    # Exudation
    C_usages['exudation'] = np.cumsum(df_axe.C_exudated.fillna(0))

    # Structural growth
    C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    C_usages['Structure_roots'] = np.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

    df_hz['C_consumption_mstruct'] = df_hz.sucrose_consumption_mstruct.fillna(0) + df_hz.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_hz['C_consumption_mstruct_tillers'] = df_hz['C_consumption_mstruct'] * df_hz['nb_replications']
    C_consumption_mstruct_shoot = df_hz.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
    C_usages['Structure_shoot'] = np.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True))

    # Non structural C
    df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS[0]
    C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

    df_elt['C_NS'] = df_elt.sucrose.fillna(0) + df_elt.fructan.fillna(0) + df_elt.starch.fillna(0) + (
                df_elt.amino_acids.fillna(0) + df_elt.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_elt['C_NS_tillers'] = df_elt['C_NS'] * df_elt['nb_replications'].fillna(1.)
    C_elt = df_elt.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_hz['C_NS'] = df_hz.sucrose.fillna(0) + df_hz.fructan.fillna(0) + (df_hz.amino_acids.fillna(0) + df_hz.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_hz['C_NS_tillers'] = df_hz['C_NS'] * df_hz['nb_replications'].fillna(1.)
    C_hz = df_hz.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO

    C_NS_autre = df_roots.C_NS.reset_index(drop=True) + C_elt.C_NS_tillers + C_hz.C_NS_tillers
    C_NS_autre_init = C_NS_autre - C_NS_autre[0]
    C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

    # Total
    C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) / \
                           C_usages.C_produced

    # ----- Graph
    fig, ax = plt.subplots()
    ax.plot(C_usages.t, C_usages.Structure_shoot / C_usages.C_produced * 100,
            label=u'Structural mass - Shoot', color='g')
    ax.plot(C_usages.t, C_usages.Structure_roots / C_usages.C_produced * 100,
            label=u'Structural mass - Roots', color='r')
    ax.plot(C_usages.t, (C_usages.NS_phloem + C_usages.NS_other) / C_usages.C_produced * 100, label=u'Non-structural C', color='darkorange')
    ax.plot(C_usages.t, (C_usages.Respi_roots + C_usages.Respi_shoot) / C_usages.C_produced * 100, label=u'C loss by respiration', color='b')
    ax.plot(C_usages.t, C_usages.exudation / C_usages.C_produced * 100, label=u'C loss by exudation', color='c')
    # ax.plot(C_usages.t, C_usages.C_budget * 100, label=u'C consumption vs. production', color='k')

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'Carbon usages : Photosynthesis (%)')
    ax.set_ylim(bottom=0, top=100.)

    fig.suptitle(u'Total cumulated usages are ' + str(round(C_usages.C_budget.tail(1).values[0] * 100, 0)) + u' % of Photosynthesis')

    plt.savefig(os.path.join(scenario_graphs_dirpath, 'C_usages_cumulated.PNG'), format='PNG', bbox_inches='tight')
    plt.close()


def graph_Conc_Nitrates(scenario, scenario_graphs_dirpath=None, scenario_postprocessing_dirpath=None):
    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_postprocessing_dirpath:
        scenario_postprocessing_dirpath = os.path.join('outputs', scenario_name, 'postprocessing')

    # --- Import simulations prostprocessings
    df_elt = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'elements_postprocessing.csv'))
    df_org = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'organs_postprocessing.csv'))
    df_axe = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'axes_postprocessing.csv'))
    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_roots = df_roots.reset_index(drop=True)

    # --- Conc_Nitrates_shoot
    nitrates_shoot = df_elt.groupby(['t'], as_index=False).agg({'nitrates': 'sum'})
    nitrates_shoot['Conc_Nitrates_DM'] = nitrates_shoot.nitrates / df_axe.sum_dry_mass_shoot

    fig, ax = plt.subplots()
    ax.plot(nitrates_shoot.t, nitrates_shoot.Conc_Nitrates_DM)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'[Nitrates] (µmol g$^{-1}$ DM)')
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'Conc_Nitrates_shoot_DM.PNG'), format='PNG', bbox_inches='tight')
    plt.close()

    # --- Conc_Nitrates_roots
    nitrates_shoot = df_elt.groupby(['t'], as_index=False).agg({'nitrates': 'sum'})
    df_roots['Conc_Nitrates_DM'] = df_roots.nitrates / df_axe.sum_dry_mass_roots

    fig, ax = plt.subplots()
    ax.plot(df_roots.t, df_roots.Conc_Nitrates_DM)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'[Nitrates] (µmol g$^{-1}$ DM)')
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'Conc_Nitrates_roots_DM.PNG'), format='PNG', bbox_inches='tight')
    plt.close()

    # --- Conc_Nitrates_plant
    nitrates_plant = nitrates_shoot.nitrates + df_roots.nitrates
    df_axe['Conc_Nitrates_DM'] = nitrates_plant / df_axe.sum_dry_mass

    fig, ax = plt.subplots()
    ax.plot(df_axe.t, df_axe.Conc_Nitrates_DM)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'[Nitrates] (µmol g$^{-1}$ DM)')
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'Conc_Nitrates_plant_DM.PNG'), format='PNG', bbox_inches='tight')
    plt.close()


def graph_Uptake_N(scenario, scenario_graphs_dirpath=None, scenario_outputs_dirpath=None):
    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_outputs_dirpath:
        scenario_outputs_dirpath = os.path.join('outputs', scenario_name, 'outputs')

    # --- Import simulations outsptus
    df_org = pd.read_csv(os.path.join(scenario_outputs_dirpath, 'organs_outputs.csv'))
    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_roots = df_roots.reset_index(drop=True)

    df_roots['Specific_Uptake_Nitrates'] = df_roots.Uptake_Nitrates / df_roots.mstruct / 3600

    fig, ax = plt.subplots()
    ax.plot(df_roots.t, df_roots.Specific_Uptake_Nitrates)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'N uptake rate (µmol g$^{-1}$ mstruct s$^{-1}$)')
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'Uptake_Nitrates_Specific.PNG'), format='PNG', bbox_inches='tight')
    plt.close()


def graph_phi_s_Devienne1994a(scenario, scenario_graphs_dirpath=None, scenario_outputs_dirpath=None):
    """Sum of NO3- reduction flux and NO3- deposition into the sylem sap"""

    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_outputs_dirpath:
        scenario_outputs_dirpath = os.path.join('outputs', scenario_name, 'outputs')

    # --- Import simulations outsptus
    df_org = pd.read_csv(os.path.join(scenario_outputs_dirpath, 'organs_outputs.csv'))
    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_roots = df_roots.reset_index(drop=True)

    df_roots['phi_s'] = (df_roots.Export_Nitrates / df_roots.mstruct / 3600) + (df_roots.S_Amino_Acids / 3600)

    fig, ax = plt.subplots()
    ax.plot(df_roots.t, df_roots.phi_s)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'NO3- root reduction + NO3- export xylem (µmol N g$^{-1}$ mstruct s$^{-1}$)')
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'phi_s.PNG'), format='PNG', bbox_inches='tight')
    plt.close()


def graph_N_dilution(scenario, scenario_graphs_dirpath=None, scenario_postprocessing_dirpath=None):
    """N shoot in function of above-ground biomass"""

    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_postprocessing_dirpath:
        scenario_postprocessing_dirpath = os.path.join('outputs', scenario_name, 'postprocessing')

    # --- Import simulations postprocessing
    df_axe = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'axes_postprocessing.csv'))

    fig, ax = plt.subplots()
    ax.plot(df_axe.sum_dry_mass_shoot, df_axe.N_content_shoot)
    ax.set_xlabel(u'Shoot biomass (g)')
    ax.set_ylabel(u'N shoot (% DM)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'N_dilution.PNG'), format='PNG', bbox_inches='tight')
    plt.close()


def graph_photosynthetic_rates(scenario, blade_number, scenario_graphs_dirpath=None, scenario_outputs_dirpath=None):
    scenario_name = 'Scenario_%.4d' % scenario
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')
    if not scenario_outputs_dirpath:
        scenario_outputs_dirpath = os.path.join('outputs', scenario_name, 'outputs')

    # import simulation outputs
    df_elt = pd.read_csv(os.path.join(scenario_outputs_dirpath, 'elements_outputs.csv'))
    df_elt = df_elt[(df_elt.axis == 'MS') & (df_elt.element == 'LeafElement1')].copy()

    # ------ Graph per blade_number
    for bl in blade_number:

        df_bl = df_elt[df_elt.metamer == bl]

        fig, (ax1) = plt.subplots(1)

        ax1.plot(df_bl.t, df_bl.Ap, label="Ap")
        ax1.plot(df_bl.t, df_bl.Ac, label="Ac")
        ax1.plot(df_bl.t, df_bl.Aj, label="Aj")
        ax1.plot(df_bl.t, df_bl.Ag, label="Ag")
        # ax1.plot(df_elt.t, df_elt.Ag_before_inhibition_WSC, label="Ag_before_inhibition_WSC")

        # Formatting
        ax1.set_ylabel(u'Photosynthetic rate (µmol m$^{-2}$ s$^{-1}$)')
        ax1.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
        ax1.set_xlabel('t (Hour)')
        plt.tight_layout()
        plt.savefig(os.path.join(scenario_graphs_dirpath, 'Photosynthetic_rates_F' + str(bl) + '.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

def graph_summary(scenario, scenario_graphs_dirpath=None, graph_list=None):
    if graph_list is None:
        graph_list = ['LAI', 'sum_dry_mass_axis', 'shoot_roots_ratio_axis', 'N_content_shoot_axis', 'Conc_Amino_acids_phloem', 'Conc_Sucrose_phloem', 'leaf_Lmax',
                      'green_area_blade']
    scenario_name = 'Scenario_%.4d' % scenario
    scenario_label = scenarios_df['Scenario_label'][scenario]
    if not scenario_graphs_dirpath:
        scenario_graphs_dirpath = os.path.join('outputs', scenario_name, 'graphs')

    nb_graphs = len(graph_list)
    if nb_graphs <= 4:
        nrow = 2
        ncol = 2
    elif nb_graphs <= 6:
        nrow = 2
        ncol = 3
    elif nb_graphs <= 8:
        nrow = 2
        ncol = 4
    elif nb_graphs == 9:
        nrow = 3
        ncol = 3
    else:
        raise AttributeError('Too many graphs to diplay.')

    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 4, nrow * 3))
    plt.suptitle('Scenario ' + str(scenario) + ' ' + scenario_label, y=1)
    i = 0
    for _, ax in np.ndenumerate(axs):
        image = mpimg.imread(os.path.join(scenario_graphs_dirpath, graph_list[i] + '.PNG'))
        ax.imshow(image)
        ax.axis('off')
        i += 1
    plt.tight_layout(pad=0.05)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(top=0.95)
    # plt.show()
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'Summary.PNG'), format='PNG', bbox_inches='tight', dpi=200)
    plt.close()


if __name__ == '__main__':

    # ------- Make the graphs and tables for all the scenarios

    for sc in scenarios:
        # graph_phi_s_Devienne1994a(int(sc))
        graph_Uptake_N(int(sc))
        graph_Conc_Nitrates(int(sc))
        graph_photosynthetic_rates(int(sc), blade_number=[6, 12])
        graph_N_dilution(int(sc))
        graph_RER(int(sc))
        graph_C_usages(int(sc))
        graph_summary(int(sc), graph_list=['LAI', 'sum_dry_mass_axis', 'shoot_roots_ratio_axis', 'N_content_shoot_axis', 'Conc_Amino_acids_phloem', 'Conc_Sucrose_phloem', 'leaf_Lmax',
                                                 'green_area_blade'])
