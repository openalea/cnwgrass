# -*- coding: latin-1 -*-

# PB dans les fichiers de sorties : les EN cachés n'apparaissent que pour les pas de temps impairs (mais les résultats sont corrects)

# -----------
# - PREAMBULE
# -----------

import os
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

random.seed(1234)
np.random.seed(1234)

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

GRAPHS_DIRPATH = 'graphs'

x_name = 't'
x_label = 'Time (Hour)'

# Define culm density (culm m-2)
DENSITY = 410.
NPLANTS = 1
CULM_DENSITY = {i: DENSITY / NPLANTS for i in range(1, NPLANTS + 1)}

FILLING_INIT = 360

TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
SUCROSE_MOLAR_MASS_C_RATIO = 0.42
HEXOSE_MOLAR_MASS_C_RATIO = 0.4
NITRATES_MOLAR_MASS_N_RATIO = 0.23
AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145
phloem_shoot_root = 0.75

# Graphs format

plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize=8)  # 'x-small')
plt.rc('ytick', labelsize=8)  # 'x-small')
plt.rc('axes', labelsize=8, titlesize=10)  #
plt.rc('legend', fontsize=8, frameon=False)  #
plt.rc('lines', markersize=6)

# --------
# - GRAPHS
# --------

fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(3, 3)

# --------
# - H0
# --------

# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = os.path.join('NEMA_H0', 'outputs')
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

# post-processing directory path
POSTPROCESSING_DIRPATH = os.path.join('NEMA_H0', 'postprocessing')
AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
HIDDENZONES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'hiddenzones_postprocessing.csv')
ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

# 1) Photosynthetic organs
ph_elements_output_df = pd.read_csv(ELEMENTS_POSTPROCESSING_FILEPATH)
ph_elements_output_df = ph_elements_output_df[ph_elements_output_df.t % 2 == 1].copy()

# 2) Roots, grains and phloem
organs_output_df = pd.read_csv(ORGANS_POSTPROCESSING_FILEPATH)

# -- OBSERVED DATA

# Data NEMA
t_NEMA = [0, 408, 648, 864, 1200]
dry_mass_ph_NEMA_H0 = [1.95, 2.34, 1.97, 1.73, 1.67]  # Photosynthetic organs i.e. laminae + stems + chaff (g)
dry_mass_ph_NEMA_H0_SD = [0.11, 0.02, 0.08, 0.06, 0.14]
dry_mass_grains_NEMA_H0 = [0.58, 1.14, 1.54, 1.48]  # Grains (g)
dry_mass_grains_NEMA_H0_SD = [0.08, 0.09, 0.10, 0.04]
dry_mass_tot_NEMA_H0 = [1.95, 2.92, 3.11, 3.27, 3.15]
dry_mass_tot_NEMA_H0_SD = [0.11, 0.06, 0.17, 0.16, 0.14]

green_area_lamina1_H0 = [34.6, 32.6, 27.2, 13.3]
green_area_lamina1_H0_SD = [3.8, 2.7, 3.7, 0]
green_area_lamina2_H0 = [34, 32.8, 22.3]
green_area_lamina2_H0_SD = [2.7, 3.2, 2.1]
green_area_lamina3_H0 = [22.8, 22.6, 20.5]
green_area_lamina3_H0_SD = [1.8, 1.3, 0]
green_area_lamina4_H0 = [16, 17.9]
green_area_lamina4_H0_SD = [2.1, 5.1]

green_area_laminae_H0 = [107.4, 105.9, 60, 8]
green_area_laminae_H0_SD = [3.8, 9, 4.5, 0]

N_tot_lamina1_H0 = [6.84, 4.55, 2.29, 0.88, 0.83]
N_tot_lamina1_H0_SD = [0.81, 0.68, 0.51, 0.12, 0.04]
N_tot_lamina2_H0 = [4.08, 2.80, 1.29, 0.62, 0.66]
N_tot_lamina2_H0_SD = [0.34, 0.39, 0.15, 0.03, 0.08]
N_tot_lamina3_H0 = [1.85, 1.18, 0.40, 0.36, 0.41]
N_tot_lamina3_H0_SD = [0.15, 0.17, 0.07, 0.06, 0.05]
N_tot_lamina4_H0 = [0.51, 0.50, 0.40, 0.29, 0.25]
N_tot_lamina4_H0_SD = [0.22, 0.07, 0.08, 0.03, 0.10]
N_tot_chaff_H0 = [5.17, 3.29, 1.81, 1.53, 1.75]
N_tot_chaff_H0_SD = [0.36, 0.26, 0.16, 0.13, 0.73]
N_tot_stem_H0 = [11.47, 8.78, 6.15, 3.44, 2.95]
N_tot_stem_H0_SD = [0.98, 0.60, 0.85, 0.26, 0.29]

N_mass_ph_NEMA_H0 = [29.91, 21.09, 12.34, 7.11, 6.84]  # Photosynthetic organs i.e. laminae + stems + chaff (mg)
N_mass_ph_NEMA_H0_SD = [2.14, 0.87, 1.27, 0.21, 1.17]
N_mass_grains_NEMA_H0 = [9.15, 16.39, 24.68, 26.12]  # Grains (mg)
N_mass_grains_NEMA_H0_SD = [1.18, 1.58, 1.51, 0.33]
N_tot_NEMA_H0 = [29.91, 30.24, 28.73, 31.79, 32.97]
N_tot_NEMA_H0_SD = [2.14, 1.66, 2.63, 1.71, 1.32]

DM_tot_lamina1_H0 = [0.18, 0.18, 0.17, 0.13, 0.13]
DM_tot_lamina1_H0_SD = [0.02, 0.01, 0.01, 0.00, 0.00]
DM_tot_lamina2_H0 = [0.13, 0.13, 0.12, 0.09, 0.08]
DM_tot_lamina2_H0_SD = [0.01, 0.01, 0.00, 0.00, 0.00]
DM_tot_lamina3_H0 = [0.08, 0.08, 0.06, 0.06, 0.05]
DM_tot_lamina3_H0_SD = [0, 0.01, 0.00, 0.00, 0.00]
DM_tot_lamina4_H0 = [0.03, 0.06, 0.05, 0.04, 0.03]
DM_tot_lamina4_H0_SD = [0.01, 0.01, 0.01, 0.00, 0.01]
DM_tot_stem_H0 = [1.27, 1.59, 1.30, 1.11, 1.02]
DM_tot_stem_H0_SD = [0.07, 0.06, 0.04, 0.04, 0.01]
DM_tot_chaff_H0 = [0.26, 0.31, 0.27, 0.29, 0.36]
DM_tot_chaff_H0_SD = [0.02, 0.03, 0.03, 0.02, 0.12]

# 1) Photosynthetic area
ax0 = plt.subplot(gs[0], clip_on=False)

# - Lamina
laminae_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'blade')].groupby('t').agg({'green_area': 'sum'})

ax0.plot(laminae_model.index, laminae_model.green_area * 10000, label=r'Laminae', linestyle='-', color='g')
ax0.errorbar(t_NEMA[:-1], green_area_laminae_H0, yerr=green_area_laminae_H0_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Stem + Chaff
stem_model = ph_elements_output_df[(ph_elements_output_df['organ'] != 'blade') & (ph_elements_output_df['organ'] != 'ear') &
                                   (ph_elements_output_df['element'] == 'StemElement')].groupby('t').agg({'green_area': 'sum'})
chaff_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'ear')].groupby('t').agg({'green_area': 'sum'})
stem_chaff_area = stem_model.green_area.add(chaff_model.green_area, fill_value=0)
ax0.plot(stem_model.index, stem_chaff_area * 10000, label=r'$\sum$ stem + chaff', linestyle='-', color='b')

# - Total aerial
total_green_area = laminae_model.green_area.add(stem_chaff_area, fill_value=0)
ax0.plot(total_green_area.index, total_green_area * 10000, label='Total aerial', linestyle='-', color='r')

# - Formatting
ax0.set_ylim(0, 200)
ax0.set_yticks([0, 50, 100, 150, 200])
ax0.set_xlim(0, 1300)
ax0.axvline(FILLING_INIT, color='k', linestyle='--')
ax0.get_yaxis().set_label_coords(-0.13, 0.5)
ax0.set_ylabel(u'Photosynthetic area (cm$^{2}$)')
plt.setp(ax0.get_xticklabels(), visible=False)

# 4) Total dry mass accumulation
ax3 = plt.subplot(gs[3], sharex=ax0, clip_on=False)

# - Photosynthetic organs
org_ph = ph_elements_output_df.groupby('t').sum()
sum_dry_mass_org_ph = ((org_ph['triosesP'] * 1E-6 * 12) / TRIOSESP_MOLAR_MASS_C_RATIO +
                       (org_ph['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (org_ph['starch'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['fructan'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                       (org_ph['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       (org_ph['proteins'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       [org_ph['mstruct'][1]] * len(org_ph.index))

ax3.plot(sum_dry_mass_org_ph.index, sum_dry_mass_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax3.errorbar(t_NEMA, dry_mass_ph_NEMA_H0, yerr=dry_mass_ph_NEMA_H0_SD, marker='o', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_dry_mass_roots = ((roots['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                      (roots['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                      (roots['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                      roots['mstruct'])

ax3.plot(sum_dry_mass_roots.index, sum_dry_mass_roots, label='Roots', linestyle='-', color='k')

# - Grains
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_dry_mass_grains = grains['Dry_Mass']

ax3.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, label='Grains', linestyle='-', color='y')
ax3.errorbar(t_NEMA[1:], dry_mass_grains_NEMA_H0, yerr=dry_mass_grains_NEMA_H0_SD, label='H3', marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_dry_mass_phloem = ((phloem['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (phloem['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO)

ax3.plot(sum_dry_mass_phloem.index, sum_dry_mass_phloem, label='Phloem', linestyle='-', color='b')
sum_dry_mass_ph_phloem = sum_dry_mass_org_ph.add(sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
ax3.plot(sum_dry_mass_ph_phloem.index, sum_dry_mass_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
sum_dry_mass_roots_phloem = sum_dry_mass_roots.add(sum_dry_mass_phloem * (1 - phloem_shoot_root), fill_value=0)
ax3.plot(sum_dry_mass_roots_phloem.index, sum_dry_mass_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Total aerial
total_dry_mass = sum_dry_mass_org_ph + sum_dry_mass_grains[sum_dry_mass_grains.index % 2 == 1] + sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root
ax3.plot(total_dry_mass.index, total_dry_mass, label='Total aerial', linestyle='-', color='r')
ax3.errorbar(t_NEMA[1:], dry_mass_tot_NEMA_H0[1:], yerr=dry_mass_tot_NEMA_H0_SD[1:], label='H3', marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Formatting
ax3.set_xlim(0, 1300)
ax3.set_ylim(0, 4)
ax3.set_yticks([0, 1, 2, 3])
ax3.axvline(FILLING_INIT, color='k', linestyle='--')
ax3.get_yaxis().set_label_coords(-0.13, 0.5)
ax3.set_ylabel('Dry mass (g)')
plt.setp(ax3.get_xticklabels(), visible=False)

# 5) Total N accumulation
ax6 = plt.subplot(gs[6], sharex=ax0, clip_on=False)

# - Photosynthetic organs
org_ph = ph_elements_output_df.groupby('t').sum()
sum_N_org_ph = ((org_ph['nitrates'] * 1E-3 * 14) +
                (org_ph['amino_acids'] * 1E-3 * 14) +
                (org_ph['proteins'] * 1E-3 * 14) +
                [org_ph['Nstruct'][1] * 1E3] * len(org_ph.index))

ax6.plot(sum_N_org_ph.index, sum_N_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax6.errorbar(t_NEMA, N_mass_ph_NEMA_H0, yerr=N_mass_ph_NEMA_H0_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_N_roots = ((roots['nitrates'] * 1E-3 * 14) +
               (roots['amino_acids'] * 1E-3 * 14) +
               (roots['Nstruct'] * 1E3))

ax6.plot(sum_N_roots.index, sum_N_roots, label='Roots', linestyle='-', color='k')

# - Grains
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_N_grains = grains['proteins'] * 1E-3 * 14

ax6.plot(sum_N_grains.index, sum_N_grains, label='Grains', linestyle='-', color='y')
ax6.errorbar(t_NEMA[1:], N_mass_grains_NEMA_H0, yerr=N_mass_grains_NEMA_H0_SD, label='H3', marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_N_phloem = phloem['amino_acids'] * 1E-3 * 14
sum_N_ph_phloem = sum_N_org_ph.add(sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
sum_N_roots_phloem = sum_N_roots.add(sum_N_phloem * (1 - phloem_shoot_root), fill_value=0)
ax6.plot(sum_N_ph_phloem.index, sum_N_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
ax6.plot(sum_N_phloem.index, sum_N_phloem, label='Phloem', linestyle='-', color='b')
ax6.plot(sum_N_roots_phloem.index, sum_N_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Total aerial
total_N_mass = sum_N_org_ph + sum_N_grains[sum_N_grains.index % 2 == 1] + sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root
ax6.plot(total_N_mass.index, total_N_mass, label='Total aerial', linestyle='-', color='r')
ax6.errorbar(t_NEMA[1:], N_tot_NEMA_H0[1:], yerr=N_tot_NEMA_H0_SD[1:], label='H3', marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Formatting
ax6.set_xlim(0, 1300)
ax6.set_xticks([0, 200, 400, 600, 800, 1000, 1200])
ax6.set_ylim(0, 60)
ax6.set_yticks([0, 15, 30, 45])
# ax6.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
ax6.axvline(FILLING_INIT, color='k', linestyle='--')
ax6.set_xlabel('Time from flowering (h)')
ax6.get_yaxis().set_label_coords(-0.13, 0.5)
ax6.set_ylabel('N mass (mg)')
xticks = ax6.xaxis.get_major_ticks()
xticks[-1].label1.set_visible(False)

# ------
# --- H3
# ------


# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = os.path.join('NEMA_H3', 'outputs')
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

# post-processing directory path
POSTPROCESSING_DIRPATH = os.path.join('NEMA_H3', 'postprocessing')
AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
HIDDENZONES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'hiddenzones_postprocessing.csv')
ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

# 1) Photosynthetic organs
ph_elements_output_df = pd.read_csv(ELEMENTS_POSTPROCESSING_FILEPATH)
ph_elements_output_df = ph_elements_output_df[ph_elements_output_df.t % 2 == 1].copy()

# 2) Roots, grains and phloem
organs_output_df = pd.read_csv(ORGANS_POSTPROCESSING_FILEPATH)

# - -- OBSERVED DATA


TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
SUCROSE_MOLAR_MASS_C_RATIO = 0.42
HEXOSE_MOLAR_MASS_C_RATIO = 0.4
NITRATES_MOLAR_MASS_N_RATIO = 0.23
AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145
phloem_shoot_root = 0.75

# Data NEMA
t_NEMA = [0, 408, 648, 864, 1200]
dry_mass_ph_NEMA_H3 = [1.95, 2.30, 2.07, 1.75, 1.57]  # Photosynthetic organs i.e. laminae + stems + chaff (g)
dry_mass_ph_NEMA_H3_SD = [0.11, 0.10, 0.02, 0.03, 0.09]
dry_mass_grains_NEMA_H3 = [0.53, 1.22, 1.46, 1.54]
dry_mass_grains_NEMA_H3_SD = [0.03, 0.09, 0.10, 0.06]
dry_mass_tot_NEMA_H3 = [1.95, 2.83, 3.29, 3.21, 3.11]
dry_mass_tot_NEMA_H3_SD = [0.11, 0.11, 0.11, 0.12, 0.12]

green_area_lamina1_H3 = [34.6, 36.6, 36.9, 18.1]
green_area_lamina1_H3_SD = [3.8, 2.0, 2.4, 1.0]
green_area_lamina2_H3 = [34, 36.5, 28.9, 18.1]
green_area_lamina2_H3_SD = [2.7, 1.2, 0.8, 0]
green_area_lamina3_H3 = [22.8, 24.0, 15.8]
green_area_lamina3_H3_SD = [1.8, 1.4, 1.7]
green_area_lamina4_H3 = [16, 16.3, 13.1]
green_area_lamina4_H3_SD = [2.1, 1.6, 0]

green_area_laminae_H3 = [107.4, 113.4, 94, 25]
green_area_laminae_H3_SD = [3.8, 2, 2.4, 9]

N_tot_lamina1_H3 = [6.84, 6.83, 4.54, 1.37, 0.97]
N_tot_lamina1_H3_SD = [0.81, 0.31, 0.21, 0.21, 0.16]
N_tot_lamina2_H3 = [4.08, 4.13, 2.04, 1.01, 0.77]
N_tot_lamina2_H3_SD = [0.34, 0.28, 0.07, 0.16, 0.06]
N_tot_lamina3_H3 = [1.85, 1.77, 0.59, 0.55, 0.47]
N_tot_lamina3_H3_SD = [0.15, 0.15, 0.05, 0.03, 0.02]
N_tot_lamina4_H3 = [0.51, 0.63, 0.37, 0.34, 0.34]
N_tot_lamina4_H3_SD = [0.22, 0.03, 0.02, 0.05, 0.02]
N_tot_chaff_H3 = [5.17, 3.69, 2.44, 1.78, 1.01]
N_tot_chaff_H3_SD = [0.36, 0.12, 0.14, 0.09, 0.33]
N_tot_stem_H3 = [11.47, 10.37, 7.65, 4.12, 3.75]
N_tot_stem_H3_SD = [0.98, 1.38, 0.66, 0.49, 0.30]

N_mass_ph_NEMA_H3 = [29.91, 27.42, 17.62, 9.18, 7.32]  # Photosynthetic organs i.e. laminae + stems + chaff (mg)
N_mass_ph_NEMA_H3_SD = [2.14, 2.03, 0.38, 0.79, 0.58]
N_mass_grains_NEMA_H3 = [10.01, 22.24, 30.83, 32.84]  # Grains (mg)
N_mass_grains_NEMA_H3_SD = [0.76, 2.36, 2.28, 2.67]
N_tot_NEMA_H3 = [29.91, 37.43, 39.87, 40.01, 40.16]
N_tot_NEMA_H3_SD = [2.14, 2.35, 2.48, 1.49, 3.14]

DM_tot_lamina1_H3 = [0.18, 0.20, 0.20, 0.16, 0.14]
DM_tot_lamina1_H3_SD = [0.02, 0.01, 0.01, 0.00, 0.00]
DM_tot_lamina2_H3 = [0.13, 0.13, 0.12, 0.10, 0.09]
DM_tot_lamina2_H3_SD = [0.01, 0.01, 0.01, 0.00, 0.00]
DM_tot_lamina3_H3 = [0.08, 0.08, 0.06, 0.06, 0.05]
DM_tot_lamina3_H3_SD = [0, 0.00, 0.01, 0.00, 0.00]
DM_tot_lamina4_H3 = [0.03, 0.05, 0.04, 0.04, 0.04]
DM_tot_lamina4_H3_SD = [0.01, 0.00, 0.00, 0.00, 0.00]
DM_tot_stem_H3 = [1.27, 1.52, 1.33, 1.11, 1.04]
DM_tot_stem_H3_SD = [0.07, 0.08, 0.02, 0.02, 0.02]
DM_tot_chaff_H3 = [0.26, 0.32, 0.31, 0.29, 0.21]
DM_tot_chaff_H3_SD = [0.02, 0.01, 0.01, 0.00, 0.07]

# 1) Photosynthetic area
ax1 = plt.subplot(gs[1], sharey=ax0, clip_on=False)

# - Lamina
laminae_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'blade')].groupby('t').agg({'green_area': 'sum'})

ax1.plot(laminae_model.index, laminae_model.green_area * 10000, label=r'Laminae', linestyle='-', color='g')
ax1.errorbar(t_NEMA[:-1], green_area_laminae_H3, yerr=green_area_laminae_H3_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Stem + Chaff
stem_model = ph_elements_output_df[(ph_elements_output_df['organ'] != 'blade') & (ph_elements_output_df['organ'] != 'ear') &
                                   (ph_elements_output_df['element'] == 'StemElement')].groupby('t').agg({'green_area': 'sum'})
chaff_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'ear')].groupby('t').agg({'green_area': 'sum'})
stem_chaff_area = stem_model.green_area.add(chaff_model.green_area, fill_value=0)
ax1.plot(stem_model.index, stem_chaff_area * 10000, label=r'$\sum$ stem + chaff', linestyle='-', color='b')

# - Total aerial
total_green_area = laminae_model.green_area.add(stem_chaff_area, fill_value=0)
ax1.plot(total_green_area.index, total_green_area * 10000, label='Total aerial', linestyle='-', color='r')

# - Formatting
ax1.set_xlim(0, 1300)
ax1.set_ylim(0, 200)
ax1.set_yticks([0, 50, 100, 150, 200])
ax1.set_xlim(0, 1300)
ax1.axvline(FILLING_INIT, color='k', linestyle='--')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)

# 4) Total dry mass accumulation
ax4 = plt.subplot(gs[4], sharex=ax1, sharey=ax3, clip_on=False)

# - Photosynthetic organs
org_ph = ph_elements_output_df.groupby('t').sum()
sum_dry_mass_org_ph = ((org_ph['triosesP'] * 1E-6 * 12) / TRIOSESP_MOLAR_MASS_C_RATIO +
                       (org_ph['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (org_ph['starch'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['fructan'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                       (org_ph['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       (org_ph['proteins'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       [org_ph['mstruct'][1]] * len(org_ph.index))

ax4.plot(sum_dry_mass_org_ph.index, sum_dry_mass_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax4.errorbar(t_NEMA, dry_mass_ph_NEMA_H3, yerr=dry_mass_ph_NEMA_H3_SD, marker='o', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_dry_mass_roots = ((roots['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                      (roots['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                      (roots['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                      roots['mstruct'])

ax4.plot(sum_dry_mass_roots.index, sum_dry_mass_roots, label='Roots', linestyle='-', color='k')

# - Grains
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_dry_mass_grains = grains['Dry_Mass']

ax4.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, label='Grains', linestyle='-', color='y')
ax4.errorbar(t_NEMA[1:], dry_mass_grains_NEMA_H3, yerr=dry_mass_grains_NEMA_H3_SD, label='H3', marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_dry_mass_phloem = ((phloem['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (phloem['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO)

ax4.plot(sum_dry_mass_phloem.index, sum_dry_mass_phloem, label='Phloem', linestyle='-', color='b')
sum_dry_mass_ph_phloem = sum_dry_mass_org_ph.add(sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
ax4.plot(sum_dry_mass_ph_phloem.index, sum_dry_mass_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
sum_dry_mass_roots_phloem = sum_dry_mass_roots.add(sum_dry_mass_phloem * (1 - phloem_shoot_root), fill_value=0)
ax4.plot(sum_dry_mass_roots_phloem.index, sum_dry_mass_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Total aerial
total_dry_mass = sum_dry_mass_org_ph + sum_dry_mass_grains[sum_dry_mass_grains.index % 2 == 1] + sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root
ax4.plot(total_dry_mass.index, total_dry_mass, label='Total aerial', linestyle='-', color='r')
ax4.errorbar(t_NEMA[1:], dry_mass_tot_NEMA_H3[1:], yerr=dry_mass_tot_NEMA_H3_SD[1:], label='H3', marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Formatting
ax4.set_xlim(0, 1300)
ax4.set_ylim(0, 4)
ax4.set_yticks([0, 1, 2, 3])
ax4.axvline(FILLING_INIT, color='k', linestyle='--')
ax4.get_yaxis().set_label_coords(-0.13, 0.5)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)

# 5) Total N accumulation
ax7 = plt.subplot(gs[7], sharex=ax1, sharey=ax6, clip_on=False)

# - Photosynthetic organs
org_ph = ph_elements_output_df.groupby('t').sum()
sum_N_org_ph = ((org_ph['nitrates'] * 1E-3 * 14) +
                (org_ph['amino_acids'] * 1E-3 * 14) +
                (org_ph['proteins'] * 1E-3 * 14) +
                [org_ph['Nstruct'][1] * 1E3] * len(org_ph.index))

ax7.plot(sum_N_org_ph.index, sum_N_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax7.errorbar(t_NEMA, N_mass_ph_NEMA_H3, yerr=N_mass_ph_NEMA_H3_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_N_roots = ((roots['nitrates'] * 1E-3 * 14) +
               (roots['amino_acids'] * 1E-3 * 14) +
               (roots['Nstruct'] * 1E3))

ax7.plot(sum_N_roots.index, sum_N_roots, label='Roots', linestyle='-', color='k')

# - Grains
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_N_grains = grains['proteins'] * 1E-3 * 14

ax7.plot(sum_N_grains.index, sum_N_grains, label='Grains', linestyle='-', color='y')
ax7.errorbar(t_NEMA[1:], N_mass_grains_NEMA_H3, yerr=N_mass_grains_NEMA_H3_SD, label='H3', marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_N_phloem = phloem['amino_acids'] * 1E-3 * 14
sum_N_ph_phloem = sum_N_org_ph.add(sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
sum_N_roots_phloem = sum_N_roots.add(sum_N_phloem * (1 - phloem_shoot_root), fill_value=0)
ax7.plot(sum_N_ph_phloem.index, sum_N_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
ax7.plot(sum_N_phloem.index, sum_N_phloem, label='Phloem', linestyle='-', color='b')
ax7.plot(sum_N_roots_phloem.index, sum_N_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Total aerial
total_N_mass = sum_N_org_ph + sum_N_grains[sum_N_grains.index % 2 == 1] + sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root
ax7.plot(total_N_mass.index, total_N_mass, label='Total aerial', linestyle='-', color='r')
ax7.errorbar(t_NEMA[1:], N_tot_NEMA_H3[1:], yerr=N_tot_NEMA_H3_SD[1:], label='H3', marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Formatting
ax7.set_xlim(0, 1300)
ax7.set_xticks([0, 200, 400, 600, 800, 1000, 1200])
ax7.set_ylim(0, 60)
ax7.set_yticks([0, 15, 30, 45])
# ax7.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
ax7.axvline(FILLING_INIT, color='k', linestyle='--')
ax7.set_xlabel('Time from flowering (h)')
xticks = ax7.xaxis.get_major_ticks()
xticks[-1].label1.set_visible(False)
plt.setp(ax7.get_yticklabels(), visible=False)

# -------
# --- H15
# -------


# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = os.path.join('NEMA_H15', 'outputs')
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

# post-processing directory path
POSTPROCESSING_DIRPATH = os.path.join('NEMA_H15', 'postprocessing')
AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
HIDDENZONES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'hiddenzones_postprocessing.csv')
ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

# 1) Photosynthetic organs
ph_elements_output_df = pd.read_csv(ELEMENTS_POSTPROCESSING_FILEPATH)
ph_elements_output_df = ph_elements_output_df[ph_elements_output_df.t % 2 == 1].copy()

# 2) Roots, grains and phloem
organs_output_df = pd.read_csv(ORGANS_POSTPROCESSING_FILEPATH)

# - -- OBSERVED DATA

TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
SUCROSE_MOLAR_MASS_C_RATIO = 0.42
HEXOSE_MOLAR_MASS_C_RATIO = 0.4
NITRATES_MOLAR_MASS_N_RATIO = 0.23
AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145
phloem_shoot_root = 0.75

# Data NEMA
t_NEMA = [0, 408, 648, 864, 1200]
dry_mass_ph_NEMA_H15 = [1.95, 2.30, 2.00, 1.71, 1.55]  # Photosynthetic organs i.e. laminae + stems + chaff (g)
dry_mass_ph_NEMA_H15_SD = [0.11, 0.02, 0.03, 0.05, 0.24]
dry_mass_grains_NEMA_H15 = [0.58, 1.29, 1.41, 1.50]
dry_mass_grains_NEMA_H15_SD = [0.02, 0.10, 0.13, 0.06]
dry_mass_tot_NEMA_H15 = [1.95, 2.88, 3.29, 3.12, 3.05]
dry_mass_tot_NEMA_H15_SD = [0.11, 0.01, 0.07, 0.15, 0.20]

green_area_lamina1_H15 = [34.6, 34.7, 35.1, 22.3]
green_area_lamina1_H15_SD = [3.8, 2.1, 1.0, 3.5]
green_area_lamina2_H15 = [34, 33.5, 31.4, 21.4]
green_area_lamina2_H15_SD = [2.7, 0.3, 1.8, 0.7]
green_area_lamina3_H15 = [22.8, 23.0, 18.8]
green_area_lamina3_H15_SD = [1.8, 0.5, 1.3]
green_area_lamina4_H15 = [16, 16.8, 13.8]
green_area_lamina4_H15_SD = [2.1, 0.1, 0.2]

green_area_laminae_H15 = [107.4, 108, 99.1, 43.7]
green_area_laminae_H15_SD = [3.8, 2.1, 7, 12]

N_tot_lamina1_H15 = [6.84, 7.35, 5.36, 1.94, 1.12]
N_tot_lamina1_H15_SD = [0.81, 0.34, 0.13, 0.03, 0.10]
N_tot_lamina2_H15 = [4.08, 4.53, 3.17, 1.17, 0.94]
N_tot_lamina2_H15_SD = [0.34, 0.15, 0.18, 0.09, 0.01]
N_tot_lamina3_H15 = [1.85, 2.14, 1.00, 0.73, 0.61]
N_tot_lamina3_H15_SD = [0.15, 0.14, 0.10, 0.06, 0.05]
N_tot_lamina4_H15 = [0.51, 0.80, 0.48, 0.46, 0.45]
N_tot_lamina4_H15_SD = [0.22, 0.16, 0.09, 0.07, 0.02]
N_tot_chaff_H15 = [5.17, 3.67, 2.38, 1.57, 0.97]
N_tot_chaff_H15_SD = [0.36, 0.28, 0.13, 0.09, 1.02]
N_tot_stem_H15 = [11.47, 12.01, 9.66, 5.94, 4.85]
N_tot_stem_H15_SD = [0.98, 1.01, 0.89, 0.31, 0.69]

N_mass_ph_NEMA_H15 = [29.91, 30.49, 22.06, 11.81, 8.94]  # Photosynthetic organs i.e. laminae + stems + chaff (mg)
N_mass_ph_NEMA_H15_SD = [2.14, 1.06, 1.22, 0.55, 1.64]
N_mass_grains_NEMA_H15 = [11.10, 25.82, 32.17, 36.33]  # Grains (mg)
N_mass_grains_NEMA_H15_SD = [0.20, 3.27, 2.86, 1.57]
N_tot_NEMA_H15 = [29.91, 41.59, 47.88, 43.99, 45.28]
N_tot_NEMA_H15_SD = [2.14, 0.89, 2.46, 2.68, 0.38]

DM_tot_lamina1_H15 = [0.18, 0.20, 0.19, 0.16, 0.14]
DM_tot_lamina1_H15_SD = [0.02, 0.01, 0.00, 0.00, 0.01]
DM_tot_lamina2_H15 = [0.13, 0.13, 0.12, 0.10, 0.09]
DM_tot_lamina2_H15_SD = [0.01, 0.00, 0.00, 0.00, 0.00]
DM_tot_lamina3_H15 = [0.08, 0.08, 0.07, 0.06, 0.05]
DM_tot_lamina3_H15_SD = [0, 0.00, 0.01, 0.00, 0.00]
DM_tot_lamina4_H15 = [0.03, 0.05, 0.04, 0.04, 0.04]
DM_tot_lamina4_H15_SD = [0.01, 0.00, 0.01, 0.01, 0.00]
DM_tot_stem_H15 = [1.27, 1.51, 1.28, 1.10, 1.06]
DM_tot_stem_H15_SD = [0.07, 0.00, 0.04, 0.04, 0.05]
DM_tot_chaff_H15 = [0.26, 0.33, 0.30, 0.26, 0.17]
DM_tot_chaff_H15_SD = [0.02, 0.02, 0.00, 0.02, 0.19]

# 1) Photosynthetic area
ax2 = plt.subplot(gs[2], sharey=ax1, clip_on=False)

laminae_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'blade')].groupby('t').agg({'green_area': 'sum'})
stem_model = ph_elements_output_df[(ph_elements_output_df['organ'] != 'blade') & (ph_elements_output_df['organ'] != 'ear') &
                                   (ph_elements_output_df['element'] == 'StemElement')].groupby('t').agg({'green_area': 'sum'})
chaff_model = ph_elements_output_df[(ph_elements_output_df['organ'] == 'ear')].groupby('t').agg({'green_area': 'sum'})
stem_chaff_area = stem_model.green_area.add(chaff_model.green_area, fill_value=0)
total_green_area = laminae_model.green_area.add(stem_chaff_area, fill_value=0)

# - Total aerial
ax2.plot(total_green_area.index, total_green_area * 10000, label='Total aerial', linestyle='-', color='r')

# - Lamina
ax2.plot(laminae_model.index, laminae_model.green_area * 10000, label=r'Laminae', linestyle='-', color='g')
ax2.errorbar(t_NEMA[:-1], green_area_laminae_H15, yerr=green_area_laminae_H15_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Stem + Chaff
ax2.plot(stem_model.index, stem_chaff_area * 10000, label=r'$\sum$ stem + chaff', linestyle='-', color='b')

# - Formatting
ax2.set_xlim(0, 1300)
ax2.set_ylim(0, 200)
ax2.set_yticks([0, 50, 100, 150, 200])
ax2.set_xlim(0, 1300)
ax2.axvline(FILLING_INIT, color='k', linestyle='--')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

# 4) Total dry mass accumulation
ax5 = plt.subplot(gs[5], sharex=ax2, sharey=ax4, clip_on=False)

org_ph = ph_elements_output_df.groupby('t').sum()
sum_dry_mass_org_ph = ((org_ph['triosesP'] * 1E-6 * 12) / TRIOSESP_MOLAR_MASS_C_RATIO +
                       (org_ph['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (org_ph['starch'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['fructan'] * 1E-6 * 12) / HEXOSE_MOLAR_MASS_C_RATIO +
                       (org_ph['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                       (org_ph['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       (org_ph['proteins'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                       [org_ph['mstruct'][1]] * len(org_ph.index))

roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_dry_mass_roots = ((roots['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                      (roots['nitrates'] * 1E-6 * 14) / NITRATES_MOLAR_MASS_N_RATIO +
                      (roots['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                      roots['mstruct'])
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_dry_mass_grains = grains['Dry_Mass']

phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_dry_mass_phloem = ((phloem['sucrose'] * 1E-6 * 12) / SUCROSE_MOLAR_MASS_C_RATIO +
                       (phloem['amino_acids'] * 1E-6 * 14) / AMINO_ACIDS_MOLAR_MASS_N_RATIO)

total_dry_mass = sum_dry_mass_org_ph + sum_dry_mass_grains[sum_dry_mass_grains.index % 2 == 1] + sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root

# - Total aerial
ax5.plot(total_dry_mass.index, total_dry_mass, label='Total aerial', linestyle='-', color='r')
ax5.errorbar(t_NEMA[1:], dry_mass_tot_NEMA_H15[1:], yerr=dry_mass_tot_NEMA_H15_SD[1:], marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Photosynthetic organs
ax5.plot(sum_dry_mass_org_ph.index, sum_dry_mass_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax5.errorbar(t_NEMA, dry_mass_ph_NEMA_H15, yerr=dry_mass_ph_NEMA_H15_SD, marker='o', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
ax5.plot(sum_dry_mass_roots.index, sum_dry_mass_roots, label='Roots', linestyle='-', color='k')

# - Grains
ax5.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, label='Grains', linestyle='-', color='y')
ax5.errorbar(t_NEMA[1:], dry_mass_grains_NEMA_H15, yerr=dry_mass_grains_NEMA_H15_SD, marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
ax5.plot(sum_dry_mass_phloem.index, sum_dry_mass_phloem, label='Phloem', linestyle='-', color='b')
sum_dry_mass_ph_phloem = sum_dry_mass_org_ph.add(sum_dry_mass_phloem[sum_dry_mass_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
ax5.plot(sum_dry_mass_ph_phloem.index, sum_dry_mass_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
sum_dry_mass_roots_phloem = sum_dry_mass_roots.add(sum_dry_mass_phloem * (1 - phloem_shoot_root), fill_value=0)
ax5.plot(sum_dry_mass_roots_phloem.index, sum_dry_mass_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Formatting
ax5.set_xlim(0, 1300)
ax5.set_ylim(0, 4)
ax5.set_yticks([0, 1, 2, 3])
ax5.axvline(FILLING_INIT, color='k', linestyle='--')
ax5.get_yaxis().set_label_coords(-0.13, 0.5)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)

# 5) Total N accumulation
ax8 = plt.subplot(gs[8], sharex=ax2, sharey=ax7, clip_on=False)

# - Photosynthetic organs
org_ph = ph_elements_output_df.groupby('t').sum()
sum_N_org_ph = ((org_ph['nitrates'] * 1E-3 * 14) +
                (org_ph['amino_acids'] * 1E-3 * 14) +
                (org_ph['proteins'] * 1E-3 * 14) +
                [org_ph['Nstruct'][1] * 1E3] * len(org_ph.index))

ax8.plot(sum_N_org_ph.index, sum_N_org_ph, label=r'$\sum$ (tp,i)', linestyle='-', color='g')
ax8.errorbar(t_NEMA, N_mass_ph_NEMA_H15, yerr=N_mass_ph_NEMA_H15_SD, marker='s', color='g', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Roots
roots = organs_output_df[organs_output_df['organ'] == 'roots'].groupby('t').sum()
sum_N_roots = ((roots['nitrates'] * 1E-3 * 14) +
               (roots['amino_acids'] * 1E-3 * 14) +
               (roots['Nstruct'] * 1E3))

ax8.plot(sum_N_roots.index, sum_N_roots, label='Roots', linestyle='-', color='k')

# - Grains
grains = organs_output_df[organs_output_df['organ'] == 'grains'].groupby('t').sum()
sum_N_grains = grains['proteins'] * 1E-3 * 14

ax8.plot(sum_N_grains.index, sum_N_grains, label='Grains', linestyle='-', color='y')
ax8.errorbar(t_NEMA[1:], N_mass_grains_NEMA_H15, yerr=N_mass_grains_NEMA_H15_SD, label='H15', marker='s', color='y', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Phloem
phloem = organs_output_df[organs_output_df['organ'] == 'phloem'].groupby('t').sum()
sum_N_phloem = phloem['amino_acids'] * 1E-3 * 14
sum_N_ph_phloem = sum_N_org_ph.add(sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root, fill_value=0)
sum_N_roots_phloem = sum_N_roots.add(sum_N_phloem * (1 - phloem_shoot_root), fill_value=0)
ax8.plot(sum_N_ph_phloem.index, sum_N_ph_phloem, label=r'$\sum$ (tp,i) + phloem', linestyle='--', color='g')
ax8.plot(sum_N_phloem.index, sum_N_phloem, label='Phloem', linestyle='-', color='b')
ax8.plot(sum_N_roots_phloem.index, sum_N_roots_phloem, label=r'$\sum$ roots + phloem', linestyle='--', color='k')

# - Total aerial
total_N_mass = sum_N_org_ph + sum_N_grains[sum_N_grains.index % 2 == 1] + sum_N_phloem[sum_N_phloem.index % 2 == 1] * phloem_shoot_root
ax8.plot(total_N_mass.index, total_N_mass, label='Total aerial', linestyle='-', color='r')
ax8.errorbar(t_NEMA[1:], N_tot_NEMA_H15[1:], yerr=N_tot_NEMA_H15_SD[1:], label='H15', marker='s', color='r', ecolor='k', linestyle='', clip_on=False, zorder=10)

# - Formatting
ax8.set_xlim(0, 1300)
ax8.set_ylim(0, 60)
ax8.set_yticks([0, 15, 30, 45])
# ax8.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
ax8.axvline(FILLING_INIT, color='k', linestyle='--')
ax8.set_xlabel('Time from flowering (h)')
plt.setp(ax8.get_yticklabels(), visible=False)

# -
# - --- Global formatting
# -

# remove vertical and horizontal gap between subplots
plt.subplots_adjust(hspace=.0, wspace=.0)

# Titles
plt.setp(ax0, title='H0')
plt.setp(ax1, title='H3')
plt.setp(ax2, title='H15')

# letters by subplots
ax0.text(0.1, 0.9, 'A', ha='center', va='center', size=9, transform=ax0.transAxes)
ax1.text(0.1, 0.9, 'B', ha='center', va='center', size=9, transform=ax1.transAxes)
ax2.text(0.1, 0.9, 'C', ha='center', va='center', size=9, transform=ax2.transAxes)
ax3.text(0.1, 0.9, 'D', ha='center', va='center', size=9, transform=ax3.transAxes)
ax4.text(0.1, 0.9, 'E', ha='center', va='center', size=9, transform=ax4.transAxes)
ax5.text(0.1, 0.9, 'F', ha='center', va='center', size=9, transform=ax5.transAxes)
ax6.text(0.1, 0.9, 'G', ha='center', va='center', size=9, transform=ax6.transAxes)
ax7.text(0.1, 0.9, 'H', ha='center', va='center', size=9, transform=ax7.transAxes)
ax8.text(0.1, 0.9, 'I', ha='center', va='center', size=9, transform=ax8.transAxes)

# legends
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)

ax5.legend(loc='lower left', bbox_to_anchor=(1, -0.4), frameon=True)

plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Comparison_sim_obs_3treat.PNG'), format='PNG', bbox_inches='tight', dpi=600)
plt.close()
