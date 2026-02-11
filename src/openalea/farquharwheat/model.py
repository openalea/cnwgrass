# -*- coding: latin-1 -*-

from __future__ import division  # use '//' to do integer division
from math import sqrt, log, exp
from farquharwheat import parameters

"""
    farquharwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of photosynthesis based on Farquhar's approach.
    The model includes the dependence of photosynthesis to organ temperature and nitrogen content.
    Internal CO2 and organ temperature are found numerically.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""


# TODO: extract all parameters and put them in farqhuar.parameters

def _organ_temperature(w, z, Zh, Ur, PAR, gsw, Ta, Ts, RH, organ_name):
    """
    Energy balance for the estimation of organ temperature

    :param float w: organ characteristic dimension (m) to be considered for heat transfer through forced convection (by wind).
             For a leaf: its width (more related to wind direction than length), for cylindric stem elements: diameter.
    :param float z: organ height from soil (m)
    :param float Zh: canopy height (m)
    :param float Ur: wind speed (m s-1) at the reference height (zr), e.g. top of the canopy + 2m (in the case of wheat, Ur can be approximated as the wind speed at 2m from soil)
    :param float PAR: absorbed PAR (µmol m-2 s-1)
    :param float gsw: stomatal conductance to water vapour (mol m-2 s-1)
    :param float Ta: air temperature (degree C)
    :param float Ts: organ temperature (degree C). Ts = Ta at the first iteration of the numeric resolution
    :param float RH: Relative humidity (decimal fraction)
    :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from cylindric organs)

    :return: Ts (organ temperature, degree C), Tr (organ transpiration rate, mm s-1)
    :rtype: (float, float)
    """

    d = parameters.Zh_d * Zh  #: Zero plane displacement height (m) #TODO: should we adapt this calculation before heading ?
    Zo = parameters.Zh_Zo * Zh  #: Roughness length (m)

    # TODO: Temporary patch to avoid div 0 error
    Ur = max(Ur, parameters.Ur_min)

    #: Wind speed
    u_star = (Ur * parameters.K) / log((parameters.ZR - d) / Zo)  #: Friction velocity (m s-1)
    Uh = (u_star / parameters.K) * log((Zh - d) / Zo)  #: Wind speed at the top of canopy (m s-1)
    u = Uh * exp(parameters.A * (z / Zh - 1))  #: Wind speed at organ height (m s-1), from Campbell and Norman (1998), second edition.

    #: Boundary layer resistance to heat (s m-1). See Finnigan J, Raupach M. 1987 and Monteith JL. 1973 for basic equations.
    if organ_name == 'blade':
        rbh = parameters.rhb_blade_A * sqrt(w / u)  #: Case of horizontal planes submitted to forced convection
    else:
        rbh = w / (parameters.rhb_other_A * ((u * w) / parameters.rhb_other_B) ** parameters.rhb_other_C)  #: Case of vertical cylinders submitted to forced convection

    #: Turbulence resistance to heat (s m-1)
    ra = 1 / (parameters.K ** parameters.ra_expo * Ur) * (log((parameters.ZR - d) / Zo)) ** parameters.ra_expo  #: Aerodynamic resistance integrated from zr to z0 + d

    #: Net absorbed radiation Rn (PAR and NIR, J m-2 s-1)
    RGa = (PAR * parameters.PARa_to_RGa) / parameters.Watt_to_PPFD  #: Global absorbed radiation by organ (J m-2 s-1).
    es_Ta = parameters.s_C * exp((parameters.s_B * Ta) / (parameters.s_A + Ta))  #: Saturated vapour pressure of the air (kPa), Ta in degree Celsius
    V = RH * es_Ta  #: Vapour pressure of the air (kPa)
    # fvap = 0.56 - 0.079*sqrt(10*V)                      #: Fraction of vapour pressure
    #
    # tau = RGa/parameters.I0                                    #: Atmospheric transmissivity (dimensionless)
    # fclear = 0.1 + 0.9*max(0, min(1, (tau-0.2)/0.5))    #: Fraction sky clearness

    Rn = RGa  # NB: this only accounts for the visible radiations.
    # General equation is Rn = RGa + epsilon*Ra + epsilon*sigma*(Ts_feuilles_voisines + parameters.KELVIN_DEGREE)**4 - epsilon*sigma*(Ts + parameters.KELVIN_DEGREE)**4
    # if Ra unavailable, use Ra = sigma*(Tair + parameters.KELVIN_DEGREE)**4*fvap*fclear

    #: Transpiration (mm s-1), Penman-Monteith
    if Ts == Ta:
        Ta_K = Ta + parameters.KELVIN_DEGREE
        s = ((parameters.s_B * parameters.s_A) / (Ta_K + parameters.s_A) ** parameters.s_expo) * es_Ta  #: Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)
    else:
        es_Tl = parameters.s_C * exp((parameters.s_B * Ts) / (parameters.s_A + Ts))  #: Saturated vapour pressure at organ level (kPa), Ts in degree Celsius
        Ts_K, Ta_K = Ts + parameters.KELVIN_DEGREE, Ta + parameters.KELVIN_DEGREE
        s = (es_Tl - es_Ta) / (Ts_K - Ta_K)  #: Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)

    VPDa = es_Ta - V
    rbw = parameters.rbh_rbw * rbh  #: Boundary layer resistance for water (s m-1)
    gsw_physic = (gsw * parameters.R * (Ts + parameters.KELVIN_DEGREE)) / parameters.PATM  #: Stomatal conductance to water in physical units (m s-1). Relation given by A. Tuzet (2003)
    rswp = 1 / gsw_physic  #: Stomatal resistance for water (s m-1)
    Tr = max(0., (s * Rn + (parameters.RHOCP * VPDa) / (rbh + ra)) / (parameters.LAMBDA * (s + parameters.GAMMA * ((rbw + ra + rswp) / (rbh + ra)))))  #: mm s-1

    #: Organ temperature
    Ts = Ta + ((rbh + ra) * (Rn - parameters.LAMBDA * Tr)) / parameters.RHOCP

    return Ts, Tr


def _stomatal_conductance(Ag, An, surfacic_nitrogen, ambient_CO2, RH):
    """
    Ball, Woodrow, and Berry model of stomatal conductance (1987)

    :param float Ag: gross assimilation rate (µmol m-2 s-1)
    :param float An: net assimilation rate (µmol m-2 s-1)
    :param float surfacic_nitrogen: surfacic nitrogen content(g m-2) including or not structural nitrogen depending on parameter.MODEL_VERSION
    :param float ambient_CO2: Air CO2 (µmol mol-1)
    :param float RH: Relative humidity (decimal fraction)

    :return: gsw (mol m-2 s-1)
    :rtype: float
    """

    Cs = ambient_CO2 - An * (parameters.K_Cs / parameters.GB)  #: CO2 concentration at organ surface (µmol mol-1 or Pa). From Prieto et al. (2012). GB in mol m-2 s-1
    m = parameters.PARAM_N['delta1'] * surfacic_nitrogen ** parameters.PARAM_N['delta2']  #: Scaling factor dependance to surfacic_nitrogen (dimensionless). This focntion is maintained
    # although I'm not conviced that it should be taken into account
    gsw = (parameters.GSMIN + m * ((Ag * RH) / Cs))  #: Stomatal conductance to water vapour (mol m-2 s-1), from Braune et al. (2009), Muller et al. (2005): using Ag rather than An.
    # Would be better with a function of VPD and with (Ci-GAMMA) instead of Cs.
    return gsw


def _calculate_Ci(ambient_CO2, An, gsw):
    """
    Calculates the internal CO2 concentration (Ci)

    :param float ambient_CO2: air CO2 (µmol mol-1)
    :param float An: net assimilation rate of CO2 (µmol m-2 s-1)
    :param float gsw: stomatal conductance to water vapour (mol m-2 s-1)

    :return: Ci (µmol mol-1)
    :rtype: float
    """
    Ci = ambient_CO2 - An * ((parameters.gsw_gs_CO2 / gsw) + (parameters.Ci_A / parameters.GB))  #: Intercellular concentration of CO2 (µmol mol-1)
    # gsw and GB in mol m-2 s-1 so that  (An * ((1.6/gs) + (1.37/parameters.GB)) is thus in µmol mol-1 as ambient_CO2

    return Ci


def _f_temperature(pname, p25, T):
    """
    Photosynthetic parameters relation to temperature

    :param str pname: name of parameter
    :param float p25: parameter value at 25 degree C
    :param float T: organ temperature (degree C)

    :return: p (parameter value at organ temperature)
    :rtype: float
    """
    Tk = T + parameters.KELVIN_DEGREE
    deltaHa = parameters.PARAM_TEMP['deltaHa'][pname]  #: Enthalpie of activation of parameter pname (kJ mol-1)
    Tref = parameters.PARAM_TEMP['Tref']

    f_activation = exp((deltaHa * (Tk - Tref)) / (parameters.R * 1E-3 * Tref * Tk))  #: Energy of activation (normalized to unity)

    if pname in ('Vc_max', 'Jmax', 'TPU'):
        deltaS = parameters.PARAM_TEMP['deltaS'][pname]  #: entropy term of parameter pname (kJ mol-1 K-1)
        deltaHd = parameters.PARAM_TEMP['deltaHd'][pname]  #: Enthalpie of deactivation of parameter pname (kJ mol-1)
        f_deactivation = (1 + exp((Tref * deltaS - deltaHd) / (Tref * parameters.R * 1E-3))) / (
                1 + exp((Tk * deltaS - deltaHd) / (Tk * parameters.R * 1E-3)))  #: Energy of deactivation (normalized to unity)
    else:
        f_deactivation = 1

    p = p25 * f_activation * f_deactivation

    return p


def _inhibition_by_NSC(NSC):
    """
    Calculates the relative diminution of Ag due to inhibition by NSC. Adapted from Azcon-Bieto 1983

    :param float NSC: Surfacic content of water-soluble carbohydrates  (µmol C m-2)

    :return: Relative diminution (dimensionless)
    :rtype: float
    """
    if NSC <= parameters.WSC_min:
        return 0
    else:
        return min(parameters.Inhibition_max * (NSC - parameters.WSC_min) / (parameters.K_Inhibition + NSC - parameters.WSC_min), 1)


def calculate_photosynthesis(PAR, surfacic_nitrogen, NSC_Retroinhibition, surfacic_NSC, Ts, Ci):
    """
    Computes photosynthesis rate following Farquhar's model with regulation by organ temperature and nitrogen content.
    In this version, most of the parameters are derived from Braune et al. (2009) on barley and Evers et al. (2010) for N dependencies.

    :param float PAR: PAR absorbed (µmol m-2 s-1)
    :param float surfacic_nitrogen: surfacic nitrogen content(g m-2) including or not structural nitrogen depending on parameter.MODEL_VERSION
    :param bool NSC_Retroinhibition: if True, Ag is inhibited by surfacic NSC (Non-Structural Carbohydrates).
    :param float surfacic_NSC: surfacic content of NSC (Non-Structural Carbohydrates) (µmol C m-2).
    :param float Ts: organ temperature (degree C)
    :param float Ci: internal CO2 (µmol mol-1), Ci = 0.7*CO2air for the first iteration

    :return: Ag (µmol m-2 s-1), An (µmol m-2 s-1), Rd (µmol m-2 s-1)
    :rtype: (float, float, float)
    """

    #: RuBisCO parameters dependance to temperature
    Kc = _f_temperature('Kc', parameters.KC25, Ts)
    Ko = _f_temperature('Ko', parameters.KO25, Ts)
    Gamma = _f_temperature('Gamma', parameters.GAMMA25, Ts)

    #: RuBisCO-limited carboxylation rate
    Sna_Vcmax25 = parameters.PARAM_N['S_surfacic_nitrogen']['Vc_max25']
    surfacic_nitrogen_min_Vcmax25 = parameters.PARAM_N['surfacic_nitrogen_min']['Vc_max25']
    Vc_max25 = Sna_Vcmax25 * (surfacic_nitrogen - surfacic_nitrogen_min_Vcmax25)  #: Relation between Vc_max25 and surfacic_nonstructural_nitrogen (µmol m-2 s-1)
    Vc_max = _f_temperature('Vc_max', Vc_max25, Ts)  #: Relation between Vc_max and temperature (µmol m-2 s-1)
    Ac = (Vc_max * (Ci - Gamma)) / (Ci + Kc * (1 + parameters.O2 / Ko))  #: Rate of assimilation under Vc_max limitation (µmol m-2 s-1)

    #: RuBP regeneration-limited carboxylation rate via electron transport
    ALPHA = parameters.PARAM_N['S_surfacic_nitrogen']['alpha'] * surfacic_nitrogen + parameters.PARAM_N['beta']  #: Relation between ALPHA and surfacic_nitrogen (mol e- mol-1 photon)
    Sna_Jmax25 = parameters.PARAM_N['S_surfacic_nitrogen']['Jmax25']
    surfacic_nitrogen_min_Jmax25 = parameters.PARAM_N['surfacic_nitrogen_min']['Jmax25']
    Jmax25 = Sna_Jmax25 * (surfacic_nitrogen - surfacic_nitrogen_min_Jmax25)  #: Relation between Jmax25 and surfacic_nitrogen (µmol m-2 s-1)
    Jmax = _f_temperature('Jmax', Jmax25, Ts)  #: Relation between Jmax and temperature (µmol m-2 s-1)

    J = ((Jmax + ALPHA * PAR) - sqrt((Jmax + ALPHA * PAR) ** parameters.J_expo - parameters.J_A * parameters.THETA * ALPHA * PAR * Jmax)) / (
            parameters.J_B * parameters.THETA)  #: Electron transport rate (Muller et al. (2005), Evers et al. (2010)) (µmol m-2 s-1)
    Aj = (J * (Ci - Gamma)) / (parameters.Aj_A * Ci + parameters.Aj_B * Gamma)  #: Rate of assimilation under RuBP regeneration limitation (µmol m-2 s-1)

    #: Triose phosphate utilisation-limited carboxylation rate --- NOT USED, calculated just for information
    Sna_TPU25 = parameters.PARAM_N['S_surfacic_nitrogen']['TPU25']
    surfacic_nitrogen_min_TPU25 = parameters.PARAM_N['surfacic_nitrogen_min']['TPU25']
    TPU25 = Sna_TPU25 * (surfacic_nitrogen - surfacic_nitrogen_min_TPU25)  #: Relation between TPU25 and surfacic_nitrogen (µmol m-2 s-1)
    TPU = _f_temperature('TPU', TPU25, Ts)  #: Relation between TPU and temperature (µmol m-2 s-1)
    Vomax = (Vc_max * Ko * Gamma) / (parameters.Vomax_A * Kc * parameters.O2)  #: Maximum rate of Vo (µmol m-2 s-1) (µmol m-2 s-1)
    Vo = (Vomax * parameters.O2) / (parameters.O2 + Ko * (1 + Ci / Kc))  #: Rate of oxygenation of RuBP (µmol m-2 s-1)
    Ap = (1 - Gamma / Ci) * (parameters.Ap_A * TPU + Vo)  #: Rate of assimilation under TPU limitation (µmol m-2 s-1).
    # I think there was a mistake in the paper of Braune et al. (2009) where they wrote Ap = (1-Gamma/Ci)*(3*TPU) + Vo
    # A more recent expression of Ap was given by S. v Caemmerer in her book (2000): AP = (3TPU * (Ci-Gamma))/(Ci-(1+3alpha)*Gamma),
    # where 0 < alpha > 1 is the fraction of glycolate carbon not returned to the chloroplast, but I couldn't find any estimation of alpha for wheat

    #: Gross assimilation rate (µmol m-2 s-1)
    if NSC_Retroinhibition:
        Ag = min(Ac, Aj) * (1 - _inhibition_by_NSC(surfacic_NSC))
    else:
        Ag = min(Ac, Aj, Ap)

    #: Mitochondrial respiration rate of organ in light Rd (processes other than photorespiration)
    Rdark25 = parameters.PARAM_N['S_surfacic_nitrogen']['Rdark25'] * (surfacic_nitrogen - parameters.PARAM_N['surfacic_nitrogen_min'][
        'Rdark25'])  #: Relation between Rdark25 (respiration in obscurity at 25 degree C) and surfacic_nitrogen (µmol m-2 s-1)
    Rdark = _f_temperature('Rdark', Rdark25, Ts)  #: Relation between Rdark and temperature (µmol m-2 s-1)
    Rd = Rdark * (parameters.Rd_A + (1 - parameters.Rd_A) * parameters.Rd_B ** (PAR / parameters.Rd_C))  # Found in Muller et al. (2005), eq. 19 (µmol m-2 s-1)

    #: Net C assimilation (µmol m-2 s-1)
    if Ag <= 0:  # Occurs when Ci is lower than Gamma or when (surfacic_nitrogen - surfacic_nitrogen_min)<0, in these cases there is no net assimilation (Farquhar, 1980; Caemmerer, 2000)
        Ag, An = 0, 0
    else:
        An = Ag - Rd

    return Ag, An, Rd


def calculate_surfacic_nitrogen(nitrates, amino_acids, proteins, Nstruct, green_area):
    """Surfacic content of nitrogen

    :param float nitrates: amount of nitrates (µmol N)
    :param float amino_acids: amount of amino_acids (µmol N)
    :param float proteins: amount of proteins (µmol N)
    :param float Nstruct: structural N (g)
    :param float green_area: green area (m-2)

    :return: Surfacic nitrogen (g m-2)
    :rtype: float
    """
    mass_N_tot = (nitrates + amino_acids + proteins) * 1E-6 * parameters.N_MOLAR_MASS + Nstruct
    return mass_N_tot / green_area


def calculate_surfacic_nonstructural_nitrogen(nitrates, amino_acids, proteins, green_area):
    """Surfacic content of non-structural nitrogen

    :param float nitrates: amount of nitrates (µmol N)
    :param float amino_acids: amount of amino_acids (µmol N)
    :param float proteins: amount of proteins (µmol N)
    :param float green_area: green area (m-2)

    :return: Surfacic non-structural nitrogen (g m-2)
    :rtype: float
    """
    mass_N_tot = (nitrates + amino_acids + proteins) * 1E-6 * parameters.N_MOLAR_MASS
    return mass_N_tot / green_area


def calculate_surfacic_photosynthetic_proteins(proteins, green_area):
    """Surfacic content of photosynthetic proteins

    :param float proteins: amount of proteins (µmol N)
    :param float green_area: green area (m-2)

    :return: Surfacic non-structural nitrogen (g m-2)
    :rtype: float
    """
    mass_N_prot = proteins * 1E-6 * parameters.N_MOLAR_MASS
    return mass_N_prot / green_area


def calculate_surfacic_nonstructural_nitrogen_Farquhar(surfacic_photosynthetic_proteins):
    """Estimate of non structural SLN used in Farquhar

    :param float surfacic_photosynthetic_proteins: surfacic proteins content (µmol N m-2)

    :return: Surfacic non-structural nitrogen (g m-2)
    :rtype: float
    """
    return surfacic_photosynthetic_proteins * parameters.Psurf_to_SLNnonstruct


def calculate_surfacic_WSC(sucrose, starch, fructan, green_area):
    """Surfacic content of water soluble carbohydrates  # TODO: rename by non structural carbohydrates because starch is not water-soluble.

    :param float sucrose: amount of sucrose (µmol C)
    :param float starch: amount of starch (µmol C)
    :param float fructan: amount of fructan (µmol C)
    :param float green_area: green area (m-2)

    :return: Surfacic content of water soluble carbohydrates  (µmol C m-2)
    :rtype: float
    """
    return (sucrose + starch + fructan) / green_area


def run(surfacic_nitrogen, NSC_Retroinhibition, surfacic_NSC, width, height, PAR, Ta, ambient_CO2, RH, Ur, organ_name, height_canopy):
    """
    Computes the photosynthesis of a photosynthetic element. The photosynthesis is computed by using the biochemical FCB model (Farquhar et al., 1980) coupled to the semiempirical
    BWB model of stomatal conductance (Ball, 1987).

    :param float surfacic_nitrogen: surfacic nitrogen content of organs (g m-2), could include total N or proteic N only depending on parameter.SurfacicProteins
           Properly speaking, photosynthesis should be related to proteins (RubisCO), but parameters of most Farquhar models are calibrated on total N measurements (DUMAS method).
           We use only non-structural nitrogen to overcome issues in the case of extrem scenarios (high SLN for thick leaves under low nitrogen conditions).
           If None, surfacic_nitrogen = :attr:`NA_0`
    :param bool NSC_Retroinhibition: if True, Ag is inhibited by surfacic NSC (Non-Structural Carbohydrates).
    :param float surfacic_NSC: surfacic content of NSC (Non-Structural Carbohydrates) (µmol C m-2).
    :param float width: width of the organ (or diameter for stem organ) (m),
           characteristic dimension to be considered for heat transfer through forced convection (by wind).
    :param float height: height of the organ from soil (m)
    :param float PAR: absorbed PAR (µmol m-2 s-1)
    :param float Ta: air temperature (°C)
    :param float ambient_CO2: air CO2 (µmol mol-1)
    :param float RH: relative humidity (decimal fraction)
    :param float Ur: wind at the reference height (zr) (m s-1), e.g. top of the canopy + 2m
           (in the case of wheat, Ur can be approximated as the wind speed at 2m from soil)
    :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from cylindric organs)
    :param float height_canopy: total canopy height (m)

    :return: Ag (µmol m-2 s-1), An (µmol m-2 s-1), Rd (µmol m-2 s-1),
        Tr (mmol m-2 s-1), Ts (°C) and  gsw (mol m-2 s-1)
    :rtype: (float, float, float, float, float, float)
    """

    if surfacic_nitrogen is None:
        surfacic_nitrogen = parameters.NA_0

    # Iterations to find organ temperature and Ci #
    Ci, Ts = parameters.Ci_init_ratio * ambient_CO2, Ta  # Initial values
    count = 0

    while True:
        prec_Ci, prec_Ts = Ci, Ts
        Ag, An, Rd = calculate_photosynthesis(PAR, surfacic_nitrogen, NSC_Retroinhibition, surfacic_NSC, Ts, Ci)
        # Stomatal conductance to water
        gsw = _stomatal_conductance(Ag, An, surfacic_nitrogen, ambient_CO2, RH)

        # New value of Ci
        Ci = _calculate_Ci(ambient_CO2, An, gsw)

        # New value of Ts
        Ts, Tr = _organ_temperature(width, height, height_canopy, Ur, PAR, gsw, Ta, Ts, RH, organ_name)
        count += 1

        if count >= 30:  # TODO: test a faire? Semble prendre du tps de calcul
            if abs((Ci - prec_Ci) / prec_Ci) >= parameters.DELTA_CONVERGENCE:
                print('{}, Ci cannot converge, prec_Ci= {}, Ci= {}'.format(organ_name, prec_Ci, Ci))
            if prec_Ts != 0 and abs((Ts - prec_Ts) / prec_Ts) >= parameters.DELTA_CONVERGENCE:
                print('{}, Ts cannot converge, prec_Ts= {}, Ts= {}'.format(organ_name, prec_Ts, Ts))
            break
        if abs((Ci - prec_Ci) / prec_Ci) < parameters.DELTA_CONVERGENCE and ((prec_Ts == 0 and (Ts - prec_Ts) == 0) or abs((Ts - prec_Ts) / prec_Ts) < parameters.DELTA_CONVERGENCE):
            break

    #: Conversion of Tr from mm s-1 to mmol m-2 s-1 (more suitable for further use of Tr)
    Tr = (Tr * 1E6) / parameters.MM_WATER  # Using 1 mm = 1kg m-2
    #: Decrease efficency of non-lamina organs
    if organ_name != 'blade':
        Ag = Ag * parameters.EFFICENCY_STEM
    return Ag, An, Rd, Tr, Ts, gsw
