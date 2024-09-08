#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@illinois.edu>
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Andrew Koehler <koehler@mines.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from chaospy import distributions as shape
from qsdsan.utils import auom, DictAttrSetter

__all__ = ('create_model',)

# kg/m3
biocrude_density = 980
# kg/m3, this is for sludge with a moisture content higher than 80%,
# google 'Design of wastewater treatment sludge thickeners Iowa State University'
sludge_density = 1000
# kg/L
water_density = 1

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_MMgal_to_L = auom('gal').conversion_factor('L')*1000000
_oil_barrel_to_L = auom('oil_barrel').conversion_factor('L')
_mile_to_km = auom('mile').conversion_factor('km')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_ft3 = auom('m3').conversion_factor('ft3')
_ton_to_kg = auom('ton').conversion_factor('kg')
_MJ_to_MMBTU = auom('MJ').conversion_factor('MMBTU')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index), [2]
GDPCTPI = {2007: 86.352,
           2008: 87.977,
           2009: 88.557,
           2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

def create_model(system=None,
                 exclude_sludge_compositions=False,
                 include_HTL_yield_as_metrics=True,
                 include_other_metrics=False,
                 include_check=True):
    '''
    Create a model based on the given system
    (or create the system based on the given configuration).
    
    Parameters
    ----------
    system : obj or str
        Can be either a :class:`System` objective for which the model will be created,
        or one of the allowed configurations ("baseline", "no_P", "PSA").
    '''
    
    if system == None:
        raise RuntimeError('must create a system first')
    else:
        sys = system
    
    flowsheet = sys.flowsheet
    cmps = qs.get_components()
    unit = flowsheet.unit
    stream = flowsheet.stream
    model = qs.Model(sys)
    param = model.parameter
    metric = model.metric
    
    WWTP = unit.WWTP
    P1 = unit.P1
    H1 = unit.H1
    HTL = unit.HTL
    H2SO4_Tank = unit.H2SO4_Tank
    SP1 = unit.SP1
    StruPre = unit.StruPre
    CHG = unit.CHG
    F1 = unit.F1
    MemDis = unit.MemDis
    P2 = unit.P2
    HT = unit.HT
    H2 = unit.H2
    F2 = unit.F2
    H3 = unit.H3
    D1 = unit.D1
    D2 = unit.D2
    D3 = unit.D3
    P3 = unit.P3
    HC = unit.HC
    H4 = unit.H4
    F3 = unit.F3
    D4 = unit.D4
    H5 = unit.H5
    H6 = unit.H6
    GasolineTank = unit.GasolineTank
    DieselTank = unit.DieselTank
    HXN = unit.HXN
    CHP = unit.CHP
    CT = unit.CT
    CWP = unit.CWP
    HT_HC_units = (P2, HT, H2, F2, H3,
                   D1, D2, D3, P3, HC,
                   H4, F3, D4, H5, H6,
                   GasolineTank,
                   DieselTank)
    
    raw_wastewater = stream.raw_wastewater
    NaOH_HTL = stream.NaOH_HTL
    HCl = stream.HCl
    H2SO4 = stream.H2SO4
    MgCl2 = stream.MgCl2
    NH4Cl = stream.NH4Cl
    MgO = stream.MgO
    struvite = stream.struvite
    CHG_catalyst_in = stream.virgin_CHG_catalyst
    NaOH_MemDis = stream.NaOH_MemDis
    Membrane_in = stream.Membrane_in
    ammonium_sulfate = stream.ammonium_sulfate
    hydrogen_gas = stream.hydrogen_gas
    CoMo_alumina_HT = stream.CoMo_alumina_HT
    CoMo_alumina_HC = stream.CoMo_alumina_HC
    gasoline = stream.gasoline
    diesel = stream.diesel
    natural_gas = stream.natural_gas
    solid_ash = stream.solid_ash
    cooling_tower_chemicals = stream.cooling_tower_chemicals
    
    feedstock = sys.flowsheet.WWTP.feedstock
    high_IRR = sys.flowsheet.WWTP.high_IRR
    
    tea = sys.TEA
    lca = sys.LCA
    
    # =========================================================================
    # WWTP
    # =========================================================================
    
    if not exclude_sludge_compositions:
        # TODO: update values for sludge
        if feedstock == 'sludge':
            
            dist = shape.Uniform(0.6,0.8)
            @param(name='sludge_moisture',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=0.7,
                   distribution=dist)
            def set_WWTP_sludge_moisture(i):
                WWTP.moisture=i
            
            dist = shape.Triangle(0.174,0.257,0.414)
            @param(name='sludge_dw_ash',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=0.257,
                   distribution=dist)
            def set_sludge_dw_ash(i):
                WWTP.dw_ash=i
            
            dist = shape.Triangle(0.08,0.204,0.308)
            @param(name='sludge_afdw_lipid',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=0.204,
                   distribution=dist)
            def set_sludge_afdw_lipid(i):
                WWTP.afdw_lipid=i
            
            dist = shape.Triangle(0.38,0.463,0.51)
            @param(name='sludge_afdw_protein',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=0.463,
                   distribution=dist)
            def set_sludge_afdw_protein(i):
                WWTP.afdw_protein=i
                
            dist = shape.Triangle(0.1944,0.3927,0.5556)
            @param(name='sludge_N_2_P',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=0.3927,
                   distribution=dist)
            def set_sludge_N_2_P(i):
                WWTP.N_2_P=i
        
        # TODO: update values for biosolid
        if feedstock == 'biosolid':
           
           dist = shape.Uniform(0.6,0.8)
           @param(name='biosolid_moisture',
                  element=WWTP,
                  kind='coupled',
                  units='-',
                  baseline=0.7,
                  distribution=dist)
           def set_WWTP_biosolid_moisture(i):
               WWTP.moisture=i
           
           dist = shape.Triangle(0.174,0.257,0.414)
           @param(name='sludge_dw_ash',
                  element=WWTP,
                  kind='coupled',
                  units='-',
                  baseline=0.257,
                  distribution=dist)
           def set_biosolid_dw_ash(i):
               WWTP.dw_ash=i
           
           dist = shape.Triangle(0.08,0.204,0.308)
           @param(name='sludge_afdw_lipid',
                  element=WWTP,
                  kind='coupled',
                  units='-',
                  baseline=0.204,
                  distribution=dist)
           def set_biosolid_afdw_lipid(i):
               WWTP.afdw_lipid=i
           
           dist = shape.Triangle(0.20,0.25,0.30)
           @param(name='sludge_afdw_protein',
                  element=WWTP,
                  kind='coupled',
                  units='-',
                  baseline=0.463,
                  distribution=dist)
           def set_biosolid_afdw_protein(i):
               WWTP.afdw_protein=i
               
           dist = shape.Triangle(0.1944,0.3927,0.5556)
           @param(name='biosolid_N_2_P',
                  element=WWTP,
                  kind='coupled',
                  units='-',
                  baseline=0.3927,
                  distribution=dist)
           def set_biosolid_N_2_P(i):
               WWTP.N_2_P=i

    dist = shape.Uniform(0.675,0.825)
    @param(name='lipid_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.75,
           distribution=dist)
    def set_lipid_2_C(i):
        WWTP.lipid_2_C=i
    
    dist = shape.Uniform(0.4905,0.5995)
    @param(name='protein_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.545,
           distribution=dist)
    def set_protein_2_C(i):
        WWTP.protein_2_C=i
    
    dist = shape.Uniform(0.36,0.44)
    @param(name='carbo_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.4,
           distribution=dist)
    def set_carbo_2_C(i):
        WWTP.carbo_2_C=i
    
    dist = shape.Uniform(0.1125,0.1375)
    @param(name='lipid_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.125,
           distribution=dist)
    def set_lipid_2_H(i):
        WWTP.lipid_2_H=i
    
    dist = shape.Uniform(0.0614,0.075)
    @param(name='protein_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.0682,
           distribution=dist)
    def set_protein_2_H(i):
        WWTP.protein_2_H=i
    
    dist = shape.Uniform(0.06,0.0733)
    @param(name='carbo_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.0667,
           distribution=dist)
    def set_carbo_2_H(i):
        WWTP.carbo_2_H=i

    dist = shape.Uniform(0.1432,0.1750)
    @param(name='protein_2_N',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.159,
           distribution=dist)
    def set_protein_2_N(i):
        WWTP.protein_2_N=i

    dist = shape.Triangle(7392,7920,8448)
    @param(name='operation_hour',
           element=WWTP,
           kind='coupled',
           units='hr/yr',
           baseline=7920,
           distribution=dist)
    def set_operation_hour(i):
        WWTP.operation_hours=sys.operating_hours=i
    
    # =========================================================================
    # HTL
    # =========================================================================
    dist = shape.Uniform(0.0170348,0.0227131)
    @param(name='enforced heating transfer coefficient',
           element=H1,
           kind='coupled',
           units='kW/m2/K',
           baseline=0.0198739,
           distribution=dist)
    def set_U(i):
        H1.U=i
    
    # TODO: if we add heat loss, consider add that as an uncertainty parameter here
    
    dist = shape.Triangle(0.692,0.846,1)
    @param(name='lipid_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.846,
           distribution=dist)
    def set_lipid_2_biocrude(i):
        HTL.lipid_2_biocrude=i
    
    dist = shape.Normal(0.445,0.030)
    @param(name='protein_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.445,
           distribution=dist)
    def set_protein_2_biocrude(i):
        HTL.protein_2_biocrude=i
    
    dist = shape.Normal(0.205,0.050)
    @param(name='carbo_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.205,
           distribution=dist)
    def set_carbo_2_biocrude(i):
        HTL.carbo_2_biocrude=i
    
    dist = shape.Normal(0.074,0.020)
    @param(name='protein_2_gas',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.074,
           distribution=dist)
    def set_protein_2_gas(i):
        HTL.protein_2_gas=i
    
    dist = shape.Normal(0.418,0.030)
    @param(name='carbo_2_gas',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.418,
           distribution=dist)
    def set_carbo_2_gas(i):
        HTL.carbo_2_gas=i
        
    dist = shape.Normal(-8.370,0.939)
    @param(name='biocrude_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=-8.370,
           distribution=dist)
    def set_biocrude_C_slope(i):
        HTL.biocrude_C_slope=i
        
    dist = shape.Normal(68.55,0.367)
    @param(name='biocrude_C_intercept',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=68.55,
           distribution=dist)
    def set_biocrude_C_intercept(i):
        HTL.biocrude_C_intercept=i
        
    dist = shape.Normal(0.133,0.005)
    @param(name='biocrude_N_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.133,
           distribution=dist)
    def set_biocrude_N_slope(i):
        HTL.biocrude_N_slope=i
        
    dist = shape.Normal(-2.610,0.352)
    @param(name='biocrude_H_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=-2.610,
           distribution=dist)
    def set_biocrude_H_slope(i):
        HTL.biocrude_H_slope=i
    
    dist = shape.Normal(8.200,0.138)
    @param(name='biocrude_H_intercept',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=8.200,
           distribution=dist)
    def set_biocrude_H_intercept(i):
        HTL.biocrude_H_intercept=i
        
    dist = shape.Normal(478,18.878)
    @param(name='HTLaqueous_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=478,
           distribution=dist)
    def set_HTLaqueous_C_slope(i):
        HTL.HTLaqueous_C_slope=i
        
    dist = shape.Triangle(0.715,0.764,0.813)
    @param(name='TOC_TC',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.764,
           distribution=dist)
    def set_TOC_TC(i):
        HTL.TOC_TC=i
    
    dist = shape.Normal(1.750,0.122)
    @param(name='hydrochar_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=1.750,
           distribution=dist)
    def set_hydrochar_C_slope(i):
        HTL.hydrochar_C_slope=i
    
    dist = shape.Triangle(0.035,0.063,0.102)
    @param(name='biocrude_moisture_content',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.063,
           distribution=dist)
    def set_biocrude_moisture_content(i):
        HTL.biocrude_moisture_content=i
    
    dist = shape.Uniform(0.84,0.88)
    @param(name='hydrochar_P_recovery_ratio',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.86,
           distribution=dist)
    def set_hydrochar_P_recovery_ratio(i):
        HTL.hydrochar_P_recovery_ratio=i
    
    # =========================================================================
    # AcidEx
    # =========================================================================
    if AcidEx:= unit.search('AcidEx'):
        dist = shape.Uniform(4,10)
        @param(name='acid_vol',
               element=AcidEx,
               kind='coupled',
               units='-',
               baseline=7,
               distribution=dist)
        def set_acid_vol(i):
            AcidEx.acid_vol=i
        
        dist = shape.Uniform(0.7,0.9)
        @param(name='P_acid_recovery_ratio',
               element=AcidEx,
               kind='coupled',
               units='-',
               baseline=0.8,
               distribution=dist)
        def set_P_recovery_ratio(i):
            AcidEx.P_acid_recovery_ratio=i
    
    # =========================================================================
    # StruPre
    # =========================================================================
    dist = shape.Uniform(8.5,9.5)
    @param(name='target_pH',
           element=StruPre,
           kind='coupled',
           units='-',
           baseline=9,
           distribution=dist)
    def set_StruPre_target_pH(i):
        StruPre.target_pH=i
        
    dist = shape.Triangle(0.7,0.828,0.95)
    @param(name='P_pre_recovery_ratio',
           element=StruPre,
           kind='coupled',
           units='-',
           baseline=0.828,
           distribution=dist)
    def set_P_pre_recovery_ratio(i):
        StruPre.P_pre_recovery_ratio=i
    
    # =========================================================================
    # CHG
    # =========================================================================
    dist = shape.Triangle(2.86,3.562,3.99)
    @param(name='WHSV',
           element=CHG,
           kind='coupled',
           units='kg/hr/kg',
           baseline=3.562,
           distribution=dist)
    def set_CHG_WHSV(i):
        CHG.WHSV=i
    
    dist = shape.Triangle(3960,7920,15840)
    @param(name='catalyst_lifetime',
           element=CHG,
           kind='coupled',
           units='hr',
           baseline=7920,
           distribution=dist)
    def set_CHG_catalyst_lifetime(i):
        CHG.catalyst_lifetime=i
    
    dist = shape.Triangle(0.1893,0.5981,0.7798)
    @param(name='gas_C_2_total_C',
           element=CHG,
           kind='coupled',
           units='-',
           baseline=0.5981,
           distribution=dist)
    def set_gas_C_2_total_C(i):
        CHG.gas_C_2_total_C=i
    
    # =========================================================================
    # MemDis
    # =========================================================================
    dist = shape.Uniform(7.91,8.41)
    @param(name='influent_pH',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=8.16,
           distribution=dist)
    def set_influent_pH(i):
        MemDis.influent_pH=i
    
    dist = shape.Uniform(10,11.8)
    @param(name='target_pH',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=10,
           distribution=dist)
    def set_MemDis_target_pH(i):
        MemDis.target_pH=i
    
    dist = shape.Uniform(0.00075,0.000917)
    @param(name='m2_2_m3',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=1/1200,
           distribution=dist)
    def set_m2_2_m3(i):
        MemDis.m2_2_m3=i
    
    dist = shape.Uniform(0.0000205,0.0000251)
    @param(name='Dm',
           element=MemDis,
           kind='coupled',
           units='m2/s',
           baseline=0.0000228,
           distribution=dist)
    def set_Dm(i):
        MemDis.Dm=i
    
    dist = shape.Uniform(0.81,0.99)
    @param(name='porosity',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=0.9,
           distribution=dist)
    def set_porosity(i):
        MemDis.porosity=i
    
    dist = shape.Uniform(0.000063,0.000077)
    @param(name='thickness',
           element=MemDis,
           kind='coupled',
           units='m',
           baseline=0.00007,
           distribution=dist)
    def set_thickness(i):
        MemDis.thickness=i
    
    dist = shape.Uniform(1.08,1.32)
    @param(name='tortuosity',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=1.2,
           distribution=dist)
    def set_tortuosity(i):
        MemDis.tortuosity=i
    
    dist = shape.Uniform(0.0000158,0.0000193)
    @param(name='Ka',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=0.0000175,
           distribution=dist)
    def set_Ka(i):
        MemDis.Ka=i
    
    dist = shape.Uniform(5.409,6.611)
    @param(name='capacity',
           element=MemDis,
           kind='coupled',
           units='-',
           baseline=6.01,
           distribution=dist)
    def set_capacity(i):
        MemDis.capacity=i
    
    # =========================================================================
    # HT
    # =========================================================================
    dist = shape.Uniform(0.5625,0.6875)
    @param(name='WHSV',
           element=HT,
           kind='coupled',
           units='kg/hr/kg',
           baseline=0.625,
           distribution=dist)
    def set_HT_WHSV(i):
        HT.WHSV=i
    
    dist = shape.Triangle(7920,15840,39600)
    @param(name='catalyst_lifetime',
           element=HT,
           kind='coupled',
           units='hr',
           baseline=15840,
           distribution=dist)
    def set_HT_catalyst_lifetime(i):
        HT.catalyst_lifetime=i
    
    dist = shape.Uniform(0.0414,0.0506)
    @param(name='hydrogen_rxned_to_biocrude',
           element=HT,
           kind='coupled',
           units='-',
           baseline=0.046,
           distribution=dist)
    def set_HT_hydrogen_rxned_to_biocrude(i):
        HT.hydrogen_rxned_to_biocrude=i
    
    if 'PSA' in flowsheet.ID:
        dist = shape.Uniform(0.8,0.9)
        @param(name='PSA_efficiency',
               element=HT,
               kind='coupled',
               units='-',
               baseline=0.9,
               distribution=dist)
        def set_PSA_efficiency(i):
            HT.PSA_efficiency=i
    
    dist = shape.Uniform(2.4,3.6)
    @param(name='hydrogen_excess',
           element=HT,
           kind='coupled',
           units='-',
           baseline=3,
           distribution=dist)
    def set_HT_hydrogen_excess(i):
        HT.hydrogen_excess=i
    
    dist = shape.Uniform(0.7875,0.9625)
    @param(name='hydrocarbon_ratio',
           element=HT,
           kind='coupled',
           units='-',
           baseline=0.875,
           distribution=dist)
    def set_HT_hydrocarbon_ratio(i):
        HT.hydrocarbon_ratio=i
    
    # =========================================================================
    # HC
    # =========================================================================
    dist = shape.Uniform(0.5625,0.6875)
    @param(name='WHSV',
           element=HC,
           kind='coupled',
           units='kg/hr/kg',
           baseline=0.625,
           distribution=dist)
    def set_HC_WHSV(i):
        HC.WHSV=i
    
    dist = shape.Uniform(35640,43560)
    @param(name='catalyst_lifetime',
           element=HC,
           kind='coupled',
           units='hr',
           baseline=39600,
           distribution=dist)
    def set_HC_catalyst_lifetime(i):
        HC.catalyst_lifetime=i
    
    dist = shape.Uniform(0.010125,0.012375)
    @param(name='hydrogen_rxned_to_heavy_oil',
           element=HC,
           kind='coupled',
           units='-',
           baseline=0.01125,
           distribution=dist)
    def set_HC_hydrogen_rxned_to_heavy_oil(i):
        HC.hydrogen_rxned_to_heavy_oil=i
    
    dist = shape.Uniform(4.4448,6.6672)
    @param(name='hydrogen_excess',
           element=HC,
           kind='coupled',
           units='-',
           baseline=5.556,
           distribution=dist)
    def set_HC_hydrogen_excess(i):
        HC.hydrogen_excess=i
    
    dist = shape.Uniform(0.9,1)
    @param(name='hydrocarbon_ratio',
           element=HC,
           kind='coupled',
           units='-',
           baseline=1,
           distribution=dist)
    def set_HC_hydrocarbon_ratio(i):
        HC.hydrocarbon_ratio=i
    
    # =========================================================================
    # TEA
    # =========================================================================
    dist = shape.Triangle(0.6,1,1.4)
    @param(name='HTL_TIC_factor',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1,
           distribution=dist)
    def set_HTL_TIC_factor(i):
        HTL.TIC_factor=i
        
    dist = shape.Triangle(0.6,1,1.4)
    @param(name='CHG_TIC_factor',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1,
           distribution=dist)
    def set_CHG_TIC_factor(i):
        CHG.TIC_factor=i
    
    dist = shape.Triangle(0.6,1,1.4)
    @param(name='HT_TIC_factor',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1,
           distribution=dist)
    def set_HT_TIC_factor(i):
        HT.TIC_factor=i
    
    dist = shape.Uniform(980,1470)
    @param(name='unit_TIC',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1225,
           distribution=dist)
    def set_unit_TIC(i):
        CHP.unit_TIC=i
    
    if high_IRR:
        dist = shape.Triangle(0.05,0.1,0.15)
        @param(name='IRR',
               element='TEA',
               kind='isolated',
               units='-',
               baseline=0.1,
               distribution=dist)
        def set_IRR(i):
            tea.IRR=i
            
    else:
        dist = shape.Triangle(0,0.03,0.05)
        @param(name='IRR',
               element='TEA',
               kind='isolated',
               units='-',
               baseline=0.03,
               distribution=dist)
        def set_IRR(i):
            tea.IRR=i
    
    H2SO4_price = (0.043*1+0.0002*(93/5-1))/(93/5)/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(H2SO4_price*0.9,H2SO4_price*1.1)
    @param(name='5% H2SO4 price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=H2SO4_price,
           distribution=dist)
    def set_H2SO4_price(i):
        H2SO4.price=i
        
    NaOH_price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(NaOH_price*0.9,NaOH_price*1.1)
    @param(name='NaOH price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=NaOH_price,
           distribution=dist)
    def set_NaOH_HTL_price(i):
        NaOH_HTL.price=NaOH_MemDis.price=i
    
    magnesium_chloride_price_min = 525/_ton_to_kg
    magnesium_chloride_price_max = 580/_ton_to_kg
    magnesium_chloride_price_ave = 535/_ton_to_kg
    dist = shape.Uniform(magnesium_chloride_price_min,magnesium_chloride_price_max)
    @param(name='MgCl2 price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=magnesium_chloride_price_ave,
           distribution=dist)
    def set_MgCl2_price(i):
        MgCl2.price=i
    
    ammonium_chloride_price_min = 0.37
    ammonium_chloride_price_max = 0.53
    ammonium_chloride_price_ave = 0.4475
    dist = shape.Triangle(ammonium_chloride_price_min,ammonium_chloride_price_ave,ammonium_chloride_price_max)
    @param(name='NH4Cl price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=ammonium_chloride_price_ave,
           distribution=dist)
    def set_NH4Cl_price(i):
        NH4Cl.price=i
    
    dist = shape.Uniform(0.1,0.3)
    @param(name='MgO price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=0.2,
           distribution=dist)
    def set_MgO_price(i):
        MgO.price=i
    
    struvite_price_min = 0.19/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    struvite_price_max = 0.55/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    struvite_price_ave = 0.30/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Triangle(struvite_price_min,struvite_price_ave,struvite_price_max)
    @param(name='struvite price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=struvite_price_ave,
           distribution=dist)
    def set_struvite_price(i):
        struvite.price=i        
    
    hydrogen_gas_price = 0.7306/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(hydrogen_gas_price*0.9,hydrogen_gas_price*1.1)
    @param(name='hydrogen gas price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=hydrogen_gas_price,
           distribution=dist)
    def set_hydrogen_gas_price(i):
        hydrogen_gas.price=i
    
    ammonium_sulfate_price_min = 0.22
    ammonium_sulfate_price_max = 0.32
    ammonium_sulfate_price_ave = 0.2825
    dist = shape.Triangle(ammonium_sulfate_price_min,ammonium_sulfate_price_ave,ammonium_sulfate_price_max)
    @param(name='ammonium sulfate price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=ammonium_sulfate_price_ave,
           distribution=dist)
    def set_ammonium_sulfate_price(i):
        ammonium_sulfate.price=i
    
    membrane_price = 90/GDPCTPI[2008]*GDPCTPI[2022]
    dist = shape.Uniform(membrane_price*0.9,membrane_price*1.1)
    @param(name='membrane price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=membrane_price,
           distribution=dist)
    def set_membrane_price(i):
        MemDis.membrane_price=i
        Membrane_in.price=i
    
    CHG_catalyst_price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2022]
    dist = shape.Triangle(CHG_catalyst_price*0.5,CHG_catalyst_price,CHG_catalyst_price*2)
    @param(name='CHG catalyst price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=CHG_catalyst_price,
           distribution=dist)
    def set_catalyst_price(i):
        CHG_catalyst_in.price=i
    
    natural_gas_price_min = 6.18/1000*_m3_to_ft3/0.657
    natural_gas_price_max = 9.58/1000*_m3_to_ft3/0.657
    natural_gas_price_ave = 7.66/1000*_m3_to_ft3/0.657
    dist = shape.Triangle(natural_gas_price_min,natural_gas_price_ave,natural_gas_price_max)
    @param(name='natural gas price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=natural_gas_price_ave,
           distribution=dist)
    def set_CH4_price(i):
        natural_gas.price=i
    
    HT_HC_catalyst_price = 15.5/_lb_to_kg/GDPCTPI[2007]*GDPCTPI[2022]
    dist = shape.Uniform(0.9*HT_HC_catalyst_price,1.1*HT_HC_catalyst_price)
    @param(name='HT & HC catalyst price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=HT_HC_catalyst_price,
           distribution=dist)
    def set_HT_HC_catalyst_price(i):
        CoMo_alumina_HT.price=CoMo_alumina_HC.price=i
    
    gasoline_price_min = 0.97
    gasoline_price_max = 1.43
    gasoline_price_ave = 1.177
    dist = shape.Triangle(gasoline_price_min,gasoline_price_ave,gasoline_price_max)
    @param(name='gasoline price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=gasoline_price_ave,
           distribution=dist)
    def set_gasoline_price(i):
        gasoline.price=i
    
    diesel_price_min = 0.99
    diesel_price_max = 1.53
    diesel_price_ave = 1.331
    dist = shape.Triangle(diesel_price_min,diesel_price_ave,diesel_price_max)
    @param(name='diesel price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=diesel_price_ave,
           distribution=dist)
    def set_diesel_price(i):
        diesel.price=i
    
    ash_disposal_price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(ash_disposal_price*1.1,ash_disposal_price*0.9)
    @param(name='ash disposal price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=ash_disposal_price,
           distribution=dist)
    def set_ash_disposal_price(i):
        solid_ash.price=i
    
    cooling_chemicals_price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(cooling_chemicals_price*0.9,cooling_chemicals_price*1.1)
    @param(name='cooling tower chemicals price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=cooling_chemicals_price,
           distribution=dist)
    def set_cooling_tower_chemicals_price(i):
        cooling_tower_chemicals.price=i
    
    electricity_price_min = qs.PowerUtility.price*0.9
    electricity_price_max = qs.PowerUtility.price*1.1
    electricity_price_ave = qs.PowerUtility.price
    dist = shape.Uniform(electricity_price_min,electricity_price_max)
    @param(name='electricity price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=electricity_price_ave,
           distribution=dist)
    def set_electricity_price(i):
        qs.PowerUtility.price=i
    
    # =========================================================================
    # LCA (unifrom ± 10%)
    # =========================================================================
    CF = 'GlobalWarming'
    # don't get joint distribution for multiple times, since the baselines for LCA will change.
    qs.ImpactItem.get_all_items().pop('raw_wastewater')
    for item in qs.ImpactItem.get_all_items().keys():
        if qs.ImpactItem.get_item(item).CFs:
            if qs.ImpactItem.get_item(item).CFs[CF] != 0:
                abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
                abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
                dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
                @param(name=f'{item}_{CF}',
                       setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', CF),
                       element='LCA',
                       kind='isolated',
                       units=qs.ImpactIndicator.get_indicator(CF).unit,
                       baseline=qs.ImpactItem.get_item(item).CFs[CF],
                       distribution=dist)
                def set_LCA(i):
                    qs.ImpactItem.get_item(item).CFs[CF]=i
    
    # =========================================================================
    # metrics
    # =========================================================================
    if include_other_metrics:
    
        # element metrics
        
        @metric(name='C_afdw',units='%',element='Sankey')
        def get_C_afdw():
            return WWTP.C*24/1000/(raw_wastewater.F_mass/157262.48454459725)/(1-WWTP.dw_ash)
        
        @metric(name='N_afdw',units='%',element='Sankey')
        def get_N_afdw():
            return WWTP.N*24/1000/(raw_wastewater.F_mass/157262.48454459725)/(1-WWTP.dw_ash)
        
        @metric(name='P_afdw',units='%',element='Sankey')
        def get_P_afdw():
            return WWTP.P*24/1000/(raw_wastewater.F_mass/157262.48454459725)/(1-WWTP.dw_ash)
        
        @metric(name='sludge_C',units='kg/hr',element='Sankey')
        def get_sludge_C():
            return WWTP.C
        
        @metric(name='sludge_N',units='kg/hr',element='Sankey')
        def get_sludge_N():
            return WWTP.N
        
        @metric(name='sludge_P',units='kg/hr',element='Sankey')
        def get_sludge_P():
            return WWTP.P
        
        @metric(name='HTLaqueous_C',units='kg/hr',element='Sankey')
        def get_HTLaqueous_C():
            return HTL.HTLaqueous_C
        
        @metric(name='HTLaqueous_N',units='kg/hr',element='Sankey')
        def get_HTLaqueous_N():
            return HTL.HTLaqueous_N
        
        @metric(name='HTLaqueous_P',units='kg/hr',element='Sankey')
        def get_HTLaqueous_P():
            return HTL.HTLaqueous_P
        
        @metric(name='biocrude_C',units='kg/hr',element='Sankey')
        def get_biocrude_C():
            return HTL.biocrude_C
        
        @metric(name='biocrude_N',units='kg/hr',element='Sankey')
        def get_biocrude_N():
            return HTL.biocrude_N
        
        @metric(name='offgas_C',units='kg/hr',element='Sankey')
        def get_offgas_C():
            return HTL.offgas_C
        
        @metric(name='hydrochar_C',units='kg/hr',element='Sankey')
        def get_hydrochar_C():
            return HTL.hydrochar_C
        
        @metric(name='hydrochar_P',units='kg/hr',element='Sankey')
        def get_hydrochar_P():
            return HTL.hydrochar_P
        
        @metric(name='HT_gasoline_C',units='kg/hr',element='Sankey')
        def get_HT_gasoline_C():
            return sum(D2.outs[0].mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HT_gasoline_N',units='kg/hr',element='Sankey')
        def get_HT_gasoline_N():
            return sum(D2.outs[0].mass*[cmp.i_N for cmp in cmps])
        
        @metric(name='HT_diesel_C',units='kg/hr',element='Sankey')
        def get_HT_diesel_C():
            return sum(D3.outs[0].mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HT_heavy_oil_C',units='kg/hr',element='Sankey')
        def get_HT_heavy_oil_C():
            return sum(D3.outs[1].mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HT_gas_C',units='kg/hr',element='Sankey')
        def get_HT_gas_C():
            mass = F2.outs[0].mass + D1.outs[0].mass
            return sum(mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HT_gas_N',units='kg/hr',element='Sankey')
        def get_HT_gas_N():
            mass = F2.outs[0].mass + D1.outs[0].mass
            return sum(mass*[cmp.i_N for cmp in cmps])
    
        @metric(name='HT_ww_C',units='kg/hr',element='Sankey')
        def get_HT_ww_C():
            mass = sum((i.outs[0].mass for i in (D1, D2, D3, F2)), D3.outs[1].mass)
            return HTL.biocrude_C - sum(mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HT_ww_N',units='kg/hr',element='Sankey')
        def get_HT_ww_N():
            mass = sum((i.outs[0].mass for i in (D1, D2, D3, F2)), D3.outs[1].mass)
            return HTL.biocrude_N - sum(mass*[cmp.i_N for cmp in cmps])
        
        @metric(name='HC_gasoline_C',units='kg/hr',element='Sankey')
        def get_HC_gasoline_C():
            return sum(D4.outs[0].mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HC_diesel_C',units='kg/hr',element='Sankey')
        def get_HC_diesel_C():
            return sum(D4.outs[1].mass*[cmp.i_C for cmp in cmps])
        
        @metric(name='HC_gas_C',units='kg/hr',element='Sankey')
        def get_HC_gas_C():
            return sum(F3.outs[0].mass*[cmp.i_C for cmp in cmps])
        
        if AcidEx:
            @metric(name='extracted_P',units='kg/hr',element='Sankey')
            def get_extracted_P():
                return AcidEx.outs[1].imass['P']
            
            @metric(name='residual_P',units='kg/hr',element='Sankey')
            def get_residual_P():
                return HTL.hydrochar_P-AcidEx.outs[1].imass['P']
        
        @metric(name='residual_C',units='kg/hr',element='Sankey')
        def get_residual_C():
            return HTL.hydrochar_C
        
        @metric(name='struvite_N',units='kg/hr',element='Sankey')
        def get_struvite_N():
            return StruPre.struvite_N
        
        @metric(name='struvite_P',units='kg/hr',element='Sankey')
        def get_struvite_P():
            return StruPre.struvite_P
        
        @metric(name='CHG_feed_C',units='kg/hr',element='Sankey')
        def get_CHG_feed_C():
            return StruPre.outs[1].imass['C']
        
        @metric(name='CHG_feed_N',units='kg/hr',element='Sankey')
        def get_CHG_feed_N():
            return StruPre.outs[1].imass['N']
        
        @metric(name='CHG_feed_P',units='kg/hr',element='Sankey')
        def get_CHG_feed_P():
            return StruPre.outs[1].imass['P']
        
        @metric(name='CHG_out_C',units='kg/hr',element='Sankey')
        def get_CHG_out_C():
            return CHG.CHGout_C
        
        @metric(name='CHG_out_N',units='kg/hr',element='Sankey')
        def get_CHG_out_N():
            return CHG.CHGout_N
        
        @metric(name='CHG_out_P',units='kg/hr',element='Sankey')
        def get_CHG_out_P():
            return CHG.CHGout_P
        
        @metric(name='CHG_gas_C',units='kg/hr',element='Sankey')
        def get_CHG_gas_C():
            return sum(F1.outs[0].mass*[cmp.i_C for cmp in cmps])
    
        @metric(name='ammoniumsulfate_N',units='kg/hr',element='Sankey')
        def get_ammoniumsulfate_N():
            return MemDis.outs[0].F_mass*14.0067*2/132.14
        
        @metric(name='MemDis_ww_C',units='kg/hr',element='Sankey')
        def get_MemDis_ww_C():
            return MemDis.outs[1].imass['C']
        
        @metric(name='MemDis_ww_N',units='kg/hr',element='Sankey')
        def get_MemDis_ww_N():
            return MemDis.outs[1].imass['N']
        
        @metric(name='MemDis_ww_P',units='kg/hr',element='Sankey')
        def get_MemDis_ww_P():
            return MemDis.outs[1].imass['P']
        
        # energy metrics
        
        @metric(name='sludge_HHV',units='MJ/kg',element='Sankey')
        def get_sludge_HHV():
            return WWTP.HHV
        
        @metric(name='sludge_E',units='GJ/hr',element='Sankey')
        def get_sludge_E():
            return (WWTP.outs[0].F_mass-WWTP.outs[0].imass['H2O'])*WWTP.HHV/1000
        
        @metric(name='biocrude_E',units='GJ/hr',element='Sankey')
        def get_biocrude_E():
            return HTL.biocrude_HHV*HTL.outs[2].imass['Biocrude']/1000
        
        @metric(name='offgas_E',units='GJ/hr',element='Sankey')
        def get_offgas_E():
            return HTL.outs[3].HHV/1000000
        
        @metric(name='HT_gasoline_E',units='GJ/hr',element='Sankey')
        def get_HT_gasoline_E():
            return D2.outs[0].HHV/1000000
        
        @metric(name='HT_diesel_E',units='GJ/hr',element='Sankey')
        def get_HT_diesel_E():
            return D3.outs[0].HHV/1000000
        
        @metric(name='HT_heavy_oil_E',units='GJ/hr',element='Sankey')
        def get_HT_heavy_oil_E():
            return D3.outs[1].HHV/1000000
        
        @metric(name='HT_gas_E',units='GJ/hr',element='Sankey')
        def get_HT_gas_E():
            return (F2.outs[0].HHV+D1.outs[0].HHV)/1000000
        
        @metric(name='HC_gasoline_E',units='GJ/hr',element='Sankey')
        def get_HC_gasoline_E():
            return D4.outs[0].HHV/1000000
        
        @metric(name='HC_diesel_E',units='GJ/hr',element='Sankey')
        def get_HC_diesel_E():
            return D4.outs[1].HHV/1000000
        
        @metric(name='HC_gas_E',units='GJ/hr',element='Sankey')
        def get_HC_gas_E():
            return F3.outs[0].HHV/1000000
        
        @metric(name='HT_hydrogen_gas_E',units='GJ/hr',element='Sankey')
        def get_HT_hydrogen_gas_E():
            return HT.ins[1].HHV/1000000
        
        @metric(name='HC_hydrogen_gas_E',units='GJ/hr',element='Sankey')
        def get_HC_hydrogen_gas_E():
            return HC.ins[1].HHV/1000000
        
        @metric(name='CHG_gas_E',units='GJ/hr',element='Sankey')
        def get_CHG_gas_E():
            return F1.outs[0].HHV/1000000
        
        # CAPEX metrics (as total installed cost, TIC)
        @metric(name='TIC',units='$',element='TEA')
        def get_TIC():
            return sys.installed_equipment_cost
        
        try:
            SluC = unit.SluC
            @metric(name='HTL_TIC',units='$',element='TEA')
            def get_HTL_TIC():
                return sum(i.installed_cost for i in (SluC, P1, H1, HTL))
        except AttributeError:
            @metric(name='HTL_TIC',units='$',element='TEA')
            def get_HTL_TIC():
                return sum(i.installed_cost for i in (P1, H1, HTL))
        
        if AcidEx:
            if SP1.F_mass_out != 0:
                def get_Phosphorus_TIC():
                    return (AcidEx.installed_cost +
                            StruPre.installed_cost +
                            H2SO4_Tank.installed_cost*SP1.outs[0].F_mass/SP1.F_mass_out)
            else:
                def get_Phosphorus_TIC():
                    return 0
        else:
            if SP1.F_mass_out != 0:
                def get_Phosphorus_TIC():
                    return (StruPre.installed_cost +
                            H2SO4_Tank.installed_cost*SP1.outs[0].F_mass/SP1.F_mass_out)
            else:
                def get_Phosphorus_TIC():
                    return 0
        metric(getter=get_Phosphorus_TIC, name='Phosphorus_TIC',units='$',element='TEA')
        
        @metric(name='CHG_TIC',units='$',element='TEA')
        def get_CHG_TIC():
            return CHG.installed_cost + F1.installed_cost
        
        if SP1.F_mass_out != 0:
            def get_Nitrogen_TIC():
                return MemDis.installed_cost+H2SO4_Tank.installed_cost*SP1.outs[1].F_mass/SP1.F_mass_out
        else:
            def get_Nitrogen_TIC():
                return 0
        metric(getter=get_Nitrogen_TIC, name='Nitrogen_TIC',units='$',element='TEA')
        
        @metric(name='HT_HC_TIC',units='$',element='TEA')
        def get_HT_HC_TIC():
            return sum(i.installed_cost for i in HT_HC_units)
        
        @metric(name='HXN_TIC',units='$',element='TEA')
        def get_HXN_TIC():
            return HXN.installed_cost
        
        @metric(name='CHP_TIC',units='$',element='TEA')
        def get_CHP_TIC():
            return CHP.installed_cost
        
        @metric(name='CT_TIC',units='$',element='TEA')
        def get_CT_TIC():
            return CT.installed_cost
        
        @metric(name='CWP_TIC',units='$',element='TEA')
        def get_CWP_TIC():
            return CWP.installed_cost
        
        @metric(name='AOC',units='$/yr',element='TEA')
        def get_AOC():
            return tea.AOC
        
        @metric(name='FOC',units='$/yr',element='TEA')
        def get_FOC():
            return tea.FOC
        
        @metric(name='VOC',units='$/yr',element='TEA')
        def get_VOC():
            return tea.VOC
        
        @metric(name='material_VOC',units='$/yr',element='TEA')
        def get_material_VOC():
            return sys.material_cost
        
        @metric(name='NaOH_HTL_VOC',units='$/yr',element='TEA')
        def get_NaOH_HTL_VOC():
            return NaOH_HTL.cost*sys.operating_hours
        
        @metric(name='HCl_VOC',units='$/yr',element='TEA')
        def get_HCl_VOC():
            return HCl.cost*sys.operating_hours
        
        @metric(name='H2SO4_VOC',units='$/yr',element='TEA')
        def get_H2SO4_VOC():
            return H2SO4.cost*sys.operating_hours
        
        @metric(name='hydrogen_gas_VOC',units='$/yr',element='TEA')
        def get_hydrogen_gas_VOC():
            return hydrogen_gas.cost*sys.operating_hours
        
        @metric(name='CHG_catalyst_VOC',units='$/yr',element='TEA')
        def get_CHG_catalyst_VOC():
            return CHG_catalyst_in.cost*sys.operating_hours
        
        @metric(name='MgCl2_VOC',units='$/yr',element='TEA')
        def get_MgCl2_VOC():
            return MgCl2.cost*sys.operating_hours
        
        @metric(name='HT_catalyst_VOC',units='$/yr',element='TEA')
        def get_HT_catalyst_VOC():
            return CoMo_alumina_HT.cost*sys.operating_hours
        
        @metric(name='HC_catalyst_VOC',units='$/yr',element='TEA')
        def get_HC_catalyst_VOC():
            return CoMo_alumina_HC.cost*sys.operating_hours
        
        @metric(name='MgO_VOC',units='$/yr',element='TEA')
        def get_MgO_VOC():
            return MgO.cost*sys.operating_hours
        
        @metric(name='NaOH_MemDis_VOC',units='$/yr',element='TEA')
        def get_NaOH_MemDis_VOC():
            return NaOH_MemDis.cost*sys.operating_hours
        
        @metric(name='membrane_VOC',units='$/yr',element='TEA')
        def get_membrane_VOC():
            return Membrane_in.cost*sys.operating_hours
        
        @metric(name='natural_gas_VOC',units='$/yr',element='TEA')
        def get_natural_gas_VOC():
            return natural_gas.cost*sys.operating_hours
        
        @metric(name='cooling_tower_chemicals_VOC',units='$/yr',element='TEA')
        def get_cooling_tower_chemicals_VOC():
            return cooling_tower_chemicals.cost*sys.operating_hours
        
        @metric(name='utility_VOC',units='$/yr',element='TEA')
        def get_utility_VOC():
            return sys.utility_cost
        
        try:
            SluC = unit.SluC
            @metric(name='HTL_utility_VOC',units='$/yr',element='TEA')
            def get_HTL_utility_VOC():
                return (SluC.utility_cost+P1.utility_cost+H1.utility_cost+HTL.utility_cost)*sys.operating_hours
        except AttributeError:
            @metric(name='HTL_utility_VOC',units='$/yr',element='TEA')
            def get_HTL_utility_VOC():
                return (P1.utility_cost+H1.utility_cost+HTL.utility_cost)*sys.operating_hours
        
        @metric(name='CHG_utility_VOC',units='$/yr',element='TEA')
        def get_CHG_utility_VOC():
            return (CHG.utility_cost+F1.utility_cost)*sys.operating_hours
        
        @metric(name='HT_HC_utility_VOC',units='$/yr',element='TEA')
        def get_HT_HC_utility_VOC():
            return sum(i.utility_cost for i in HT_HC_units)*sys.operating_hours
        
        @metric(name='HXN_utility_VOC',units='$/yr',element='TEA')
        def get_HXN_utility_VOC():
            return HXN.utility_cost*sys.operating_hours
        
        @metric(name='CHP_utility_VOC',units='$/yr',element='TEA')
        def get_CHP_utility_VOC():
            return CHP.utility_cost*sys.operating_hours
        
        @metric(name='CT_utility_VOC',units='$/yr',element='TEA')
        def get_CT_utility_VOC():
            return CT.utility_cost*sys.operating_hours
        
        @metric(name='CWP_utility_VOC',units='$/yr',element='TEA')
        def get_CWP_utility_VOC():
            return CWP.utility_cost*sys.operating_hours
        
        # LCA metrics
        
        @metric(name='stream_GWP',units='kg CO2 eq',element='LCA')
        def get_stream_GWP():
            return lca.get_stream_impacts()['GlobalWarming']
        
        @metric(name='other_GWP',units='kg CO2 eq',element='LCA')
        def get_other_GWP():
            return lca.get_other_impacts()['GlobalWarming']
        
        @metric(name='HTL_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_HTL_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['NaOH_HTL']+table_stream['HCl']
        
        @metric(name='CHG_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_CHG_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['CHG_catalyst_out']
        
        @metric(name='HT_HC_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_HT_HC_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['hydrogen_gas']+table_stream['HT_catalyst_out']+table_stream['HC_catalyst_out']+table_stream['gasoline']+table_stream['diesel']
        
        @metric(name='nutrient_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_nutrient_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['H2SO4']+table_stream['MgCl2']+table_stream['MgO']+table_stream['struvite']+\
                   table_stream['NH4Cl']+table_stream['Membrane_in']+table_stream['NaOH_MemDis']+table_stream['ammonium_sulfate']
                   
        @metric(name='CHP_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_CHP_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['natural_gas']
        
        @metric(name='CT_stream_GWP',units='kg CO2 eq',element='LCA')
        def get_CT_stream_GWP():
            table_stream = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_stream['CT_chemicals [kg]']
        
        try:
            SluC = unit.SluC
            @metric(name='HTL_utility_GWP',units='kg CO2 eq',element='LCA')
            def get_HTL_utility_GWP():
                table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
                return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                       (SluC.power_utility.consumption+P1.power_utility.consumption)*sys.operating_hours
        
        except AttributeError:
            @metric(name='HTL_utility_GWP',units='kg CO2 eq',element='LCA')
            def get_HTL_utility_GWP():
                table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
                return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                       (P1.power_utility.consumption)*sys.operating_hours    
    
        @metric(name='CHG_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_CHG_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                   CHG.power_utility.consumption*sys.operating_hours
        
        @metric(name='HT_HC_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_HT_HC_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                   (P2.power_utility.consumption+P3.power_utility.consumption)*sys.operating_hours
        
        @metric(name='CHP_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_CHP_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                   (-CHP.power_utility.production)*sys.operating_hours
        
        @metric(name='CT_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_CT_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                   CT.power_utility.consumption*sys.operating_hours
                   
        @metric(name='CWP_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_CWP_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']/(sys.get_electricity_consumption()-sys.get_electricity_production())*\
                   CWP.power_utility.consumption*sys.operating_hours
        
        @metric(name='electricity_GWP',units='kg CO2 eq',element='LCA')
        def get_electricity_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity [kWh]']
        
        @metric(name='HTL_cooling_percentage',units='-',element='utilities')
        def get_HTL_cooling_percentage():
            HTL_cool=0
            for i in range(len(unit.HTL.heat_utilities)):
                if unit.HTL.heat_utilities[i].duty < 0:
                    HTL_cool+=unit.HTL.heat_utilities[i].duty
                    
            sys_cool=0
            for i in range(len(sys.heat_utilities)):
                if sys.heat_utilities[i].duty < 0:
                    sys_cool+=sys.heat_utilities[i].duty
                    
            HXN_cool=0
            for i in range(len(unit.HXN.heat_utilities)):
                if unit.HXN.heat_utilities[i].duty > 0:
                    HXN_cool+=unit.HXN.heat_utilities[i].duty
            
            return HTL_cool/(sys_cool-HXN_cool)
        
        @metric(name='CHG_cooling_percentage',units='-',element='utilities')
        def get_CHG_cooling_percentage():
            CHG_cool=0
            for i in range(len(unit.CHG.heat_utilities)):
                if unit.CHG.heat_utilities[i].duty < 0:
                    CHG_cool+=unit.CHG.heat_utilities[i].duty
                    
            sys_cool=0
            for i in range(len(sys.heat_utilities)):
                if sys.heat_utilities[i].duty < 0:
                    sys_cool+=sys.heat_utilities[i].duty
                    
            HXN_cool=0
            for i in range(len(unit.HXN.heat_utilities)):
                if unit.HXN.heat_utilities[i].duty > 0:
                    HXN_cool+=unit.HXN.heat_utilities[i].duty
            
            return CHG_cool/(sys_cool-HXN_cool)
        
        @metric(name='HTL_heating_percentage',units='-',element='utilities')
        def get_HTL_heating_percentage():
            HTL_heat=0
            for i in range(len(unit.H1.heat_utilities)):
                if unit.H1.heat_utilities[i].duty > 0:
                    HTL_heat+=unit.H1.heat_utilities[i].duty
                    
            sys_heat=0
            for i in range(len(sys.heat_utilities)):
                if sys.heat_utilities[i].duty > 0:
                    sys_heat+=sys.heat_utilities[i].duty
                    
            HXN_heat=0
            for i in range(len(unit.HXN.heat_utilities)):
                if unit.HXN.heat_utilities[i].duty < 0:
                    HXN_heat+=unit.HXN.heat_utilities[i].duty
                    
            CHP_heat=0
            for i in range(len(unit.CHP.heat_utilities)):
                if unit.CHP.heat_utilities[i].duty < 0:
                    CHP_heat+=unit.CHP.heat_utilities[i].duty
            
            return HTL_heat/(sys_heat-HXN_heat-CHP_heat)
        
        @metric(name='CHG_heating_percentage',units='-',element='utilities')
        def get_CHG_heating_percentage():
            CHG_heat=0
            for i in range(len(unit.CHG.heat_utilities)):
                if unit.CHG.heat_utilities[i].duty > 0:
                    CHG_heat+=unit.CHG.heat_utilities[i].duty
                    
            sys_heat=0
            for i in range(len(sys.heat_utilities)):
                if sys.heat_utilities[i].duty > 0:
                    sys_heat+=sys.heat_utilities[i].duty
                    
            HXN_heat=0
            for i in range(len(unit.HXN.heat_utilities)):
                if unit.HXN.heat_utilities[i].duty < 0:
                    HXN_heat+=unit.HXN.heat_utilities[i].duty
                    
            CHP_heat=0
            for i in range(len(unit.CHP.heat_utilities)):
                if unit.CHP.heat_utilities[i].duty < 0:
                    CHP_heat+=unit.CHP.heat_utilities[i].duty
            
            return CHG_heat/(sys_heat-HXN_heat-CHP_heat)
        
        @metric(name='HXN_heat_offset',units='-',element='utilities')
        def get_HXN_heat_offset():
            return 1-unit.HXN.actual_heat_util_load/unit.HXN.original_heat_util_load
        
        @metric(name='HXN_cool_offset',units='-',element='utilities')
        def get_HXN_cool_offset():
            return 1-unit.HXN.actual_cool_util_load/unit.HXN.original_cool_util_load
    
    if include_HTL_yield_as_metrics:
        @metric(name='HTL_biocrude_yield',units='-',element='HTL')
        def get_HTL_biocrude_yield():
            return HTL.biocrude_yield
        
        @metric(name='HTL_aqueous_yield',units='-',element='HTL')
        def get_HTL_aqueous_yield():
            return HTL.aqueous_yield
        
        @metric(name='HTL_hydrochar_yield',units='-',element='HTL')
        def get_HTL_hydrochar_yield():
            return HTL.hydrochar_yield
        
        @metric(name='HTL_gas_yield',units='-',element='HTL')
        def get_HTL_gas_yield():
            return HTL.gas_yield
    
    # key metrics
    @metric(name='MDSP',units='$/gal diesel',element='key_metric')
    def get_MDSP():
        diesel_gal_2_kg=3.220628346
        return tea.solve_price(diesel)*diesel_gal_2_kg
    
    @metric(name='sludge_management_price',units='$/tonne dry sludge',element='key_metric')
    def get_sludge_treatment_price():
        return -tea.solve_price(raw_wastewater)*_MMgal_to_L/WWTP.ww_2_dry_sludge
    
    @metric(name='GWP_diesel',units='kg CO2/MMBTU diesel',element='key_metric')
    def get_GWP_diesel():
        return lca.get_total_impacts(exclude=(diesel,))['GlobalWarming']/diesel.F_mass/sys.operating_hours/lca.lifetime/45.6/_MJ_to_MMBTU
    # gasoline : 46.4 MJ/kg
    # diesel: 45.6 MJ/kg
    
    @metric(name='GWP_sludge',units='kg CO2/tonne dry sludge',element='key_metric')
    def get_GWP_sludge():
        return lca.get_total_impacts(exclude=(raw_wastewater,))['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)/lca.lifetime
        
    if include_check:
        @metric(name='sludge_afdw_carbohydrate',units='-',element='test')
        def get_sludge_afdw_carbohydrate():
            return unit.WWTP.afdw_carbo
    
    return model