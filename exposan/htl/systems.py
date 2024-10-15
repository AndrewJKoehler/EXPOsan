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

References:
[1] https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20).
[2] https://data.bls.gov/cgi-bin/srgate (accessed 2024-08-06).
[3] Snowden-Swan, L. J.; Li, S.; Thorson, M. R.; Schmidt, A. J.; Cronin, D. J.;
    Zhu, Y.; Hart, T. R.; Santosa, D. M.; Fox, S. P.; Lemmon, T. L.; Swita, M. S.
    Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2022 State of Technology; PNNL-33622; Pacific Northwest National Lab. (PNNL),
    Richland, WA (United States), 2022. https://doi.org/10.2172/1897670.
[4] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. Process Design and Economics for
    the Conversion of Algal Biomass to Hydrocarbons: Whole Algae Hydrothermal
    Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; p PNNL--23227, 1126336.
    https://doi.org/10.2172/1126336.
[5] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.
[6] Davis, R.; Hawkins, T.; Coleman, A.; Gao, S.; Klein, B.; Wiatrowski, M.;
    Zhu, Y.; Xu, Y.; Snowden-Swan, L.; Valdez, P.; Zhang, J.; Singh, U.; Ou, L.
    Economic, Greenhouse Gas, and Resource Assessment for Fuel and Protein
    Production from Microalgae: 2022 Algae Harmonization Update. 2024.
    https://doi.org/10.2172/2318964.
[7] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.; Beckham, G. T.;
    Humbird, D.; Thompson, D. N.; Roni, M. S. Process Design and Economics for the
    Conversion of Lignocellulosic Biomass to Hydrocarbon Fuels and Coproducts:
    2018 Biochemical Design Case Update; Biochemical Deconstruction and Conversion
    of Biomass to Fuels and Products via Integrated Biorefinery Pathways;
    NREL/TP--5100-71949, 1483234; 2018; p NREL/TP--5100-71949, 1483234.
    https://doi.org/10.2172/1483234.
[8] https://www.chemanalyst.com/Pricing-data/magnesium-chloride-1403 (accessed 2024-09-06)
[9] https://businessanalytiq.com/procurementanalytics/index/ammonium-chloride-price-index/
    (accessed 2024-09-06)
[10] https://www.alibaba.com/showroom/magnesium-oxide-price.html (accessed 2024-09-06)
[11] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
    Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
    2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
    https://doi.org/10.2172/1855835.
[12] Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
     Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
     Thermal Efficiency, Sensitivity Study and Cost Estimation. Journal of Membrane
     Science 2008, 323 (1), 85–98. https://doi.org/10.1016/j.memsci.2008.06.006.
[13] https://businessanalytiq.com/procurementanalytics/index/ammonium-sulfate-index/
     (accessed 2024-08-03).
[14] https://businessanalytiq.com/procurementanalytics/index/gasoline-price-index/
    (accessed 2024-09-06)
[15] https://businessanalytiq.com/procurementanalytics/index/diesel-price-index/
    (accessed 2024-09-06)
[16] https://en.wikipedia.org/wiki/Methane (accessed 2024-08-03).
[17] https://www.plinovodi.si/en/transmission-system/environment-and-safety/
     about-natural-gas/ (accessed 2024-08-08).
[18] https://www.eia.gov/dnav/ng/hist/n3035us3A.htm accessed (2024-08-03).
[19] Snowden-Swan, L. J.; Zhu, Y.; Bearden, M. D.; Seiple, T. E.; Jones, S. B.;
     Schmidt, A. J.; Billing, J. M.; Hallen, R. T.; Hart, T. R.; Liu, J.;
     Albrecht, K. O.; Fox, S. P.; Maupin, G. D.; Elliott, D. C.
     Conceptual Biorefinery Design and Research Targeted for 2022:
     Hydrothermal Liquefacation Processing of Wet Waste to Fuels; PNNL-27186;
     Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2017.
     https://doi.org/10.2172/1415710.
[20] Stewart, D. W.; Cortés-Peña, Y. R.; Li, Y.; Stillwell, A. S.; Khanna, M.;
     Guest, J. S. Implications of Biorefinery Policy Incentives and
     Location-SpecificEconomic Parameters for the Financial Viability of Biofuels.
     Environ. Sci. Technol. 2023. https://doi.org/10.1021/acs.est.2c07936.
'''

import qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from biosteam.units import IsenthalpicValve
from qsdsan.utils import auom, clear_lca_registries
from exposan.htl import (
    _load_components,
    _load_process_settings,
    create_tea,
    )
from exposan.htl import _sanunits as su
from biosteam import settings

__all__ = ('create_system',)

# kg/m3
water_density = 1000

_MMgal_to_L = auom('gal').conversion_factor('L')*1000000
_m3_to_L = auom('m3').conversion_factor('L')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_ft3 = auom('m3').conversion_factor('ft3')
_ton_to_kg = auom('ton').conversion_factor('kg')

# TODO: which cost year do we want to use?
# TODO: also update cost year in the process setting
# use 2022 for now

# GDPCTPI (Gross Domestic Product: Chain-type Price Index), [1]
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

# the labor index can be found in [2] with the series id CEU3232500008,
# remember to select 'include annual average'
labor_index = {2014: 21.49,
               2015: 21.76,
               2016: 22.71,
               2017: 24.29,
               2018: 25.46,
               2019: 25.46,
               2020: 26.03,
               2021: 26.69,
               2022: 27.36,
               2023: 29.77}

def create_system(configuration='baseline',
                  feedstock='sludge',
                  HTL_model='kinetics',
                  capacity=100,
                  rxn_time=60,
                  rxn_temp=350,
                  set_moisture=0.8,
                  NaOH_M=3,
                  # TODO: update the values for waste_cost and waste_GWP
                  waste_cost=500,
                  waste_GWP=200,
                  high_IRR=False,
                  HCl_neutralize=False):
    
    configuration = configuration or 'baseline'
    if configuration not in ('baseline','no_P','PSA'):
        raise ValueError('`configuration` can only be "baseline", '
                         '"no_P" (without acid extraction and P recovery), '
                         'or "PSA" (with H2 recovery through pressure swing adsorption), '
                         f'not "{configuration}".')
    flowsheet_ID = f'htl_{configuration}'
    
    if feedstock not in ['sludge','biosolid']:
        raise ValueError("invalid feedstock, select from 'sludge' and 'biosolid'")
    
    if HTL_model not in ['MCA','kinetics']:
        raise ValueError("invalid feedstock, select from 'MCA' and 'kinetics'")
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    _load_process_settings()
    
    # TODO: change lipids and proteins based on average (Jan 23, 2024)
    # TODO: find values for sludge and biosolids - CITE SOURCES (Jan 30, 2024)
    
    # sludge compositions from [3]
    if feedstock == 'sludge':
        moisture_content = 0.975 #value from Metcalf & Eddy, 2003
        #range from 1.5-4, average is 2.5
        #TODO: check the difference in LCA, TEA if value is changed
        dw_ash_content = 0.231
        afdw_lipid_content = 0.206
        afdw_protein_content = 0.456      
        N_2_P_value = 0.3927
        VSS_reduction = 0
        before_AD_dw_ash = 0.231
        assert dw_ash_content == before_AD_dw_ash
    # TODO: update biosolids compositions
    else:
        moisture_content = 0.8
        dw_ash_content = 0.45
        afdw_lipid_content = 0.097
        afdw_protein_content = 0.324
        N_2_P_value = 0.3927
        # TODO: update these values with Metcalfe and Eddy values
        # value = loss of VSS, average from lit, see excel
        VSS_reduction = 0.385664
        # TODO: make sure before_AD_dw_ash has the same value as the dw_ash_content of sludge
        before_AD_dw_ash = 0.231
 
    if HTL_model == 'kinetics':
        # TODO: is the lignin content the same for sludge and biosolids?
        # TODO: if yes, we do not need the if statement below
        if feedstock == 'sludge':
            afdw_lignin_content = 0.02
        else:
            afdw_lignin_content = 0.231
    # MCA model does not consider lignin
    else:
        afdw_lignin_content = 0
    
    # set H2O equal to the total raw wastewater into the WWTP
    raw_wastewater = qs.WasteStream(ID='raw_wastewater', H2O=capacity, units='MGD', T=25+273.15)
    
    # =========================================================================
    # pretreatment
    # =========================================================================
    # assume methane from AD is sent to CHP
    WWTP = su.WWTP(ID='WWTP', ins=raw_wastewater, outs=('sludge','treated_water','methane'),
                   # how much metric ton/day sludge can be produced by 1 MGD of ww
                   ww_2_dry_sludge=1,
                   moisture=moisture_content,
                   dw_ash=dw_ash_content, 
                   afdw_lipid=afdw_lipid_content, 
                   afdw_lignin=afdw_lignin_content,
                   afdw_protein=afdw_protein_content,
                   N_2_P=N_2_P_value,
                   feedstock = feedstock,
                   VSS_reduction=VSS_reduction,
                   before_AD_dw_ash=before_AD_dw_ash,
                   operation_hours=7920,
                   high_IRR=high_IRR)
    
    raw_wastewater.price = -WWTP.ww_2_dry_sludge*waste_cost/_MMgal_to_L*_m3_to_L/water_density
    
    if WWTP.moisture <= set_moisture:
        # assume make_up water is the effluent of the WWTP, therefore, no cost and no CI
        Humidifier = su.Humidifier(ID='Humidifier', ins=(WWTP-0, 'makeup_water'), outs='HTL_influent', set_moisture = set_moisture)
        
        # 3049.7 psia, [4]
        P1 = qsu.SludgePump(ID='P1', ins=Humidifier-0, outs='press_sludge', P=3049.7*6894.76,
                  init_with='Stream')
        P1.include_construction = False
    else:    
        SluC = qsu.SludgeCentrifuge(ID='SluC', ins=WWTP-0,
                                    outs=('supernatant','compressed_sludge'),
                                    init_with='Stream',
                                    solids=('Sludge_lipid','Sludge_protein',
                                            'Sludge_carbo','Sludge_ash','Sludge_lignin'),
                                    sludge_moisture=set_moisture)
        
        # 3049.7 psia, [4]
        P1 = qsu.SludgePump(ID='P1', ins=SluC-1, outs='press_sludge', P=3049.7*6894.76,
                            init_with='Stream')
        P1.include_construction = False
    
    # =========================================================================
    # HTL
    # =========================================================================
    
    # TODO: add AD as a sanunit if necessary
    # TODO: from Jeremy: add heat loss for AD; add here by setting heat_transfer_efficiency?
    # U: 3, 3.5, 4 BTU/hr/ft2/F as minimum, baseline, and maximum (case B in [5])
    # unit conversion can be found in the original HTL model
    H1 = qsu.HXutility(ID='H1', include_construction=False,
                       ins=P1-0, outs='heated_sludge', T=rxn_temp+273.15,
                       U=0.0198739, init_with='Stream', rigorous=True,
                       heat_transfer_efficiency=1)
  
    HCl_Tank = qsu.StorageTank(ID='HCl_Tank', ins='HCl', outs='HCl_out',
                              init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.49 2020$/lb, [6]
    HCl_Tank.ins[0].price = 0.49/_lb_to_kg/GDPCTPI[2020]*GDPCTPI[2022]
    
    # create a NaOH stream as 'WasteStream' to enable LCA (this is a fake unit with no construction)
    NaOH_WS = su.preWSConverter(ID='NaOH_WS', ins='NaOH_HTL', outs='NaOH_out', init_with='WasteStream')
    # 0.2384 2016$/lb, [7]
    NaOH_WS.ins[0].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    HTL = qsu.HydrothermalLiquefaction(ID='HTL',
                                       ins=(H1-0, NaOH_WS-0, 'PFAS_in', HCl_Tank-0),
                                       outs=('hydrochar','HTL_aqueous','biocrude','offgas_HTL'),
                                       HTL_model=HTL_model,
                                       rxn_moisture=set_moisture,
                                       NaOH_M=NaOH_M,
                                       feedstock=feedstock,
                                       rxn_time=rxn_time,
                                       rxn_temp=rxn_temp,
                                       HCl_neut=HCl_neutralize,
                                       mositure_adjustment_exist_in_the_system=True)
    
    # =========================================================================
    # CHG
    # =========================================================================
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # H2SO4 price from [7]
    H2SO4_Tank.ins[0].price = (0.043*1+0.0002*(93/5-1))/(93/5)/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]  
    
    # must put after AcidEx and MemDis in path during simulation to ensure the input is not empty
    SP1 = qsu.ReversedSplitter(ID='SP1',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    
    if configuration == 'no_P':
        M1_ins1 = ''
    else:
        # assume the residuals after AcidEx is neutral for TEA and LCA (therefore, not included)
        AcidEx = su.AcidExtraction(ID='AcidEx', ins=(HTL-0, SP1-0),
                                   outs=('residual','extracted'))
        
        M1_ins1 = AcidEx.outs[1]
    
    M1 = su.HTLmixer(ID='M1', ins=(HTL-1, M1_ins1), outs='mixture')
    
    StruPre = su.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                       outs=('struvite','CHG_feed'))
    # 2022 medium price for magnesium chloride: $535/ton, [8]
    StruPre.ins[1].price = 535/_ton_to_kg
    # 2022 price for ammonium chloride: 0.4475 $/kg, [9]
    StruPre.ins[2].price = 0.4475
    # estimated price: 0.2/kg, [10]
    StruPre.ins[3].price = 0.2
    # 0.30 2016$/lb, [11]
    StruPre.outs[0].price = 0.30/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]

    CHG = qsu.CatalyticHydrothermalGasification(ID='CHG',
                                                ins=(StruPre-1, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out', 'used_CHG_catalyst'))
    # 60 2011$/lb, [4]
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2022]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.include_construction = False
    
    # assume the pH of the influent of MemDis is independent from HTL as there are CHG and a flash unit in between
    # and the NaOH in MemDis is relatively small (~0.01%) compared to even 1 M of NaOH in HTL
    # assume no value/cost and no environmental benefit/impact associated with MemDis_ww (the system is in a WWTP) and solution
    MemDis = qsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH_MemDis', 'Membrane_in'),
                                      outs=('ammonium_sulfate','MemDis_ww','Membrane_out','solution'), init_with='WasteStream')
    # 0.2384 2016$/lb, [7]
    MemDis.ins[2].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    # RO membrane price: [12] (likely 2008$)
    MemDis.ins[3].price = 90/GDPCTPI[2008]*GDPCTPI[2022]
    # ammonium sulfate (2022 average): [13]
    MemDis.outs[0].price = 0.2825
    MemDis.include_construction = False
    
    # =========================================================================
    # HT
    # =========================================================================
    
    P2 = qsu.SludgePump(ID='P2', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
                        init_with='Stream')
    # 1530.0 psia, [4]
    P2.include_construction = False
    
    # reversed splitter, write before HT and HC, simulate after HT and HC
    RSP1 = qsu.ReversedSplitter(ID='RSP1', ins='hydrogen_gas', outs=('HT_H2','HC_H2'),
                                init_with='WasteStream')
    # 0.7306 2016$/lb, [7]
    RSP1.ins[0].price = 0.7306/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    include_PSA = False if 'PSA' not in configuration else True
    HT = qsu.Hydrotreating(ID='HT', ins=(P2-0, RSP1-0, 'CoMo_alumina_HT'),
                           outs=('HTout','CoMo_alumina_HT_out'), include_PSA=include_PSA)
    # 15.5 2007$/lb, [4]
    HT.ins[2].price = 15.5/_lb_to_kg/GDPCTPI[2007]*GDPCTPI[2022]
    
    V2 = IsenthalpicValve(ID='V2', ins=HT-0, outs='depressed_HT', P=717.4*6894.76, vle=True)
    
    H2 = qsu.HXutility(ID='H2', include_construction=False,
                       ins=V2-0, outs='cooled_HT', T=60+273.15,
                       init_with='Stream', rigorous=True)
    
    F2 = qsu.Flash(ID='F2', ins=H2-0, outs=('HT_fuel_gas','HT_aqueous'), T=43+273.15,
                   P=717.4*6894.76, thermo=settings.thermo.ideal())
    F2.include_construction = False
    
    V3 = IsenthalpicValve(ID='V3', ins=F2-1, outs='depressed_flash_effluent', P=55*6894.76, vle=True)
    
    # separate water and oil based on gravity
    SP2 = qsu.Splitter(ID='SP2', ins=V3-0, outs=('HT_ww','HT_oil'),
                       split={'H2O':1}, init_with='Stream')
    
    # temperature from stream #334, [4]
    # the first distillation column was removed
    H3 = qsu.HXutility(ID='H3', include_construction=False,
                       ins=SP2-1, outs='heated_oil', T=104+273.15, rigorous=True)
    
    D1 = qsu.BinaryDistillation(ID='D1', ins=H3-0,
                                outs=('HT_light','HT_heavy'),
                                LHK=('C4H10','TWOMBUTAN'), P=50*6894.76,
                                y_top=188/253, x_bot=53/162, k=2, is_divided=True)
    D1.include_construction = False
    
    D2 = qsu.BinaryDistillation(ID='D2', ins=D1-1,
                                outs=('HT_Gasoline','HT_other_oil'),
                                LHK=('C10H22','C4BENZ'), P=25*6894.76,
                                y_top=116/122, x_bot=114/732, k=2, is_divided=True)
    D2.include_construction = False
    
    D3 = qsu.BinaryDistillation(ID='D3', ins=D2-1,
                                outs=('HT_Diesel','HT_heavy_oil'),
                                LHK=('C19H40','C21H44'),P=18.7*6894.76,
                                y_top=2421/2448, x_bot=158/2448, k=2, is_divided=True)
    D3.include_construction = False
    
    # =========================================================================
    # HC
    # =========================================================================
    
    # 1034.7 psia, [4]
    P3 = qsu.SludgePump(ID='P3', ins=D3-1, outs='press_heavy_oil', P=1034.7*6894.76,
                        init_with='Stream')
    P3.include_construction = False
    
    HC = qsu.Hydrocracking(ID='HC', ins=(P3-0, RSP1-1, 'CoMo_alumina_HC'),
                           outs=('HC_out','CoMo_alumina_HC_out'))
    # 15.5 2007$/lb, [4]
    HC.ins[2].price = 15.5/_lb_to_kg/GDPCTPI[2007]*GDPCTPI[2022]
    
    H4 = qsu.HXutility(ID='H4', include_construction=False,
                       ins=HC-0, outs='cooled_HC', T=60+273.15,
                       init_with='Stream', rigorous=True)
    
    V4 = IsenthalpicValve(ID='V4', ins=H4-0, outs='cooled_depressed_HC', P=30*6894.76, vle=True)
    
    F3 = qsu.Flash(ID='F3', ins=V4-0, outs=('HC_fuel_gas','HC_aqueous'), T=60.2+273,
                   P=30*6894.76)
    F3.include_construction = False
    
    D4 = qsu.BinaryDistillation(ID='D4', ins=F3-1, outs=('HC_Gasoline','HC_Diesel'),
                                LHK=('C9H20','C10H22'), P=20*6894.76,
                                y_top=360/546, x_bot=7/708, k=2, is_divided=True)
    D4.include_construction = False
    
    # =========================================================================
    # storage & disposal
    # =========================================================================
    
    GasolineMixer = qsu.Mixer(ID='GasolineMixer', ins=(D2-0, D4-0), outs='mixed_gasoline',
                              init_with='Stream', rigorous=True)
    
    DieselMixer = qsu.Mixer(ID='DieselMixer', ins=(D3-0, D4-1), outs='mixed_diesel',
                            init_with='Stream', rigorous=True)
    
    H5 = qsu.HXutility(ID='H5', include_construction=False,
                       ins=GasolineMixer-0, outs='cooled_gasoline',
                       T=60+273.15, init_with='Stream', rigorous=True)
    
    H6 = qsu.HXutility(ID='H6', include_construction=False,
                       ins=DieselMixer-0, outs='cooled_diesel',
                       T=60+273.15, init_with='Stream', rigorous=True)
    
    PC1 = qsu.PhaseChanger(ID='PC1', ins=H5-0, outs='cooled_gasoline_liquid')
    
    PC2 = qsu.PhaseChanger(ID='PC2', ins=H6-0, outs='cooled_diesel_liquid')
    
    PC3 = qsu.PhaseChanger(ID='PC3', ins=CHG-1, outs='CHG_catalyst_out', phase='s')
    
    PC4 = qsu.PhaseChanger(ID='PC4', ins=HT-1, outs='HT_catalyst_out', phase='s')
    
    PC5 = qsu.PhaseChanger(ID='PC5', ins=HC-1, outs='HC_catalyst_out', phase='s')
    
    # store for 3 days, [4]
    GasolineTank = qsu.StorageTank(ID='GasolineTank', ins=PC1-0, outs='gasoline',
                                   tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    # 2022 average price, 1.177$/kg, [14]
    GasolineTank.outs[0].price = 1.177
    
    # store for 3 days, [4]
    DieselTank = qsu.StorageTank(ID='DieselTank', ins=PC2-0, outs=('diesel'),
                                 tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    # 2022 average price, 1.331$/kg, [15]
    DieselTank.outs[0].price = 1.331
    
    GasMixer = qsu.Mixer(ID='GasMixer', ins=(WWTP-2, HTL-3, F1-0, F2-0, D1-0, F3-0),
                         outs=('fuel_gas'), init_with='Stream')
    
    # assume the effluent of WWmixer goes back to WWTP
    try:
        WWmixer = su.WWmixer(ID='WWmixer', ins=(SluC-0, MemDis-1, SP2-0),
                             outs='wastewater', init_with='Stream')
    except UnboundLocalError:
        WWmixer = su.WWmixer(ID='WWmixer', ins=('', MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')
    
    # =========================================================================
    # facilities
    # =========================================================================
    
    # assume T_min_app to be 86 K, [4]
    qsu.HeatExchangerNetwork(ID='HXN', T_min_app=86, force_ideal_thermo=True)
    
    # assume no value/cost and no environmental benefit/impact associated with emission (they are not treated and are biogenic; only CO2 from natural gas is non-biogenic but we have included the environmental impact for it)
    # use CHP to produce electricity does not provide benefit; therefore, set supplement_power_utility=False
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=False,
                                ins=(GasMixer-0, 'natural_gas', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    # assume natural gas is CH4 (density is 0.657 kg/m3, [16])
    # the density of natural gas (0.68 kg/m3, [17]) is larger than that of CH4, therefore, the calculation here is conservative
    # 2022 US NG industrial price: 7.66 $/1000 ft3, [18]
    CHP.ins[1].price = 7.66/1000*_m3_to_ft3/0.657
    # 1.41 MM 2016$/year for 4270/4279 kg/hr ash, 7880 annual operating hours, from [7]
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2022]
    
    # construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_chemicals: 1.7842 2016$/lb, [7]
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    # CWP only have one influent (recirculated_chilled_water)
    # assume no cost and environmental impacts for this stream
    CWP = bst.facilities.ChilledWaterPackage(ID='CWP')
    
    sys = qs.System.from_units(ID=f'sys_{configuration}',
                               units=list(flowsheet.unit), 
                               operating_hours=WWTP.operation_hours)
    sys.register_alias('sys')
    
    sys.simulate()
    
    # =========================================================================
    # LCA
    # =========================================================================
    
    # for now, just include GWP, add others if necessary
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.48748859)
    
    # deionized water
    CT_chemicals = qs.ImpactItem('CT_chemicals', functional_unit='kg')
    CT_chemicals.add_indicator(GlobalWarming, 0.00042012744)
    
    impact_items = {'raw_wastewater':    [stream.raw_wastewater, -WWTP.ww_2_dry_sludge*waste_GWP/_MMgal_to_L*_m3_to_L/water_density],
                    'HCl':               [stream.HCl, 0.56248255],
                    'NaOH_HTL':          [stream.NaOH_HTL, 1.2497984],
                    'H2SO4':             [stream.H2SO4, 0.005529872568],
                    'MgCl2':             [stream.MgCl2, 0],
                    'NH4Cl':             [stream.NH4Cl, 1.5195637],
                    'MgO':               [stream.MgO, 1.1605114],
                    'struvite':          [stream.struvite, -0.415531665577197],
                    'CHG_catalyst':      [stream.CHG_catalyst_out, 471.098936962268],
                    'NaOH_MemDis':       [stream.NaOH_MemDis, 1.2497984],
                    'RO_membrane':       [stream.Membrane_in, 2.2709135],
                    'NH42SO4':           [stream.ammonium_sulfate, -1.1139067],
                    'hydrogen_gas':      [stream.hydrogen_gas, 1.5621537],
                    'HT_catalyst':       [stream.HT_catalyst_out, 5.94738828207449],
                    'HC_catalyst':       [stream.HC_catalyst_out, 5.94738828207449],
                    'gasoline':          [stream.gasoline, -0.33231451],
                    'diesel':            [stream.diesel, -0.44231486],
                    'natural_gas':       [stream.natural_gas, 2.5199928],
                    'ash_disposal':      [stream.solid_ash, 0.0082744841]}
    
    for item in impact_items.items():
        qs.StreamImpactItem(ID=item[0], linked_stream=item[1][0], GlobalWarming=item[1][1])
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
           CT_chemicals=lambda:CT.ins[2].F_mass*sys.flowsheet.WWTP.operation_hours*30)
    
    # =========================================================================
    # TEA
    # =========================================================================
    
    # based on the labor cost for the HTL plant from [19], 2014 level:
    # 1 plant manager (0.15 MM$/year)
    # 1 plant engineer (0.07 MM$/year)
    # 1 maintenance supervisor (0.06 MM$year)
    # 1 lab manager (0.06 MM$year)
    # variable cost (proportional to the sludge amount, the following is for a
    # plant of 110 dry ton [100 dry metric tonne] sludge per day):
    # 3 shift supervisors (0.14 MM$/year)
    # 1 lab technican (0.04 MM$/year)
    # 1 maintenance technician (0.04 MM$/year)
    # 4 shift operators (0.19 MM$/year)
    # 1 yard employee (0.03 MM$/year)
    # 1 clerk & secretary (0.04 MM$/year)
    # annual wage [$/year]
    wage = (0.34/labor_index[2014]*labor_index[2022]+\
            0.48/labor_index[2014]*labor_index[2022]*capacity*WWTP.ww_2_dry_sludge/100)*10**6
    
    federal_income_tax_rate_value = 0.21
    # use the mode state income tax (6.5%) from [20]
    state_income_tax_rate_value = 0.065
    income_tax = federal_income_tax_rate_value + state_income_tax_rate_value
    
    if high_IRR:
        create_tea(sys, IRR_value=0.1, income_tax_value=income_tax, finance_interest_value=0.08, labor_cost_value=wage)
    else:
        create_tea(sys, IRR_value=0.03, income_tax_value=income_tax, finance_interest_value=0.03, labor_cost_value=wage)
    
    return sys