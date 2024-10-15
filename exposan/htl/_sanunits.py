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

# import biosteam as bst
from math import ceil, log, pi
from qsdsan import SanUnit
from qsdsan.sanunits import Reactor
from qsdsan.utils import ospath, auom, data_path, load_data
from qsdsan.processes import Decay

__all__ = (
    'AcidExtraction',
    'FuelMixer',
    'HTLmixer',
    'Humidifier',
    'StruvitePrecipitation',
    'WWmixer',
    'WWTP',
    'preWSConverter'
    )

yearly_operation_hour = 7920 # Jones
_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')

# =============================================================================
# Anaerobic Digestion
# =============================================================================

ad_path = ospath.join(data_path, 'sanunit_data/_anaerobic_digestion.tsv')

class AnaerobicDigestion(SanUnit, Decay):
    '''
    Anaerobic digestion of wastes with the production of biogas based on
    `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    To enable life cycle assessment, the following impact items should be pre-constructed:
    `Concrete`, `Excavation`.

    Cost is calculated by the unit cost of the impact items and their quantities.

    Parameters
    ----------
    ins : Iterable
        Waste for treatment.
    outs : Iterable
        Treated waste, captured biogas, fugitive CH4, and fugitive N2O.
    flow_rate : float
        Total flow rate through the reactor (for sizing purpose), [m3/d].
        If not provided, will use F_vol_in.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.

    See Also
    --------
    :ref:`qsdsan.processes.Decay <processes_Decay>`
    '''
    _N_ins = 1
    _N_outs = 4
    _run = Decay._first_order_run
    _units = {
        'Volumetric flow rate': 'm3/hr',
        'Residence time': 'd',
        'Single reactor volume': 'm3',
        'Reactor diameter': 'm',
        'Reactor height': 'm'
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 include_construction=True,
                 flow_rate=None, degraded_components=('OtherSS',),
                 if_capture_biogas=True, if_N2O_emission=False,
                 **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1,
                        include_construction=include_construction,
                        degraded_components=degraded_components,
                        if_capture_biogas=if_capture_biogas,
                        if_N2O_emission=if_N2O_emission,)
        self._flow_rate = flow_rate

        data = load_data(path=ad_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
### Commented out because LCA not desired for construction - consider turning on
### if LCA for heat/energy use not accounted for
###Note: if later turned on, Self.construction must still be deleted
    # def _init_lca(self):
    #     self.construction = [
    #         Construction('concrete', linked_unit=self, item='Concrete', quantity_unit='m3'),
    #         Construction('excavation', linked_unit=self, item='Excavation', quantity_unit='m3'),
    #         ]
    def _run(self):  pass
    


    def _design(self):
        design = self.design_results
        design['Volumetric flow rate'] = Q = self.flow_rate
        design['Residence time'] = tau = self.tau
        design['Reactor number'] = N = self.N_reactor
        V_tot = Q * tau*24

        # One extra as a backup
        design['Single reactor volume'] = V_single = V_tot/(1-self.headspace_frac)/(N-1)

        # Rx modeled as a cylinder
        design['Reactor diameter'] = D = (4*V_single*self.aspect_ratio/pi)**(1/3)
        design['Reactor height'] = H = self.aspect_ratio * D

        if self.include_construction:
            constr = self.construction
            concrete =  N*self.concrete_thickness*(2*pi/4*(D**2)+pi*D*H)
            constr[0].quantity = concrete
            constr[1].quantity = V_tot # excavation

            self.add_construction()
    

    @property
    def flow_rate(self):
        '''
        [float] Total flow rate through the reactor (for sizing purpose), [m3/d].
        If not provided, will calculate based on F_vol_in.
        '''
        return self._flow_rate if self._flow_rate else self.F_vol_in*24
    @flow_rate.setter
    def flow_rate(self, i):
        self._flow_rate = i

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = i

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during treatment.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = i

    @property
    def N_reactor(self):
        '''[int] Number of reactors, float will be converted to the smallest integer.'''
        return self._N_reactor
    @N_reactor.setter
    def N_reactor(self, i):
        self._N_reactor = ceil(i)

    @property
    def aspect_ratio(self):
        '''[float] Diameter-to-height ratio of the reactor.'''
        return self._aspect_ratio
    @aspect_ratio.setter
    def aspect_ratio(self, i):
        self._aspect_ratio = i

    @property
    def headspace_frac(self):
        '''[float] Fraction of the reactor volume for headspace gas.'''
        return self._headspace_frac
    @headspace_frac.setter
    def headspace_frac(self, i):
        self._headspace_frac = i

    @property
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = i

                         
# =============================================================================
# Acid Extraction
# =============================================================================

class AcidExtraction(Reactor):
    '''
    H2SO4 is added to hydrochar from HTL to extract phosphorus.
    
    Parameters
    ----------
    ins : Iterable(stream)
        hydrochar, acid.
    outs : Iterable(stream)
        residual, extracted.
    acid_vol: float
        1 M H2SO4 to hydrochar ratio: mL/g.
    P_acid_recovery_ratio: float
        The ratio of phosphorus that can be extracted.
        
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
        Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
        2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
        https://doi.org/10.2172/1855835.
    '''
    _N_ins = 2
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', acid_vol=2, P_acid_recovery_ratio=0.8,
                 P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    def _run(self):
        
        hydrochar, acid = self.ins
        residual, extracted = self.outs
        
        self.HTL = self.ins[0]._source
        
        if hydrochar.F_mass <= 0:
            pass
        else:
            if self.HTL.hydrochar_P <= 0:
                residual.copy_like(hydrochar)
            else: 
                acid.imass['H2SO4'] = hydrochar.F_mass*self.acid_vol*1*98.079/1000
                # 1 M H2SO4 acid_vol (2 mL/1 g) hydrochar
                acid.imass['H2O'] = hydrochar.F_mass*self.acid_vol*1.05 -\
                                    acid.imass['H2SO4']
                # 0.5 M H2SO4 density: 1.05 kg/L 
                # https://us.vwr.com/store/product/21221349/null (accessed 2024-09-08)
                
                residual.imass['Residual'] = hydrochar.F_mass - self.ins[0]._source.\
                                             hydrochar_P*self.P_acid_recovery_ratio
                
                extracted.copy_like(acid)
                extracted.imass['P'] = hydrochar.F_mass - residual.F_mass
                # assume just P can be extracted
                
                residual.phase = 's'
                
                residual.T = extracted.T = hydrochar.T
                residual.P = hydrochar.P
                # H2SO4 reacts with hydrochar to release heat and temperature will increase
    
    @property
    def residual_C(self):
        return self.ins[0]._source.hydrochar_C
    
    @property
    def residual_P(self):
        return self.ins[0]._source.hydrochar_P - self.outs[1].imass['P']
        
    def _design(self):
        self.N = ceil(self.HTL.WWTP.ins[0].F_vol/788.627455/self.V)
        # 1/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 10 m3)
        self.P = self.ins[1].P
        Reactor._design(self)

# =============================================================================
# FuelMixer
# =============================================================================

class FuelMixer(SanUnit):
    '''
    Convert gasoline to diesel or diesel to gasoline based on LHV.
    
    Parameters
    ----------
    ins: Iterable(stream)
        gasoline, diesel
    outs: Iterable(stream)
        fuel
    target: str
        The target can only be 'gasoline' or 'diesel'.
    gasoline_price: float
        Gasoline price, [$/kg].
    diesel_price: float
        Diesel price, [$/kg].
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target='diesel',
                 gasoline_gal_2_kg=2.834894885,
                 diesel_gal_2_kg=3.220628346,
                 gasoline_price=0.9388,
                 diesel_price=0.9722):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target = target
        self.gasoline_gal_2_kg = gasoline_gal_2_kg
        self.diesel_gal_2_kg = diesel_gal_2_kg
        self.gasoline_price = gasoline_price
        self.diesel_price = diesel_price

    def _run(self):
        
        gasoline, diesel = self.ins
        fuel = self.outs[0]
        target = self.target
        
        gasoline_LHV_2_diesel_LHV = (gasoline.LHV/gasoline.F_mass)/(diesel.LHV/diesel.F_mass)
        # KJ/kg gasoline:KJ/kg diesel
        
        if target == 'gasoline':
            fuel.imass['Gasoline'] = gasoline.F_mass + diesel.F_mass/gasoline_LHV_2_diesel_LHV
            fuel.T = gasoline.T
            fuel.P = gasoline.P
        elif target == 'diesel':
            fuel.imass['Diesel'] = diesel.F_mass + gasoline.F_mass*gasoline_LHV_2_diesel_LHV
            fuel.T = diesel.T
            fuel.P = diesel.P
    
    def _cost(self):
        if self.target == 'gasoline':
            self.outs[0].price = self.gasoline_price
        elif self.target == 'diesel':
            self.outs[0].price = self.diesel_price
            
    @property
    def target(self):
        return self._target
    @target.setter
    def target(self, i):
        if i not in ('gasoline', 'diesel'):
            raise ValueError('`target` must be either "gasoline" or "diesel" ',
                             f'the provided "{i}" is not valid.')
        self._target = i

# =============================================================================
# HTL mixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the mixture of HTL
    aqueous and AcidEx effluent.
    
    Parameters
    ----------
    ins : Iterable(stream)
        HTLaqueous, extracted
    outs : Iterable(stream)
        mixture
        
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
        
    def _run(self):
        
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        mixture.empty()
        
        self.HTL = self.ins[0]._source
        
        mixture.imass['C'] = self.HTL.HTLaqueous_C
        mixture.imass['N'] = self.HTL.HTLaqueous_N
        mixture.imass['P'] = self.HTL.HTLaqueous_P + extracted.imass['P']
        mixture.imass['H2O'] = HTLaqueous.F_mass + extracted.F_mass -\
                               mixture.imass['C'] - mixture.imass['N'] -\
                               mixture.imass['P']
        # represented by H2O except C, N, P
    
    @property
    def pH(self):
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        base_mol_per_h = HTLaqueous.F_vol*1000*10**(self.HTL.aq_pH-14)
        acid_mol_per_h = extracted.imass['H2SO4']*1000/98*2
        
        volume_L_per_h = mixture.F_vol*1000
        
        base_M = base_mol_per_h/volume_L_per_h
        acid_M = acid_mol_per_h/volume_L_per_h
        
        if base_M > acid_M:
            hydrogen_ion_M = 10**-14/(base_M-acid_M)
        elif base_M == acid_M:
            hydrogen_ion_M = 10**(-7)
        else:
            hydrogen_ion_M = acid_M - base_M

        return -log(hydrogen_ion_M, 10)

# =============================================================================
# Humidifier
# =============================================================================

class Humidifier(SanUnit):
    '''
    A fake unit increases the moisture content of HTL feedstocks to a certain moisture content.
    
    Parameters
    ----------
    ins : Iterable(stream)
        feedstock, makeup, recycle
    outs : Iterable(stream)
        mixture
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', set_moisture = 0.8):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.set_moisture = set_moisture

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        if self.set_moisture > 1:
            raise Exception("Moisture content cannot be >1 (above 100%)")
    
        feedstock, makeup = self.ins
        mixture = self.outs[0]
        
        makeup.imass['H2O'] = (feedstock.F_mass - feedstock.imass['H2O'])/(1-self.set_moisture)*self.set_moisture - feedstock.imass['H2O']

        mixture.mix_from(self.ins)
        
    @property
    def moisture_level(self):
        return self.set_moisture

# =============================================================================
# Struvite Precipitation
# =============================================================================

class StruvitePrecipitation(Reactor):
    '''
    Extracted and HTL aqueous are mixed together before adding MgCl2 for struvite precipitation.
    If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1.
    
    Parameters
    ----------
    ins : Iterable(stream)
        mixture, supply_MgCl2, supply_NH4Cl, base.
    outs : Iterable(stream)
        struvite, effluent.
    target_pH: float
        Target pH for struvite precipitation.
    Mg_P_ratio: float
        mol(Mg) to mol(P) ratio.   
    P_pre_recovery_ratio: float
        Ratio of phosphorus that can be precipitated out.
    HTLaqueous_NH3_N_2_total_N: float
        Ratio of NH3-N to TN in HTL aqueous phase.
        
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
        Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
        2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
        https://doi.org/10.2172/1855835.
    .. [2] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 4
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', 
                  target_pH = 9,
                  Mg_P_ratio=1,
                  P_pre_recovery_ratio=0.828, # [1]
                  HTLaqueous_NH3_N_2_total_N = 0.853, # [2]
                  P=None, tau=1, V_wf=0.8, # tau: [1]
                  length_to_diameter=2, N=1, V=20, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                  wall_thickness_factor=1,
                  vessel_material='Carbon steel', # basic condition
                  vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_pH = target_pH
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.HTLaqueous_NH3_N_2_total_N = HTLaqueous_NH3_N_2_total_N
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    def _run(self):
        
        mixture, supply_MgCl2, supply_NH4Cl, base = self.ins
        struvite, effluent = self.outs
        
        self.HTLmixer = self.ins[0]._source
        
        if self.HTLmixer.outs[0].imass['P'] == 0:
            effluent.copy_like(mixture)
        else:
            old_pH = self.HTLmixer.pH
            if old_pH > 9:
                base.imass['MgO'] = 0
            elif old_pH >= 7:
                OH_M = 10**(old_pH-14)
                OH_M_needed = 10**(self.target_pH - 14) - OH_M
                base.imass['MgO'] = OH_M_needed/2*40.3044/1000*self.ins[0].F_vol*1000
            else:
                neutral_OH_M = 10**(-old_pH)
                to_target_OH_M = 10**(self.target_pH - 14)
                OH_M_needed = neutral_OH_M + to_target_OH_M
                base.imass['MgO'] = OH_M_needed/2*40.3044/1000*self.ins[0].F_vol*1000
            
            supply_MgCl2.imass['MgCl2'] = max((mixture.imass['P']/30.973762*self.Mg_P_ratio -\
                                            base.imass['MgO']/40.3044)*95.211, 0)
    
            if mixture.imass['P']/30.973762 > mixture.imass['N']*self.HTLaqueous_NH3_N_2_total_N/14.0067:
            # if P > N, add NH4Cl to make sure N ≥ P
                supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/30.973762 - mixture.imass['N']*\
                                                self.HTLaqueous_NH3_N_2_total_N/14.0067)*53.491

            struvite.imass['Struvite'] = mixture.imass['P']*\
                                          self.P_pre_recovery_ratio/\
                                          30.973762*245.41
            supply_MgCl2.phase = supply_NH4Cl.phase = base.phase = 's'
            
            effluent.copy_like(mixture)
            effluent.imass['P'] -= struvite.imass['Struvite']*30.973762/245.41
            effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14.0067/53.491 -\
                                    struvite.imass['Struvite']*14.0067/245.41
            effluent.imass['H2O'] = self.F_mass_in - struvite.F_mass -\
                                    effluent.imass['C'] - effluent.imass['N'] -\
                                    effluent.imass['P']
            struvite.phase = 's'    
                
            struvite.T = mixture.T
            effluent.T = mixture.T
        
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*30.973762/245.41

    @property
    def struvite_N(self):
        return self.struvite_P*14.0067/30.973762

    def _design(self):
        self.N = ceil(self.HTLmixer.ins[0]._source.WWTP.ins[0].F_vol*2/788.627455/self.V)
        # 2/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 20 m3)
        self.P = self.ins[0].P
        Reactor._design(self)

# =============================================================================
# WWmixer
# =============================================================================

class WWmixer(SanUnit):
    '''
    A fake unit that mixes all wastewater streams and calculates C, N, P, and H2O
    amount.
    Parameters
    ----------
    ins : Iterable(stream)
        supernatant_1, supernatant_2, memdis_ww, ht_ww
    outs : Iterable(stream)
        mixture
    '''
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
    def _run(self):
        
        supernatant, memdis_ww, ht_ww = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        
        HT = self.ins[2]._source.ins[0]._source.ins[0]._source.ins[0]._source.\
             ins[0]._source.ins[0]._source
        
        # only account for C and N from HT if they are not less than 0
        if HT.HTaqueous_C >= 0:
            mixture.imass['C'] += HT.HTaqueous_C
            mixture.imass['H2O'] -= HT.HTaqueous_C
        if HT.HTaqueous_N >=0:
            mixture.imass['N'] += HT.HTaqueous_N
            mixture.imass['H2O'] -= HT.HTaqueous_N

# =============================================================================
# WWTP
# =============================================================================

class WWTP(SanUnit):
    '''
    WWTP is a fake unit that can set up solid biochemical compositions
    and calculate solid elemental compositions.
    
    Parameters
    ----------
    ins : Iterable(stream)
        ww.
    outs : Iterable(stream)
        solid, treated.
    ww_2_dry_sludge: float
        Wastewater-to-dry-sludge conversion factor, [metric ton/day/MGD].
    moisture: float
        Solid moisture content.
    dw_ash: float
        Solid dry weight ash content.
    afdw_lipid: float
        Solid ash free dry weight lipid content.
    afdw_protein: float
        Solid ash free dry weight protein content.
    afdw_lignin: float
        Solid ash free dry weight lignin content.
    lipid_2_C: float
        Lipid to carbon factor.     
    protein_2_C: float
        Protein to carbon factor.
    carbo_2_C: float
        Carbohydrate to carbon factor.
    C_2_H: float
        Carbon to hydrogen factor.
    protein_2_N: float
        Protein to nitrogen factor.
    N_2_P: float
        Nitrogen to phosphorus factor.
    feedstock: str
        Can only be 'sludge' or 'biosolid'.
    VSS_reduction: floar
        VSS reduction ratio after anaerobic digestion.
    before_AD_dw_ash: float
        Solid dry weight ash content before AD.
    operation_hour: float
        Plant yearly operation hour, [hr/yr].
    high_IRR: 
        If True, IRR = 0.1 and finance_interest_value = 0.08,
        If False, IRR = 0.03 and finance_interest_value = 0.03.
    
    References
    ----------
    .. [1] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    '''
    _N_ins = 1
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', 
                 ww_2_dry_sludge=1,
                 moisture=0.99,
                 dw_ash=0.257, 
                 afdw_lipid=0.204,
                 afdw_protein=0.463,
                 afdw_lignin=0,
                 lipid_2_C=0.750,
                 protein_2_C=0.545,
                 carbo_2_C=0.400, 
                 lipid_2_H=0.125,
                 protein_2_H=0.068,
                 carbo_2_H=0.067, 
                 protein_2_N=0.159,
                 N_2_P=0.3927,
                 feedstock='sludge',
                 # TODO: update these values with Metcalf and Eddy values
                 # value = loss of VSS, average from lit, see excel
                 VSS_reduction=0.38,
                 before_AD_dw_ash=0.257,
                 operation_hours=yearly_operation_hour,
                 high_IRR=None):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.moisture = moisture
        self.dw_ash = dw_ash
        self.afdw_lipid = afdw_lipid
        self.afdw_lignin = afdw_lignin
        self.afdw_protein = afdw_protein
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.lipid_2_H = lipid_2_H
        self.protein_2_H = protein_2_H
        self.carbo_2_H = carbo_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.feedstock = feedstock
        self.VSS_reduction = VSS_reduction
        self.before_AD_dw_ash = before_AD_dw_ash
        self.operation_hours = operation_hours
        self.high_IRR = high_IRR
    
    def _run(self):
        
        ww = self.ins[0]
        solid, treated, methane = self.outs
        
        self.afdw_carbo = round(1 - self.afdw_protein - self.afdw_lipid - self.afdw_lignin, 5)   
        
        if self.dw_ash >= 1:
            raise Exception ('ash can not be larger than or equal to 1')
        
        if self.afdw_protein + self.afdw_lipid + self.afdw_lignin > 1:
            raise Exception ('protein and lipid exceed 1')
        
        if self.feedstock == 'biosolid':
            self.dw = ww.F_vol*_m3perh_to_MGD*self.ww_2_dry_sludge*1000/24*(1-(1-self.before_AD_dw_ash)*self.VSS_reduction)
        else:
            self.dw = ww.F_vol*_m3perh_to_MGD*self.ww_2_dry_sludge*1000/24

        solid.imass['H2O'] = self.dw/(1-self.moisture)*self.moisture
        solid.imass['Sludge_ash'] = self.dw*self.dw_ash
        afdw = self.dw*(1 - self.dw_ash)
        
        # anaerobic digestion        
        if self.feedstock == 'biosolid':
            # TODO: double check the calculation here
            # TODO: ask Andrew, double check the density of methane (0.657 kg/m3)
            # TODO: update these values with Metcalfe and Eddy values
            # 0.178108571 L methane/kg VSS, methane density 0.7513 kg methane/L of methane
            methane.imass['CH4'] = afdw*self.VSS_reduction*0.178108571/1000*0.657
            # methane.imass['CH4'] =  afdw * VSS_reduction * 0.178108571 * 0.7513
        else:
            methane.imass['CH4'] = 0   
        
        solid.imass['Sludge_lipid'] = afdw*self.afdw_lipid
        solid.imass['Sludge_protein'] = afdw*self.afdw_protein
        solid.imass['Sludge_carbo'] = afdw*self.afdw_carbo
        solid.imass['Sludge_lignin'] = afdw*self.afdw_lignin
        
        treated.imass['H2O'] = ww.F_mass - solid.F_mass
    
    @property
    def moisture_level(self):
        return self.moisture
    
    @property
    def dw_protein(self):
        return self.afdw_protein*(1-self.dw_ash)
    
    @property
    def dw_lipid(self):
        return self.afdw_lipid*(1-self.dw_ash) 
    
    @property
    def dw_carbo(self):
        return self.afdw_carbo*(1-self.dw_ash)
    
    @property
    def dw_lignin(self):
        return self.afdw_lignin*(1-self.dw_ash)
    
    # TODO: ask Jianan if lignin should be included in C_ratio
    # TODO: yes, we need to determine the element composition of lignin
    # maybe the same as carbohydrate
    # TODO: add unceratinty for lignin_2_C
    @property
    def C_ratio(self):
       return self.dw_protein*self.protein_2_C + self.dw_lipid*self.lipid_2_C +\
           self.dw_carbo*self.carbo_2_C
    
    # TODO: ask Jianan if lignin should be included in H_ratio       
    # TODO: yes, we need to determine the element composition of lignin  
    # maybe the same as carbohydrate
    # TODO: add unceratinty for lignin_2_H
    @property
    def H_ratio(self):
       return self.dw_protein*self.protein_2_H + self.dw_lipid*self.lipid_2_H +\
           self.dw_carbo*self.carbo_2_H
    
    @property
    def N_ratio(self):
       return self.dw_protein*self.protein_2_N
    
    @property
    def P_ratio(self):
       return self.N_ratio*self.N_2_P
    
    @property
    def O_ratio(self):
       return 1 - self.C_ratio - self.H_ratio -\
           self.N_ratio - self.dw_ash
    
    @property
    def C(self):
       return self.C_ratio*self.dw
    
    @property
    def H(self):
       return self.H_ratio*self.dw
    
    @property
    def N(self):
       return self.N_ratio*self.dw
    
    @property
    def P(self):
       return self.P_ratio*self.dw
    
    @property
    def O(self):
       return self.O_ratio*self.dw
    
    @property
    def AOSc(self):
       return (3*self.N_ratio/14.0067 + 2*self.O_ratio/15.999 -\
               self.H_ratio/1.00784)/(self.C_ratio/12.011)
    
    @property
    def HHV(self):
       return 100*(0.338*self.C_ratio + 1.428*(self.H_ratio -\
              self.O_ratio/8)) # [1]
    
    @property
    def H_C_eff(self):
        return (self.H/1.00784-2*self.O/15.999)/self.C*12.011
    
# =============================================================================
# WSConverter
# =============================================================================

class preWSConverter(SanUnit):
    '''
    A fake unit that enables LCA. This unit must be added before the converted stream.
    
    Parameters
    ----------
    ins : Iterable(stream)
        inlet
    outs : Iterable(stream)
        outlet
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
    def _run(self):
        
        inlet = self.ins[0]
        outlet = self.outs[0]
        
        inlet.copy_like(outlet)