# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Component, Components, WasteStream, set_thermo, sanunits as su

__all__ = ('cmps',)


# %%

# =============================================================================
# Components and streams
# =============================================================================

X = Component('X', phase='s', measured_as='COD', i_COD=0, description='Biomass',
              organic=True, particle_size='Particulate', degradability='Readily')
X_inert = Component('X_inert', phase='s', description='Inert biomass', i_COD=0,
                    organic=True, particle_size='Particulate', degradability='Undegradable')
Substrate = Component('Substrate', phase='s', i_COD=300/18.3,
                      organic=True, particle_size='Particulate', degradability='Readily')
cmps = Components([X, X_inert, Substrate])

CH4 = Component('CH4', phase='g', organic=True, particle_size='Particulate',
                degradability='Readily')
cmps = Components.append_combustion_components([*cmps, CH4], lock_state_at='')
set_thermo(cmps)


inf = WasteStream('inf')
inf.set_flow_by_concentration(flow_tot=20,
                              concentrations={'Substrate': 18.3},
                              units=('mgd', 'mg/L'))

'''
System should have the following units:
    Screening
    ActivatedSludgeProcess (aeration tanks & clarifiers)
    Sludge handling units (at least GBT, maybe also others):
        X_thickened_WAS = 25000; %[mg/L]
    CHP

'''

# @property
# def active_biomass(self):
#     '''[tuple] IDs of active biomass (particulate & organic & degradable) components.'''
#     return self.get_IDs_from_array(self.x*self.b*self.org)

# @property
# def inert_biomass(self):
#     '''[tuple] IDs of inert biomass (particulate & organic & undegradable) components.'''
#     return self.get_IDs_from_array(self.x*self.org-self.x*self.b*self.org)