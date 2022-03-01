# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the BioSTEAM package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit
from ..equipments import VerticalMixer, VacuumPump, HollowFiberMembrane, Beads
from ..utils import auom
from warnings import warn
from math import pi

__all__ = ('H2E', 'CH4E')

class H2E(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 1
    A = 265
    b = 0.513
    _Vmin = 1e4    # gal
    _Vmax = 1e6    # gal    
    _rho_bead = 265  # kg/m3
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 init_with='WasteStream', equipments=(), lifetime=20, 
                 COD_removal=0.2, H2_yield=4e-4, CH4_yield=0.0, frac_H2=0.28, 
                 tau=1, T=273.15+35, safety_factor=1.3, 
                 V_frac_beads=0.09, e_heat=0.8, p_treatment=0.24, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, 
                         equpiments=equipments, lifetime=lifetime, F_BM_default=1)
        self.COD_removal = COD_removal
        self.H2_yield = H2_yield
        self.CH4_yield = CH4_yield
        self.frac_H2 = frac_H2
        self.tau = tau
        self.T = T
        self.safety_factor = safety_factor
        self.V_frac_beads = V_frac_beads
        self.e_heat = e_heat
        self.p_treatment = p_treatment
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_equip(lifetime)

    def _init_equip(self, lifetime):
        if self.equipments:
            isa = isinstance
            equips = list(self.equipments)
            vm = None
            vp = None
            mb = None
            encap = None
            for i in equips:
                if isa(i, VerticalMixer): 
                    vm = i
                    equips.remove(i)
                elif isa(i, VacuumPump): 
                    vp = i
                    equips.remove(i)
                elif isa(i, HollowFiberMembrane):
                    mb = i
                    equips.remove(i)
                elif isa(i, Beads):
                    encap = i
                    equips.remove(i)
            if not vm:
                warn('lacking vertical mixer as equipment, will use a default one')
                vm = VerticalMixer(linked_unit=self, ID='mixer', lifetime=lifetime)
            if not vp:
                warn('lacking vacuum pump as equipment, will use a default one')
                vp = VacuumPump(linked_unit=self, ID='vacuum', lifetime=lifetime)
            if not mb:
                warn('lacking hollow fiber membrane as equipment, will use a default one')
                mb = HollowFiberMembrane(linked_unit=self, ID='membrane', lifetime=10)
            if not encap:
                warn('lacking encapsulation beads as equipment, will use a default one')
                encap = Beads(linked_unit=self, ID='beads', lifetime=0.5)
            self.equipments = (vm, vp, mb, encap, *equips)
        else:
            vm = VerticalMixer(linked_unit=self, ID='mixer', lifetime=lifetime)
            vp = VacuumPump(linked_unit=self, ID='vacuum', lifetime=lifetime)
            mb = HollowFiberMembrane(linked_unit=self, ID='membrane', lifetime=10)
            encap = Beads(linked_unit=self, ID='beads', lifetime=0.5)
            self.equipments = (vm, vp, mb, encap)

    def _run(self):
        waste, = self.ins
        eff, biogas = self.outs
        eff.copy_like(waste)
        biogas.phase = 'g'

        # COD removal
        COD_rmv = eff.imass['COD'] * self.COD_removal
        eff.imass['COD'] -= COD_rmv
        biogas.imass['H2'] = h2 =  COD_rmv*self.H2_yield
        biogas.imass['CH4'] = ch4 = COD_rmv*self.CH4_yield
        biogas.imass['N2'] = h2/self.frac_H2 - h2 - ch4

    _units = {
        'Residence time': 'd',
        'Reactor volume': 'm3',
        'Bead total volume': 'm3'
        }

    def _design(self):
        Q = self.ins[0].F_vol * 24 #m3/d
        design = self.design_results
        design['Residence time'] = t = self.tau
        design['Reactor volume'] = V = t * Q * self.safety_factor
        design['Bead total volume'] = V * self.V_frac_beads
        self.add_equipment_design()
        mixer, vacuum = self.equipments[:2]
        P_mix = mixer.power * mixer.N_mix
        P_vcm = vacuum.P_vacuum
        self.power_utility(rate=P_mix+P_vcm)
        hu = self.heat_utilities[0]
        hu.heat_transfer_efficiency = self.e_heat
        inf = self.ins[0]
        T_in = inf.T
        T_target = self.T
        unit_duty = inf.F_mass * inf.Cp * (T_target - T_in) #kJ/hr
        hu(unit_duty, T_in, T_target)
        
    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        V = D['Reactor volume'] * auom('m3').conversion_factor('gal')
        C['Reactor'] = self.A * min(self._Vmax, max(self._Vmin, V)) ** self.b   # cone-roof carbon steel storage tank
        self.add_equipment_cost()
        eff = self.outs[0]
        if eff.isproduct():
            self.add_OPEX = {'Effluent treatment': self.p_treatment * eff.imass['COD']} # USD/hr

    @property
    def COD_removal(self):
        return self._rcod
    
    @COD_removal.setter
    def COD_removal(self, r):
        if r > 1 or r < 0:
            raise ValueError(f'COD removal must be in [0, 1], not {r}')
        self._rcod = r
    
    @property
    def H2_yield(self):
        return self._yh2
    
    @H2_yield.setter
    def H2_yield(self, y):
        if y > 8 or y < 0:
            raise ValueError(f'H2 yield must be in [0, 8] g-H2/g-COD-removed, not {y}')
        self._yh2 = y
        
    @property
    def CH4_yield(self):
        return self._ych4
    
    @CH4_yield.setter
    def CH4_yield(self, y):
        if y > 4 or y < 0:
            raise ValueError(f'H2 yield must be in [0, 4] g-CH4/g-COD-removed, not {y}')
        self._ych4 = y
    
    @property
    def frac_H2(self):
        return self._fh2
    @frac_H2.setter
    def frac_H2(self, f):
        if f > 1 or f < 0:
            raise ValueError(f'H2 mass fraction in biogas must be in [0, 1], not {f}')
        self._fh2 = f
    
    @property
    def tau(self):
        return self._tau
    
    @tau.setter
    def tau(self, t):
        if t < 0: raise ValueError(f'residence time tau cannot be negative: {t}')
        self._tau = t

    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, i):
        if i < self.ins[0].T:
            warn('Operating temperature should not be lower than influent temperature {self.ins[0].T}')
            i = self.ins[0].T
        self._T = i

    @property
    def safety_factor(self):
        return self._sf
    @safety_factor.setter
    def safety_factor(self, sf):
        if sf < 1: raise ValueError(f'safety factor cannot be less than 1: {sf}')
        self._sf = sf

    @property
    def V_frac_beads(self):
        return self._f_Vb
    @V_frac_beads.setter
    def V_frac_beads(self, f):
        if f > 1 or f < 0:
            raise ValueError(f'Volume fraction of beads in the reactor must be in [0, 1], not {f}')
        self._f_Vb = f
    
    @property
    def e_heat(self):
        return self._eh
    @e_heat.setter
    def e_heat(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'Heat transfer efficiency must be in [0, 1], not {e}')
        self._eh = e

    @property
    def p_treatment(self):
        return self._ptreat     
    @p_treatment.setter
    def p_treatment(self, p):
        self._ptreat = p # USD/kg-COD

class CH4E(H2E):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 equipments=(), lifetime=20, COD_removal=0.55, CH4_yield=0.135, 
                 frac_CH4=0.62, tau=15, T=273.15+35, d_wall=0.1, d_slab=0.1, 
                 safety_factor=1.3, V_frac_beads=0.09, e_heat=0.8, 
                 p_treatment=0.24, p_wall_concrete=497.25, 
                 p_slab_concrete=267.75, p_steel=1203, **kwargs):        
        H2E.__init__(self, ID, ins, outs, thermo, init_with, equipments, lifetime,
                     COD_removal, 0, CH4_yield, 0, tau, T, safety_factor, 
                     V_frac_beads, e_heat, p_treatment, **kwargs)
        self.frac_CH4 = frac_CH4
        self.d_wall = d_wall
        self.d_slab = d_slab
        self.p_wall_concrete = p_wall_concrete
        self.p_slab_concrete = p_slab_concrete
        self.p_steel = p_steel

    def _run(self):
        waste, = self.ins
        eff, biogas = self.outs
        eff.copy_like(waste)
        biogas.phase = 'g'

        # COD removal
        COD_rmv = eff.imass['COD'] * self.COD_removal
        eff.imass['COD'] -= COD_rmv
        biogas.imass['H2'] = h2 =  COD_rmv*self.H2_yield
        biogas.imass['CH4'] = ch4 = COD_rmv*self.CH4_yield
        biogas.imass['N2'] = ch4/self.frac_CH4 - h2 - ch4

    _units = {
        'Wall concrete volume': 'm3',
        'Slab concrete volume': 'm3',
        'Steel cover area': 'm2',
        **H2E._units
        }

    def _design(self):
        H2E._design(self)
        design = self.design_results
        V = design['Reactor volume']
        design['Wall concrete volume'] = 2 * pi**(1/3) * V**(2/3) * self.d_wall
        design['Slab concrete volume'] = pi**(1/3) * V**(2/3) * self.d_slab
        design['Steel cover area'] = 5**(1/2) / 2 * pi**(1/3) * V**(2/3)

    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete volume'] * self.p_wall_concrete
        C['Slab concrete'] = D['Slab concrete volume'] * self.p_slab_concrete
        C['Steel cover'] = D['Steel cover area'] * self.p_steel
        self.add_equipment_cost()
        eff = self.outs[0]
        if eff.isproduct():
            self.add_OPEX = {'Effluent treatment': self.p_treatment * eff.imass['COD']} # USD/hr

    @property
    def frac_CH4(self):
        return self._fch4
    @frac_CH4.setter
    def frac_CH4(self, f):
        if f > 1 or f < 0:
            raise ValueError(f'CH4 mass fraction in biogas must be in [0, 1], not {f}')
        self._fch4 = f
    
    @property
    def d_wall(self):
        return self._dw

    @d_wall.setter
    def d_wall(self, d):
        if d < 0: raise ValueError(f'wall thickness cannot be negative: {d}')
        self._dw = d
    
    @property
    def d_slab(self):
        return self._ds

    @d_slab.setter
    def d_slab(self, d):
        if d < 0: raise ValueError(f'slab thickness cannot be negative: {d}')
        self._ds = d    
    
    @property
    def p_wall_concrete(self):
        return self._pwall    
    @p_wall_concrete.setter
    def p_wall_concrete(self, p):
        self._pwall = p # USD/m3
    
    @property
    def p_slab_concrete(self):
        return self._pslab    
    @p_slab_concrete.setter
    def p_slab_concrete(self, p):
        self._pslab = p # USD/m3
    
    @property
    def p_steel(self):
        return self._pstl   
    @p_steel.setter
    def p_steel(self, p):
        self._pstl = p # USD/m2