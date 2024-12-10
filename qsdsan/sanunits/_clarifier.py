# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from numpy import maximum as npmax, minimum as npmin, exp as npexp
from math import sqrt, pi
from warnings import warn
from numba import njit
from .. import SanUnit, WasteStream
import numpy as np
from ..sanunits import WWTpump
from ..sanunits._pumping import default_F_BM as default_WWTpump_F_BM
from ..sanunits import dydt_cstr

__all__ = ('FlatBottomCircularClarifier',
           'IdealClarifier',
           'PrimaryClarifierBSM2',
           'PrimaryClarifier')

F_BM_pump = 1.18*(1 + 0.007/100) # 0.007 is for miscellaneous costs

default_F_BM = {
        'Pumps': F_BM_pump,
        'Pump building': F_BM_pump,
        }

default_equipment_lifetime = {
    'Pumps': 15,
    'Pump pipe stainless steel': 15,
    'Pump stainless steel': 15,
    }

# Assign a bare module of 1 to all
default_F_BM = {
        'Wall concrete': 1.,
        'Slab concrete': 1.,
        'Wall stainless steel': 1.,
        'Scraper': 1,
        'v notch weir': 1,
        'Pumps': 1
        }
default_F_BM.update(default_WWTpump_F_BM)

#%% Takács Clarifer
@njit(cache=True)
def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp, n0):
    X_star = npmax(X-X_min, n0)
    v = npmin(v_max_practical, v_max*(npexp(-rh*X_star) - npexp(-rp*X_star)))
    return X*npmax(v, n0)

# from math import exp
# def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp, n0):
#     X_star = max(X-X_min, 0)
#     v = min(v_max_practical, v_max*(exp(-rh*X_star) - exp(-rp*X_star)))
#     return X*max(v, 0)


class FlatBottomCircularClarifier(SanUnit):
    """
    A flat-bottom circular clarifier with a simple 1-dimensional
    N-layer settling model.

    Parameters
    ----------
    ID : str
        ID for the clarifier. The default is ''.
    ins : :class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1.
    outs : :class:`WasteStream`
        Treated effluent and sludge.
    underflow : float, optional
        Designed recycling sludge flowrate (RAS), in [m^3/d]. The default is 2000.
    wastage : float, optional
        Designed wasted sludge flowrate (WAS), in [m^3/d]. The default is 385.
    surface_area : float, optional
        Surface area of the clarifier, in [m^2]. The default is 1500.
    height : float, optional
        Height of the clarifier, in [m]. The default is 4.
    N_layer : int, optional
        The number of layers to model settling. The default is 10.
    feed_layer : int, optional
        The feed layer counting from top to bottom. The default is 4.
    X_threshold : float, optional
        Threshold suspended solid concentration, in [g/m^3]. The default is 3000.
    v_max : float, optional
        Maximum theoretical (i.e. Vesilind) settling velocity, in [m/d]. The
        default is 474.
    v_max_practical : float, optional
        Maximum practical settling velocity, in [m/d]. The default is 250.
    rh : float, optional
        Hindered zone settling parameter in the double-exponential settling velocity
        function, in [m^3/g]. The default is 5.76e-4.
    rp : float, optional
        Flocculant zone settling parameter in the double-exponential settling velocity
        function, in [m^3/g]. The default is 2.86e-3.
    fns : float, optional
        Non-settleable fraction of the suspended solids, dimensionless. Must be within
        [0, 1]. The default is 2.28e-3.
    maximum_nonsettleable_solids : float, optional
        Maximum non-settleable solids concentration, in mgTSS/L. The default is None.
    downward_flow_velocity : float, optional
        Speed on the basis of which center feed diameter is designed [m/hr]. The default is 42 m/hr (0.7 m/min). [2]
    design_influent_TSS : float, optional
        The design TSS concentration [mg/L] in the influent going to the secondary clarifier. 
    design_influent_flow : float, optional
        The design influent tptal volumetric flow [m3/hr] going to the secondary clarifier. 
    design_solids_loading_rate : float, optional
        Rate of total suspended solids entering the secondary clarifier (kg/(m2*hr)). 
        The default is 5 kg/(m2*hr) [3, 4]
    
    References
    ----------
    [1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the Clarification
    -Thickening Process. Water Res. 1991, 25 (10), 1263–1271.
    https://doi.org/10.1016/0043-1354(91)90066-Y.
    
    [2] Chapter-12: Suspended-growth Treatment Processes. WEF Manual of Practice No. 8. 
    6th Edition. Virginia: McGraw-Hill, 2018. 
    
    [3] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    
    [4] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    
    [5] RECOMMENDED STANDARDS for WASTEWATER FACILITIES. 10 state standards. 2014 edition. 
    """

    _N_ins = 1
    _N_outs = 3
    
    # # Costs
    # wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    # slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    # stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    # pumps = ('ras', 'was',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', underflow=2000, wastage=385,
                 surface_area=1500, height=4, N_layer=10, feed_layer=4,
                 X_threshold=3000, v_max=474, v_max_practical=250,
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3, 
                 maximum_nonsettleable_solids=None,
                 F_BM_default=default_F_BM, isdynamic=True,
                 downward_flow_velocity=42, design_influent_TSS = None, design_influent_flow = None,
                 design_solids_loading_rate = 6, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic, F_BM_default=1)
        self._h = height
        self._Qras = underflow
        self._Qwas = wastage
        self._sludge = WasteStream()
        
        if surface_area != None:
            self._A = surface_area
        elif design_influent_TSS != None and design_influent_flow != None:
            self._A = (design_influent_TSS*design_influent_flow)/(design_solids_loading_rate*1000) # 1000 in denominator for unit conversion
        else:
            RuntimeError('Either surface_area, or design_influent_TSS and design_influent_flow expected from user')
        
        self._V = self._A * height
        self._hj = height/N_layer
        self._N_layer = N_layer
        self.feed_layer = feed_layer
        self._v_max = v_max
        self._v_max_p = v_max_practical
        self._X_t = X_threshold
        self._rh = rh
        self._rp = rp
        self._fns = fns
        self.maximum_nonsettleable_solids = maximum_nonsettleable_solids
        self._solids = None
        self._solubles = None
        self._X_comp = np.zeros(len(self.components))
        self._dX_comp = self._X_comp.copy()
        
        self._downward_flow_velocity = downward_flow_velocity # in m/hr (converted from 12 mm/sec)
        self._design_tss = design_influent_TSS
        self._design_flow = design_influent_flow
        self._slr = design_solids_loading_rate
        
        self._mixed = WasteStream(f'{ID}_mixed')
        header = self._state_header
        self._state_header = list(header) + [f'TSS{i+1} [mg/L]' for i in range(N_layer)]
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
        self._inf = self.ins[0].copy(f'{ID}_inf')
        self._ras = self.outs[1].copy(f'{ID}_ras')
        self._was = self.outs[2].copy(f'{ID}_was')
        
    @property
    def height(self):
        '''[float] Height of the clarifier in m.'''
        return self._h

    @height.setter
    def height(self, h):
        self._h = h

    @property
    def underflow(self):
        '''[float] The designed recycling sludge flow rate in m3/d.'''
        return self._Qras

    @underflow.setter
    def underflow(self, ras):
        self._Qras = ras

    @property
    def wastage(self):
        '''[float] The designed wasted sludge flow rate in m3/d.'''
        return self._Qwas

    @wastage.setter
    def wastage(self, was):
        self._Qwas = was

    @property
    def V_settle(self):
        '''[float] Total volume modeled for settling in m^3, calculated based on surface area and height.'''
        return self._V

    @property
    def A_settle(self):
        '''[float] The surface area for settling in m^2, i.e., the area of the clarifier's flat bottom.'''
        return self._A

    @A_settle.setter
    def A_settle(self, A):
        self._A = A

    @property
    def h_layer(self):
        '''[float] The height of each layer in the settling model, in m.'''
        return self._hj

    @h_layer.setter
    def h_layer(self, h):
        self._hj = h

    @property
    def N_layer(self):
        '''[int] The number of layers into which the clarifier is divided in the settling model.'''
        return self._N_layer

    @N_layer.setter
    def N_layer(self, N):
        self._N_layer = int(N)

    @property
    def feed_layer(self):
        '''[int] The feed layer counting from top to bottom.'''
        return self._feed_layer

    @feed_layer.setter
    def feed_layer(self, jf):
        jf = int(jf)
        if jf > self._N_layer or jf < 1:
            raise ValueError(f'feed layer {self._feed_layer} is out of range.'
                             f'must be an integer between 1 and {self._N_layer}.')
        self._feed_layer = jf

    @property
    def v_max(self):
        '''[float] Maximum theoretical (i.e. Vesilind) settling velocity, in m/d'''
        return self._v_max

    @v_max.setter
    def v_max(self, vm):
        if vm < self._v_max_p:
            raise ValueError('v_max must be greater than or equal to v_max_p.')
        self._v_max = vm

    @property
    def v_max_p(self):
        '''[float] Maximum practical settling velocity, in m/d'''
        return self._v_max_p

    @v_max_p.setter
    def v_max_p(self, vmp):
        if vmp > self._v_max or vmp <= 0:
            raise ValueError('v_max_p must be within (0, v_max].')
        self._v_max_p = vmp

    @property
    def X_t(self):
        '''[float] Threshold suspended solid concentration, in g/m^3.'''
        return self._X_t

    @X_t.setter
    def X_t(self, xt):
        if xt < 0: raise ValueError('X_t must be positive.')
        self._X_t = xt

    @property
    def rh(self):
        '''[float] Hindered zone settling parameter in the double-exponential settling velocity function, in m^3/g.'''
        return self._rh

    @rh.setter
    def rh(self, rh):
        if rh > self._rp: raise ValueError('rh must be less than or equal to rp.')
        self._rh = rh

    @property
    def rp(self):
        '''[float] Flocculant zone settling parameter in the double-exponential settling velocity function, in m^3/g.'''
        return self._rp

    @rp.setter
    def rp(self, rp):
        if rp < self._rh: raise ValueError('rp must be greater than or equal to rp.')
        self._rp = rp

    @property
    def fns(self):
        '''[float] Non-settleable fraction of the suspended solids'''
        return self._fns

    @fns.setter
    def fns(self, fns):
        if fns < 0 or fns > 1: raise ValueError('fns must be within [0,1].')
        self._fns = fns
    
    @property
    def maximum_nonsettleable_solids(self):
        '''[float] Maximum non-settleable solids concentration, in mgTSS/L.'''
        return self._max_ns
    @maximum_nonsettleable_solids.setter
    def maximum_nonsettleable_solids(self, ns):
        self._max_ns = ns
    
    @property
    def solids_loading_rate(self):
        '''solids_loading_rate is the loading in the clarifier'''
        return self._slr
        
    @solids_loading_rate.setter
    def solids_loading_rate(self, slr):
        if slr is not None:
            self._slr = slr
        else: 
            raise ValueError('solids_loading_rate of the clarifier expected from user')
            
    def set_init_solubles(self, **kwargs):
        '''set the initial concentrations [mg/L] of solubles in the clarifier.'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        x = self.components.x
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        self._solubles = Cs*(1-x)

    def set_init_sludge_solids(self, **kwargs):
        '''set the initial concentrations [mg/L] of solids in the underflow sludge.'''
        cmpx = self.components.index
        x = self.components.x
        Xs = np.zeros_like(x)
        for k, v in kwargs.items(): Xs[cmpx(k)] = v
        Xs *= x
        tss = sum(Xs * self.components.i_mass)
        if tss != 0: self._X_comp = Xs / tss

    def set_init_TSS(self, arr):
        '''set the initial TSS [mg/L] in each layer of the clarifier.'''
        if len(arr) != self._N_layer:
            raise ValueError(f'expects an iterable of length {self._N_layer}, not {len(arr)}')
        self._solids = np.asarray(arr, dtype=float)

    def _init_state(self):
        n = self._N_layer
        x = self.components.x
        imass = self.components.i_mass
        QCs = self._ins_QC[0]
        Q = QCs[-1]
        Z = self._solubles if self._solubles is not None \
            else QCs[:-1]*(1-x)
        TSS_in = sum(QCs[:-1] * x * imass)
        TSS = self._solids if self._solids is not None \
            else TSS_in*(20**np.linspace(-1,1,n))
        ZQs = np.append(Z, Q)
        self._state = np.append(ZQs, TSS)
        self._dstate = self._state * 0.
        if TSS_in != 0: self._X_comp = QCs[:-1] * x / TSS_in

    def _update_state(self):
        arr = self._state
        x = self.components.x
        n = self._N_layer
        arr[-(1+n)] = Q_in = self._ins_QC[0, -1]
        Q_e = Q_in - self._Qras - self._Qwas
        # Q_e = arr[-(1+n)] - self._Qras - self._Qwas
        Z = arr[:len(x)]
        inf, = self.ins
        imass = self.components.i_mass
        C_in = inf.state[:-1]
        X_composition = self._X_comp = C_in*x/sum(C_in*imass*x)
        X_e = arr[-n] * X_composition
        C_s = Z + arr[-1] * X_composition
        eff, ras, was = self._outs
        if eff.isproduct() and eff.state is None:
            eff.state = np.append(Z+X_e, Q_e)
        else:
            eff.state[:-1] = Z+X_e   # not sure if this works for a setter
            eff.state[-1] = Q_e
        #!!! might need to enable dynamic sludge volume flows
        if ras.isproduct() and ras.state is None:
            ras.state = np.append(C_s, self._Qras)
        else:
            ras.state[:-1] = C_s
            ras.state[-1] = self._Qras
        if was.isproduct() and was.state is None:
            was.state = np.append(C_s, self._Qwas)
        else:
            was.state[:-1] = C_s
            was.state[-1] = self._Qwas

    def _update_dstate(self):
        arr = self._dstate
        x = self.components.x
        n = self._N_layer
        dQ = arr[-(1+n)]
        dZ = arr[:len(x)]
        TSS_e, TSS_s = self._state[-n], self._state[-1]
        X_composition = self._X_comp # (m, ), mg COD/ mg TSS
        dX_composition = self._dX_comp
        dC_e = dZ + arr[-n] * X_composition + dX_composition * TSS_e
        dC_s = dZ + arr[-1] * X_composition + dX_composition * TSS_s
        eff, ras, was = self._outs
        if eff.isproduct() and eff.dstate is None:
            eff.dstate = np.append(dC_e, dQ)
        else:
            eff.dstate[:-1] = dC_e # not sure if this works for a setter
            eff.dstate[-1] = dQ
        #!!! might need to enable dynamic sludge volume flows
        if ras.isproduct() and ras.dstate is None:
            ras.dstate = np.append(dC_s, 0.)
        else:
            ras.dstate[:-1] = dC_s
        if was.isproduct() and was.dstate is None:
            was.dstate = np.append(dC_s, 0.)
        else:
            was.dstate[:-1] = dC_s

    def _run(self):
        '''only to converge volumetric flows.'''
        inf, = self.ins
        sludge = self._sludge
        Q_in = inf.get_total_flow('m3/d')
        eff, ras, was = self.outs
        Q_ras = self._Qras
        Q_was = self._Qwas
        s_e = 1 - (Q_ras+Q_was)/Q_in
        inf.split_to(eff, sludge, s_e)
        sludge.split_to(ras, was, Q_ras/(Q_ras+Q_was))

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        tss = self._state[-self._N_layer:].mean()
        mass = cmps.i_mass * self._X_comp * tss
        return self._V * mass[cmps.indices(biomass_IDs)].sum()

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        n = self._N_layer
        jf = self._feed_layer - 1
        x = self.components.x
        m = len(x)
        imass = self.components.i_mass
        fns = self._fns
        max_ns = self._max_ns or 1e6
        Q_s = self._Qras + self._Qwas

        dQC = self._dstate
        dX_comp = self._dX_comp
        _update_dstate = self._update_dstate

        nzeros = np.zeros(n)
        Q_jout = nzeros.copy()
        X_rolled = nzeros.copy()
        X_min_arr = nzeros.copy()
        settle_out = nzeros.copy()
        settle_in = nzeros.copy()

        # Make these constants into arrays so it'll be faster in `dy_dt`
        # vmax_arr = np.full_like(nzeros, self._v_max)
        # vmaxp_arr = np.full_like(nzeros, self._v_max_p)
        # rh_arr = np.full_like(nzeros, self._rh)
        # rp_arr = np.full_like(nzeros, self._rp)
        vmax_arr = self._v_max
        vmaxp_arr = self._v_max_p
        rh_arr = self._rh
        rp_arr = self._rp
        func_vx = lambda x_arr, xmin_arr : _settling_flux(x_arr, vmax_arr, vmaxp_arr, xmin_arr, rh_arr, rp_arr, nzeros)
        
        A, hj, V = self._A, self._hj, self._V
        # A_arr = np.full_like(nzeros, A)
        # hj_arr = np.full_like(nzeros, hj)
        J = np.zeros(n-1)
        # X_t_arr = np.full(jf, self._X_t)
        X_t = self._X_t
        # Q_in_arr = np.zeros(m)
        # V_arr = np.full(m, V)

        def dy_dt(t, QC_ins, QC, dQC_ins):
            # dQC[-(n+1)] = dQC_ins[0,-1]
            dQC[-(n+1)] = 0.
            Q_in = QC_ins[0,-1]
            Q_e = Q_in - Q_s
            C_in = QC_ins[0,:-1]
            dC_in = dQC_ins[0,:-1]
            Z_in = C_in*(1-x)
            X_in = sum(C_in*imass*x)           # influent TSS
            dX_in = sum(dC_in*imass*x)
            X_min_arr[:] = min(X_in * fns, max_ns)
            X = QC[-n:]                        # (n, ), TSS for each layer
            Z = QC[:m] * (1-x)
            #***********TSS*************
            Q_jout[:jf] = Q_e
            Q_jout[jf] = Q_in
            Q_jout[jf+1:] = Q_s
            flow_out = X * Q_jout
            X_rolled[:jf] = X[1: jf+1]
            X_rolled[jf] = X_in
            X_rolled[jf+1:] = X[jf: -1]
            flow_in = X_rolled * Q_jout
            VX = func_vx(X, X_min_arr)
            J[:] = npmin(VX[:-1], VX[1:])
            condition = (X_rolled[:jf]<X_t)
            J[:jf][condition] = VX[:jf][condition]
            settle_out[:-1] = J
            settle_in[1:] = J
            dQC[-n:] = ((flow_in - flow_out)/A + settle_in - settle_out)/hj       # (n,)
            #*********solubles**********
            dQC[:m] = Q_in/V*(Z_in - Z)
            # instrumental variables
            dX_comp[:] = (dC_in * X_in - dX_in * C_in) * x / X_in**2
            _update_dstate()

        self._ODE = dy_dt
    
    #!!! should consolidate design & costing equations between primary & secondary clarifiers
    # _units = {
    #     'Number of clarifiers': 'ea',
    #     'Volumetric flow': 'm3/day',
    #     'Clarifier depth': 'm',
    #     'Surface area': 'm2',
    #     'Clarifier diameter': 'm',
    #     'Clarifier volume': 'm3',
    #     'Design solids loading rate': 'kg/m2/hr',
    #     'Surface overflow rate': 'm3/day/m2',
    #     'Hydraulic Retention Time': 'hr', 
    #     'Center feed depth': 'm',
    #     'Downward flow velocity': 'm/hr',
    #     'Center feed diameter': 'm',
    #     'Volume of concrete wall': 'm3',
    #     'Stainless steel': 'kg',
    #     'Pump pipe stainless steel' : 'kg',
    #     'Pump stainless steel': 'kg',
    #     'Number of pumps': 'ea'
    # }
     
    # def _design_pump(self):
    #     ID, pumps = self.ID, self.pumps
    
    #     self._ras.copy_like(self.outs[1])
    #     self._was.copy_like(self.outs[2])
        
    #     ins_dct = {
    #         'ras': self._ras,
    #         'was': self._was,
    #         }
        
    #     D = self.design_results
        
    #     ras_flow = self._ras.get_total_flow('m3/hr')
    #     was_flow = self._was.get_total_flow('m3/hr')
        
    #     ras_flow_u = ras_flow/D['Number of clarifiers']*0.00634
    #     was_flow_u = was_flow/D['Number of clarifiers']*0.00634
        
    #     Q_mgd = {
    #         'ras': ras_flow_u,
    #         'was': was_flow_u,
    #         }
        
    #     type_dct = dict.fromkeys(pumps, 'sludge')
    #     inputs_dct = dict.fromkeys(pumps, (1,))
       
    #     for i in pumps:
    #         if hasattr(self, f'{i}_pump'):
    #             p = getattr(self, f'{i}_pump')
    #             setattr(p, 'add_inputs', inputs_dct[i])
    #         else:
    #             ID = f'{ID}_{i}'
    #             capacity_factor=1
    #             pump = WWTpump(
    #                 ID=ID, ins=ins_dct[i], thermo = self.thermo, pump_type=type_dct[i],
    #                 Q_mgd=Q_mgd[i], add_inputs=inputs_dct[i],
    #                 capacity_factor=capacity_factor,
    #                 include_pump_cost=True,
    #                 include_building_cost=False,
    #                 include_OM_cost=True,
    #                 )
    #             setattr(self, f'{i}_pump', pump)

    #     pipe_ss, pump_ss = 0., 0.
    #     for i in pumps:
    #         p = getattr(self, f'{i}_pump')
    #         p.simulate()
    #         p_design = p.design_results
    #         pipe_ss += p_design['Pump pipe stainless steel']
    #         pump_ss += p_design['Pump stainless steel']
    #     return pipe_ss, pump_ss
     
    # def _design(self):
        
    #     self._mixed.mix_from(self.ins)
    #     mixed = self._mixed
    #     D = self.design_results
        
    #     # Number of clarifiers based on tentative suggestions by Jeremy 
    #     # (would be verified through collaboration with industry)
    #     total_flow = (mixed.get_total_flow('m3/hr')*24)/3785 # in MGD
    #     if total_flow <= 3:
    #         D['Number of clarifiers'] = 2
    #     elif total_flow > 3 and total_flow <= 8:
    #         D['Number of clarifiers'] = 3
    #     elif total_flow > 8 and total_flow <=20:
    #         D['Number of clarifiers'] = 4
    #     else:
    #         D['Number of clarifiers'] = 4
    #         total_flow -= 20
    #         D['Number of clarifiers'] += np.ceil(total_flow/20)
                
    #     D['Volumetric flow'] =  (mixed.get_total_flow('m3/hr')*24)/D['Number of clarifiers'] #m3/day
        
    #     # Sidewater depth of a cylindrical clarifier lies between 4-5 m (MOP 8)
    #     D['Clarifier depth'] = self._h # in m
        
    #     # Area of clarifier 
    #     # D['Surface area'] = solids_clarifier/D['Solids loading rate'] #m2
    #     D['Surface area'] = self._A/D['Number of clarifiers']
    #     D['Clarifier diameter'] = np.sqrt(4*D['Surface area']/np.pi) # in m
    #     D['Clarifier volume'] = D['Surface area']*D['Clarifier depth'] # in m3
        
    #     # Checks on SLR,, SOR, and HRT 
        
    #     D['Design solids loading rate'] = self._slr # kg/(m2*hr)
        
    #     total_solids = mixed.get_TSS()*mixed.get_total_flow('m3/hr')/1000 # in kg/hr (mg/l * m3/hr)
    #     solids_clarifier = total_solids/D['Number of clarifiers'] # in kg/hr
    #     simulated_slr = solids_clarifier/D['Surface area'] # in kg/(m2*hr)
        
    #     # Consult Joy on the margin or error
    #     if simulated_slr < 0.8*D['Design solids loading rate'] or simulated_slr > 1.2*D['Design solids loading rate']:
    #         design_slr = D['Design solids loading rate']
    #         warn(f'Solids loading rate = {simulated_slr} is not within 20% of the recommended design level of {design_slr} kg/hr/m2')
        
    #     # Check on SLR [3, 4, 5] 
    #     if simulated_slr > 14:
    #         warn(f'Solids loading rate = {simulated_slr} is above recommended level of 14 kg/hr/m2')
        
    #     # Check on SOR [3, 4, 5]
    #     D['Surface overflow rate'] = D['Volumetric flow']/D['Surface area']  # in m3/m2/hr
    #     if D['Surface overflow rate'] > 49:
    #         sor = D['Surface overflow rate']
    #         warn(f'Surface overflow rate = {sor} is above recommended level of 49 m3/day/m2')
        
    #     # HRT
    #     D['Hydraulic Retention Time'] = D['Clarifier volume']*24/D['Volumetric flow'] # in hr
        
    #     # Clarifiers can be center feed or peripheral feed. The design here is for the more commonly deployed center feed.
    #     # Depth of the center feed lies between 30-75% of sidewater depth. [2]
    #     D['Center feed depth'] = 0.5*D['Clarifier depth']
    #     # Criteria for downward velocity of flow determine 
    #     D['Downward flow velocity'] = self._downward_flow_velocity # in m/hr
    #     Center_feed_area = (D['Volumetric flow']/24)/D['Downward flow velocity'] # in m2
    #     D['Center feed diameter'] = np.sqrt(4*Center_feed_area/np.pi)

    #     #Sanity check: Diameter of the center feed lies between 20-25% of tank diameter [2]
    #     if D['Center feed diameter'] < 0.20*D['Clarifier diameter'] or D['Center feed diameter']  > 0.25*D['Clarifier diameter']:
    #         cf_dia = D['Center feed diameter'] 
    #         tank_dia = D['Clarifier diameter']
    #         warn(f'Diameter of the center feed does not lie between 20-25% of tank diameter. It is {cf_dia*100/tank_dia} % of tank diameter')
            
    #     # Amount of concrete required
    #     D_tank = D['Clarifier depth']*39.37 # m to inches 
    #     # Thickness of the wall concrete, [m]. Default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft. (Brian's code)
    #     thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
    #     inner_diameter = D['Clarifier diameter']
    #     outer_diameter = inner_diameter + 2*thickness_concrete_wall
    #     D['Volume of concrete wall']  = (np.pi*D['Clarifier depth']/4)*(outer_diameter**2 - inner_diameter**2)
        
    #     # Concrete slab thickness, [ft], default to be 2 in thicker than the wall thickness. (Brian's code)
    #     thickness_concrete_slab = thickness_concrete_wall + (2/12)*0.3048 # from inch to m
    #     # From Brian's code
    #     D['Volume of concrete slab']  = (thickness_concrete_slab + thickness_concrete_wall)*D['Surface area']
        
    #     # Amount of metal required for center feed
    #     thickness_metal_wall = 0.3048 # equal to 1 feet, in m (!! NEED A RELIABLE SOURCE !!)
    #     inner_diameter_center_feed = D['Center feed diameter']
    #     outer_diameter_center_feed = inner_diameter_center_feed + 2*thickness_metal_wall
    #     volume_center_feed = (np.pi*D['Center feed depth']/4)*(outer_diameter_center_feed**2 - inner_diameter_center_feed **2)
    #     density_ss = 7930 # kg/m3, 18/8 Chromium
    #     D['Stainless steel'] = volume_center_feed*density_ss # in kg
       
    #     # Pumps
    #     pipe, pumps = self._design_pump()
    #     D['Pump pipe stainless steel'] = pipe
    #     D['Pump stainless steel'] = pumps
        
    #     # For secondary clarifier
    #     D['Number of pumps'] = 2*D['Number of clarifiers']
        
    # def _cost(self):
       
    #     D = self.design_results
    #     C = self.baseline_purchase_costs
       
    #     # Construction of concrete and stainless steel walls
    #     C['Wall concrete'] = D['Number of clarifiers']*D['Volume of concrete wall']*self.wall_concrete_unit_cost
        
    #     C['Slab concrete'] = D['Number of clarifiers']*D['Volume of concrete slab']*self.slab_concrete_unit_cost
        
    #     C['Wall stainless steel'] = D['Number of clarifiers']*D['Stainless steel']*self.stainless_steel_unit_cost
        
    #     # Cost of equipment 
        
    #     # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
    #     # Scraper 
    #     # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
    #     # base_cost_scraper = 2500
    #     # base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
    #     clarifier_flow =  D['Volumetric flow']/24 # in m3/hr
    #     # C['Scraper'] = D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
    #     # base_power_scraper = 2.75 # in kW
    #     # THE EQUATION BELOW IS NOT CORRECT TO SCALE SCRAPER POWER REQUIREMENTS 
    #     # scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
    #     # v notch weir
    #     # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
    #     base_cost_v_notch_weir = 6888
    #     base_flow_v_notch_weir = 10 # in m3/hr
    #     C['v notch weir'] = D['Number of clarifiers']*base_cost_v_notch_weir*(clarifier_flow/base_flow_v_notch_weir)**0.6
       
    #     # Pump (construction and maintainance)
    #     pumps = self.pumps
    #     add_OPEX = self.add_OPEX
    #     pump_cost = 0.
    #     building_cost = 0.
    #     opex_o = 0.
    #     opex_m = 0.
       
    #     # i would be 0 and 1 for RAS and WAS respectively 
    #     for i in pumps:
    #         p = getattr(self, f'{i}_pump')
    #         p_cost = p.baseline_purchase_costs
    #         p_add_opex = p.add_OPEX
    #         pump_cost += p_cost['Pump']
    #         building_cost += p_cost['Pump building']
    #         opex_o += p_add_opex['Pump operating']
    #         opex_m += p_add_opex['Pump maintenance']
            
    #     # All costs associated with pumping need to be multiplied by number of clarifiers 
    #     C['Pumps'] = pump_cost*D['Number of clarifiers']
    #     C['Pump building'] = building_cost*D['Number of clarifiers']
    #     add_OPEX['Pump operating'] = opex_o*D['Number of clarifiers']
    #     add_OPEX['Pump maintenance'] = opex_m*D['Number of clarifiers']
       
    #     # Power
    #     pumping = 0.
    #     for ID in self.pumps:
    #         p = getattr(self, f'{ID}_pump')
    #         if p is None:
    #             continue
    #         pumping += p.power_utility.rate
        
    #     pumping = pumping*D['Number of clarifiers']
        
    #     self.power_utility.rate += pumping
    #     # self.power_utility.consumption += scraper_power
        
# %% 
   
class IdealClarifier(SanUnit):    
    """
    Ideal clarifier that settles suspended solids by specified efficiency. Has
    no design or costing algorithm. Governing equations are
    
    .. math::
        Q_i X_i = Q_e X_e + Q_s X_s
        
        Q_i = Q_e + Q_s
    
        X_e = X_i * (1-e_rmv)
    
    where subscripts 'i', 'e', 's' represent influent, overflow effluent, and
    underflow sludge, respectively. 'Q' indicates volumetric flowrate and 'X' 
    indicates suspended solids concentration.

    Parameters
    ----------
    sludge_flow_rate : float, optional
        Underflow sludge flowrate [m3/d]. The default is 2000.
    solids_removal_efficiency : float, optional
        Removal efficiency (concentration basis) of suspended solids, unitless. 
        The default is 0.995.
    sludge_MLSS : float, optional
        Underflow MLSS [mg/L]. Used only when either `solids_removal_efficiency`
        or `sludge_flow_rate` is unspecified. The default is None.

    """

    _N_ins = 1
    _N_outs = 2  # [0] effluent overflow, [1] sludge underflow
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=2000, solids_removal_efficiency=0.995,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream',
                 F_BM_default=None, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.sludge_flow_rate = sludge_flow_rate
        self.solids_removal_efficiency = solids_removal_efficiency
        self.sludge_MLSS = sludge_MLSS
        self._mixed = WasteStream()
        self._f_uf = None
        self._f_of = None

    @property
    def sludge_flow_rate(self):
        '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
        return self._Qs

    @sludge_flow_rate.setter
    def sludge_flow_rate(self, Qs):
        self._Qs = Qs

    @property
    def solids_removal_efficiency(self):
        return self._e_rmv

    @solids_removal_efficiency.setter
    def solids_removal_efficiency(self, f):
        if f is not None and (f > 1 or f < 0):
            raise ValueError(f'solids removal efficiency must be within [0, 1], not {f}')
        self._e_rmv = f

    @property
    def sludge_MLSS(self):
        return self._MLSS

    @sludge_MLSS.setter
    def sludge_MLSS(self, MLSS):
        if MLSS is not None:
            warn(f'sludge MLSS {MLSS} mg/L is only used to estimate '
                 f'sludge flowrate or solids removal efficiency, when either '
                 f'one of them is unspecified.')
        self._MLSS = MLSS

    def _run(self):
        inf = self._mixed
        inf.mix_from(self.ins)
        of, uf = self.outs
        TSS_in = inf.get_TSS()
        if TSS_in <= 0:
            uf.empty()
            of.copy_like(inf)
        else:
            Q_in = inf.F_vol * 24 # m3/d
            x = inf.components.x
            Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
            if Qs and e_rmv:
                f_Qu = Qs/Q_in
                f_Xu = e_rmv + (1-e_rmv) * f_Qu
            elif Qs and mlss:
                f_Qu = Qs/Q_in
                f_Xu = f_Qu*mlss/TSS_in
            elif e_rmv and mlss:
                f_Qu = e_rmv / (mlss/TSS_in - (1-e_rmv))
                f_Xu = e_rmv + (1-e_rmv) * f_Qu
            split_to_uf = (1-x)*f_Qu + x*f_Xu
            if any(split_to_uf > 1): split_to_uf = 1
            inf.split_to(uf, of, split_to_uf)

    def _init_state(self):
        inf = self._mixed
        C_in = inf.conc
        Q_in = inf.F_vol * 24
        self._state = np.append(C_in, Q_in)
        self._dstate = self._state * 0.
        
    def _update_state(self):
        arr = self._state
        Cs = arr[:-1]
        Qi = arr[-1]
        Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
        x = self.components.x
        i_tss = x * self.components.i_mass

        of, uf = self.outs
        if uf.state is None: uf.state = np.zeros(len(x)+1)
        if of.state is None: of.state = np.zeros(len(x)+1)

        if Qs:
            Qe = Qi - Qs
            if e_rmv:
                fuf = e_rmv * Qi/Qs + (1-e_rmv)
                fof = 1-e_rmv
            elif mlss:
                tss_in = sum(Cs * i_tss)
                tss_e = (Qi * tss_in - Qs * mlss)/Qe
                fuf = mlss/tss_in
                fof = tss_e/tss_in
        elif e_rmv and mlss:
            tss_in = sum(Cs * i_tss)
            Qs = Qi * e_rmv / (mlss/tss_in - (1-e_rmv))
            Qe = Qi - Qs
            fuf = mlss/tss_in
            fof = 1-e_rmv
        else:
            raise RuntimeError('missing parameter')
            
        if Qs >= Qi: 
            uf.state[:] = arr
            of.state[:] = 0.
        else:
            self._f_uf = fuf
            self._f_of = fof
            uf.state[:-1] = Cs * ((1-x) + x*fuf)
            uf.state[-1] = Qs
            of.state[:-1] = Cs * ((1-x) + x*fof)
            of.state[-1] = Qe

    def _update_dstate(self):
        of, uf = self.outs
        x = self.components.x
        if uf.dstate is None: uf.dstate = np.zeros(len(x)+1)
        if of.dstate is None: of.dstate = np.zeros(len(x)+1)
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):        
        _state = self._state
        # _dstate = self._dstate
        _update_state = self._update_state
        # _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            # dQ_ins = dQC_ins[:, -1]
            # dC_ins = dQC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            # Q_dot = dQ_ins.sum()
            # C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            # _dstate[-1] = Q_dot
            # _dstate[:-1] = C_dot
            _update_state()
            # _update_dstate()
        self._AE = yt
   
    
# %%
# Total COD removal efficiency
nCOD = lambda f_corr, fx, HRT: f_corr*(2.88*fx - 0.118)*(1.45 + 6.15*np.log(HRT*24*60))

def calc_f_i(fx, f_corr, HRT):
    '''calculates the effluent-to-influent ratio of solid concentrations'''
    nX = nCOD(f_corr, fx, HRT)/fx
    if nX > 100: nX = 100
    if nX < 0: nX = 0
    return 1-(nX/100)
        
class PrimaryClarifierBSM2(SanUnit):
   
    """
    A Primary clarifier based on the Otterpohl model [1] in BSM2 [2]. 

    Parameters
    ----------
    ID : str
        ID for the clarifier.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 3.
    outs : class:`WasteStream`
        Sludge (uf) and treated effluent (of).
    volume : float, optional
        Clarifier volume, in m^3. The default is 900.
    ratio_uf : float
        The volumetric ratio of sludge to primary influent. The default is 0.007, 
        based on IWA report.[2]
    mean_f_x : float, optional
        The average fraction of particulate COD out of total COD in primary influent. 
        The default is 0.85.
    f_corr : float
        Dimensionless correction factor for removal efficiency in the primary clarifier.[2]    

    Examples
    --------
    >>> from qsdsan import set_thermo, Components, WasteStream
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.sanunits import PrimaryClarifierBSM2
    >>> PC = PrimaryClarifierBSM2(ID='PC', ins= (ws,), outs=('eff', 'sludge'),
    ...                           isdynamic=False)
    >>> PC.simulate()
    >>> of, uf = PC.outs
    >>> uf.imass['X_OHO']/ws.imass['X_OHO'] # doctest: +ELLIPSIS
    0.598...
    >>> PC.show()
    PrimaryClarifierBSM2: PC
    ins...
    [0] ws
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1e+04
                    S_NH4  2e+04
                    X_OHO  1.5e+04
                    H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         Alkalinity : 2.5 mmol/L
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
         TSS        : 11124.4 mg/L
    outs...
    [0] eff
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    9.93e+03
                    S_NH4  1.99e+04
                    X_OHO  6.03e+03
                    H2O    9.93e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 15436.7 mg/L
         BOD        : 10190.8 mg/L
         TC         : 5208.2 mg/L
         TOC        : 5208.2 mg/L
         TN         : 19890.1 mg/L
         TP         : 206.9 mg/L
         TK         : 27.8 mg/L
         TSS        : 4531.6 mg/L
    [1] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    70
                    S_NH4  140
                    X_OHO  8.97e+03
                    H2O    7e+03
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 693717.5 mg/L
         BOD        : 393895.7 mg/L
         TC         : 253653.4 mg/L
         TOC        : 253653.4 mg/L
         TN         : 57923.7 mg/L
         TP         : 13132.3 mg/L
         TK         : 3282.0 mg/L
         TSS        : 534594.0 mg/L
   
    References
    ----------
    [1] Otterpohl R. and Freund M. (1992). Dynamic Models for clarifiers of activated sludge 
    plants with dry and wet weather flows. Water Sci. Technol., 26(5-6), 1391-1400. 
    
    [2] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    """
    
    _N_ins = 3
    _N_outs = 2   # [0] effluent; [1] underflow
    _ins_size_is_fixed = False

    t_m = 0.125 # Smoothing time constant for qm calculation
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  isdynamic=True, init_with='WasteStream', 
                  volume=900, ratio_uf=0.007, mean_f_x=0.85, f_corr=0.65, 
                  F_BM=default_F_BM, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                         init_with=init_with)
        self._mixed = self.ins[0].copy(f'{ID}_mixed')
        self._sludge = np.ones(len(self.components)+1)
        self._effluent = np.ones(len(self.components)+1)
        self.V = volume
        self.ratio_uf = ratio_uf
        self.f_x = mean_f_x
        self.f_corr = f_corr
        self.F_BM.update(default_F_BM)
        self._concs = None

    @property
    def ratio_uf(self):
        return self._r
    @ratio_uf.setter
    def ratio_uf(self, r):
        if r > 1 or r < 0:
            raise ValueError(f'Sludge to influent ratio must be within [0, 1], not {r}')
        self._r = r
        self._sludge[-1] = r
        self._effluent[-1] = 1-r

    @property
    def f_x(self):
        '''[float] Fraction of particulate COD [-].'''
        if self._f_x: return self._f_x
        else:
            concs = self._mixed.conc
            cmps = self._mixed.components
            cod_concs = concs*cmps.i_COD
            if sum(cod_concs) == 0: return
            return sum(cod_concs*cmps.x)/sum(cod_concs)
    @f_x.setter
    def f_x(self, f):
        if isinstance(f, (float, int)) and (f < 0 or f > 1): 
            raise ValueError('f_x must be within [0,1]')
        self._f_x = f
                
    def _run(self):
        of, uf = self.outs
        mixed = self._mixed
        mixed.mix_from(self.ins)
        x = self.components.x
        r = self._r
        f_i = calc_f_i(self.f_x, self.f_corr, self.t_m)
        split_to_uf = (1-x)*r + x*(1-(1-r)*f_i)
        mixed.split_to(uf, of, split_to_uf)

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L].'''
        self._concs = self.components.kwarray(kwargs)

    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._concs is not None: Cs = self._concs
        else: Cs = mixed.conc
        self._state = np.append(Cs, Q).astype('float64')
        self._dstate = self._state * 0.
    
    def _update_parameters(self):
        x = self.components.x
        r = self._r
        Q = self._state[-1]
        f_i = calc_f_i(self.f_x, self.f_corr, self.V/Q)
        self._sludge[:-1] = x * ((1-f_i)/r+f_i) + (1-x)
        self._effluent[:-1] = x * f_i + (1-x)
    
    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Primary Clarifier'''
        of, uf = self.outs
        uf.state = self._sludge * self._state
        of.state = self._effluent * self._state
    
    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Primary Clarifier'''
        of, uf = self.outs
        uf.dstate = self._sludge * self._dstate
        of.dstate = self._effluent * self._dstate
       
    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE
    
    def _compile_ODE(self):
        _dstate = self._dstate
        _update_parameters = self._update_parameters
        _update_dstate = self._update_dstate
        V = self.V
        t_m = self.t_m
        def dy_dt(t, QC_ins, QC, dQC_ins):
            dydt_cstr(QC_ins, QC, V, _dstate)
            _dstate[-1] = (sum(QC_ins[:,-1])-QC[-1])/t_m
            _update_parameters()
            _update_dstate()
        self._ODE = dy_dt
        
#%%
# Assign a bare module of 1 to all
default_F_BM = {
        'Wall concrete': 1.,
        'Slab concrete': 1.,
        'Wall stainless steel': 1.,
        'Scraper': 1,
        'v notch weir': 1,
        'Pumps': 1
        }
default_F_BM.update(default_WWTpump_F_BM)

class PrimaryClarifier(IdealClarifier):
    
    """
    Primary clarifier with an ideal settling process model.
    
    Parameters
    ----------
    surface_overflow_rate : float
        Surface overflow rate in the clarifier in [(m3/day)/m2]. [1]
        Design SOR value for clarifier is 41 (m3/day)/m2 if it does not receive WAS.
        Design SOR value for clarifier is 29 (m3/day)/m2 if it receives WAS.
        Typically SOR lies between 30-50 (m3/day)/m2. 
        Here default value of 41 (m3/day)/m2 is used.
    depth_clarifier : float
        Depth of clarifier. Typical depths range from 3 m to 4.9 m [1], [2]. 
        Default value of 4.5 m would be used here. 
    downward_flow_velocity : float, optional
        Speed on the basis of which center feed diameter is designed [m/hr]. [3]
        The default is 36 m/hr. (10 mm/sec)
    F_BM : dict
        Equipment bare modules.
        
    Examples
    --------
    >>> from qsdsan import set_thermo, Components, WasteStream, System
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.sanunits import PrimaryClarifier
    >>> PC = PrimaryClarifier(ID='PC', ins=ws, outs=('effluent', 'sludge'),
    ...                       solids_removal_efficiency=0.6, 
    ...                       sludge_flow_rate=ws.F_vol*24*0.3,
    ...                       isdynamic=True)
    >>> sys = System('sys', path=(PC,))
    >>> sys.simulate(t_span=(0,10), method='BDF')  # doctest: +ELLIPSIS
    >>> PC.show() # doctest: +ELLIPSIS
    PrimaryClarifier: PC
    ins...
    [0] ws
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1e+04
                    S_NH4  2e+04
                    X_OHO  1.5e+04
                    H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         Alkalinity : 2.5 mmol/L
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
         TSS        : 11124.4 mg/L
    outs...
    [0] effluent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    7e+03
                    S_NH4  1.4e+04
                    X_OHO  4.2e+03
                    H2O    7.04e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 15278.7 mg/L
         BOD        : 10093.4 mg/L
         TC         : 5152.8 mg/L
         TOC        : 5152.8 mg/L
         TN         : 19776.2 mg/L
         TP         : 204.4 mg/L
         TK         : 27.3 mg/L
         TSS        : 4449.8 mg/L
    [1] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    3e+03
                    S_NH4  6e+03
                    X_OHO  1.08e+04
                    H2O    2.96e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 43926.3 mg/L
         BOD        : 26326.2 mg/L
         TC         : 15637.8 mg/L
         TOC        : 15637.8 mg/L
         TN         : 21732.8 mg/L
         TP         : 748.7 mg/L
         TK         : 163.9 mg/L
         TSS        : 26698.6 mg/L
    
    References
    ----------
    [1] Chapter-10: Primary Treatment. Design of water resource recovery facilities. 
    WEF Manual of Practice No. 8. 6th Edition. Virginia: McGraw-Hill, 2018. 
    
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    
    [3] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    """
    
    _N_ins = 1
    _N_outs = 2  # [0] overflow effluent [1] underflow sludge
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True
    
    peak_flow_safety_factor = 2.5 # assumed based on average and maximum velocities

    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    def __init__(self, ID='', ins=None, outs=(), 
                 sludge_flow_rate=2000, solids_removal_efficiency=0.6,
                 sludge_MLSS=None, thermo=None, isdynamic=False, 
                 init_with='WasteStream', 
                 surface_overflow_rate=41, depth_clarifier=4.5,
                 downward_flow_velocity=36, F_BM=default_F_BM, **kwargs):
        super().__init__(ID, ins, outs, thermo,
                         sludge_flow_rate=sludge_flow_rate, 
                         solids_removal_efficiency=solids_removal_efficiency,
                         sludge_MLSS=sludge_MLSS,
                         isdynamic=isdynamic, 
                         init_with=init_with)

        self.surface_overflow_rate = surface_overflow_rate
        self.depth_clarifier = depth_clarifier
        self.downward_flow_velocity = downward_flow_velocity
        self.F_BM.update(F_BM)
        self._sludge = uf = WasteStream(f'{ID}_sludge')
        pump_id = f'{ID}_sludge_pump'
        self.sludge_pump = WWTpump(
            ID=pump_id, ins=uf, thermo=thermo, pump_type='sludge', 
            prefix='Sludge',
            include_pump_cost=True,
            include_building_cost=False,
            include_OM_cost=True, 
            )
         
        # self.auxiliary_unit_names = tuple({*self.auxiliary_unit_names, pump_id})
            
    # @property
    # def solids_loading_rate(self):
    #     '''solids_loading_rate is the loading in the clarifier'''
    #     return self._slr

    # @solids_loading_rate.setter
    # def solids_loading_rate(self, slr):
    #     if slr is not None:
    #         self._slr = slr
    #     else: 
    #         raise ValueError('solids_loading_rate of the clarifier expected from user')
            
    def _design_pump(self):
        D = self.design_results
        N = D['Number of pumps']
        pump = self.sludge_pump
        self._sludge.copy_like(self.outs[1])
        self._sludge.scale(1/N)
        pump.simulate()
        D.update(pump.design_results)
    
    _units = {
        'Number of clarifiers': 'ea',
        'SOR': 'm3/day/m2',
        'Volumetric flow': 'm3/day',
        'Surface area': 'm2',
        'Cylindrical diameter': 'm',
        'Conical radius': 'm',
        'Conical depth': 'm',
        'Clarifier depth': 'm',
        'Cylindrical depth': 'm',
        'Cylindrical volume': 'm3',
        'Conical volume': 'm3',
        'Volume': 'm3',
        'Hydraulic Retention Time': 'hr', 
        'Center feed depth': 'm',
        'Downward flow velocity': 'm/hr',
        'Center feed diameter': 'm',
        'Volume of concrete wall': 'm3',
        'Volume of concrete slab': 'm3',
        'Stainless steel': 'kg',
        # 'Pump pipe stainless steel' : 'kg',
        # 'Pump stainless steel': 'kg',
        'Number of pumps': 'ea'
    }
    
    density_ss = 7930 # kg/m3, 18/8 Chromium

    def _design(self):    
        mixed = self._mixed
        mixed.mix_from(self.ins)
        D = self.design_results
        
        # Number of clarifiers based on tentative suggestions by Jeremy 
        # (would be verified through collaboration with industry)
        Q_mgd = mixed.get_total_flow('MGD')
        if Q_mgd <= 3: N = 2
        elif Q_mgd <= 8: N = 3
        elif Q_mgd <= 20: N = 4
        else: N = 3 + int(Q_mgd / 20)
        D['Number of clarifiers'] = D['Number of pumps'] = N
        
        SOR = D['SOR'] = self.surface_overflow_rate # in (m3/day)/m2
        Q = D['Volumetric flow'] = mixed.get_total_flow('m3/d')/N # m3/day
        A = D['Surface area'] = Q/SOR # in m2
        dia = D['Cylindrical diameter'] = sqrt(4*A/pi) #in m
        
        # Check on cylindrical diameter d [2, 3]
        if dia < 3 or dia > 60:
            warn(f'Cylindrical diameter = {dia:.2f} is not between 3 m and 60 m')
        
        rad = D['Conical radius'] = dia/2
        # The slope of the bottom conical floor lies between 1:10 to 1:12 [3, 4]
        h_cone = D['Conical height'] = rad/12 
        h = D['Clarifier depth'] = self.depth_clarifier # in m 
        h_cyl = D['Cylindrical height'] = h - h_cone
        
        # Check on cylindrical and conical depths 
        if h_cyl < h_cone:
            warn(f'Cylindrical highet = {h_cyl} is lower than conical height = {h_cone}')
        
        V_cyl = D['Cylindrical volume'] = A * h_cyl     # in m3
        V_cone = D['Conical volume'] = A * h_cone / 3   # in m3
        V = D['Volume'] = V_cyl + V_cone                # in m3
        
        HRT = D['Hydraulic Retention Time'] = V/(Q/24)        # in hrs
        
        # Check on cylinderical HRT [3]
        if HRT < 1.5 or HRT > 2.5:
            warn(f'HRT = {HRT} is not between 1.5 and 2.5 hrs')
        
        # The design here is for center feed of clarifier.
        
        # Depth of the center feed lies between 30-75% of sidewater depth. [3, 4]
        h_cf = D['Center feed depth'] = 0.5*h_cyl
        # Typical conventional feed wells are designed for an average downflow velocity
        # of 10-13 mm/s and maximum velocity of 25-30 mm/s. [4]
        v_down = D['Downward flow velocity'] = self.downward_flow_velocity*self.peak_flow_safety_factor # in m/hr
        
        A_cf = (Q/24)/v_down # in m2
        
        dia_cf = D['Center feed diameter'] = sqrt(4*A_cf/pi)

        #Sanity check: Diameter of the center feed lies between 15-25% of tank diameter [4]
        #The lower limit of this check has been modified to 10% based on advised range of down flow velocity in [4]. 
        if dia_cf < 0.10*dia or dia_cf  > 0.25*dia:
            warn(f'Diameter of the center feed does not lie between 15-25% of tank diameter. It is {dia_cf*100/dia:.2f}% of tank diameter')

        # Amount of concrete required
        # D_tank = D['Cylindrical depth']*39.37 # m to inches
        h_ft = h*3.2808398950131235         # m to feet
        # Thickness of the wall concrete [m]. Default to be minimum of 1 feet with 1 inch added for every feet of depth over 12 feet.
        # thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
        d_wall = (1 + max(h_ft-12, 0)/12) * 0.3048  # feet to m
        OD = dia + 2*d_wall
        D['Volume of concrete wall'] = pi*h_cyl/4*(OD**2 - dia**2)  # m3
        
        # Concrete slab thickness, [ft], default to be 2 in thicker than the wall thickness. (Brian's code)
        d_slab = d_wall + (2/12)*0.3048 # from inch to m
        # outer_diameter_cone = inner_diameter + 2*(thickness_concrete_wall + thickness_concrete_slab)
        OD_cone = dia + 2*d_slab
        # volume_conical_wall = (np.pi/(3*4))*(((D['Conical depth'] + thickness_concrete_wall + thickness_concrete_slab)*(outer_diameter_cone**2)) - (D['Conical depth']*(inner_diameter)**2))
        # D['Volume of concrete slab'] = volume_conical_wall
        D['Volume of concrete slab'] = pi/3*((h_cone + d_slab)*(OD_cone/2)**2 - h_cone*(dia/2)**2)
        
        # Amount of metal required for center feed
        #!!! consider empirical estimation of steel volume for all equipment (besides center feed, e.g., scrapper, support column, EDI, skimmer, walkway etc.)
        d_wall_cf = 0.3048 # equal to 1 feet, in m (!! NEED A RELIABLE SOURCE !!)
        OD_cf = dia_cf + 2*d_wall_cf
        volume_center_feed = (pi*h_cf/4)*(OD_cf**2 - dia_cf**2)
        D['Stainless steel'] = volume_center_feed*self.density_ss # in kg
       
        # Pumps
        self._design_pump()
        
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        N = D['Number of clarifiers']
        
        # Construction of concrete and stainless steel walls
        C['Wall concrete'] = N*D['Volume of concrete wall']*self.wall_concrete_unit_cost
        C['Slab concrete'] = N*D['Volume of concrete slab']*self.slab_concrete_unit_cost
        C['Wall stainless steel'] = N*D['Stainless steel']*self.stainless_steel_unit_cost
       
        # Cost of equipment 
        
        # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
        # Scraper 
        # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
        # base_cost_scraper = 2500
        # base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
        Q = D['Volumetric flow']/24
        
        # C['Scraper'] =  D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # base_power_scraper = 2.75 # in kW
        # THE EQUATION BELOW IS NOT CORRECT TO SCALE SCRAPER POWER REQUIREMENTS 
        # scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # v notch weir
        # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
        base_cost_v_notch_weir = 6888
        base_flow_v_notch_weir = 10 # in m3/hr
        C['v notch weir'] = N*base_cost_v_notch_weir*(Q/base_flow_v_notch_weir)**0.6
       
        # Pump (construction and maintainance)
        pump = self.sludge_pump
        add_OPEX = self.add_OPEX
        add_OPEX.update({k: v*N for k,v in pump.add_OPEX.items()})
        C.update({k: v*N for k,v in pump.baseline_purchase_costs.items()})
        self.power_utility.rate += pump.power_utility.rate*N
        # self.power_utility.rate += scraper_power
