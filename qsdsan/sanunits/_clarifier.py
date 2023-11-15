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
from warnings import warn
from .. import SanUnit, WasteStream
import numpy as np
from ..sanunits import WWTpump
from ..sanunits._pumping import default_F_BM as default_WWTpump_F_BM

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

def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp, n0):
    X_star = npmax(X-X_min, n0)
    v = npmin(v_max_practical, v_max*(npexp(-rh*X_star) - npexp(-rp*X_star)))
    return X*npmax(v, n0)

class FlatBottomCircularClarifier(SanUnit):
    """
    A flat-bottom circular clarifier with a simple 1-dimensional
    N-layer settling model. [1]_

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
    .. [1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the Clarification
        -Thickening Process. Water Res. 1991, 25 (10), 1263–1271.
        https://doi.org/10.1016/0043-1354(91)90066-Y.
    .. [2] Chapter-12: Suspended-growth Treatment Processes. WEF Manual of Practice No. 8. 
        6th Edition. Virginia: McGraw-Hill, 2018. 
    .. [3] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    .. [4] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    .. [5] RECOMMENDED STANDARDS for WASTEWATER FACILITIES. 10 state standards. 2014 edition. 
    """

    _N_ins = 1
    _N_outs = 3
    
    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    pumps = ('ras', 'was',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', underflow=2000, wastage=385,
                 surface_area=1500, height=4, N_layer=10, feed_layer=4,
                 X_threshold=3000, v_max=474, v_max_practical=250,
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3, F_BM_default=default_F_BM, isdynamic=True,
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
        Q_e = arr[-(1+n)] - self._Qras - self._Qwas
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
        vmax_arr = np.full_like(nzeros, self._v_max)
        vmaxp_arr = np.full_like(nzeros, self._v_max_p)
        rh_arr = np.full_like(nzeros, self._rh)
        rp_arr = np.full_like(nzeros, self._rp)
        func_vx = lambda x_arr, xmin_arr : _settling_flux(x_arr, vmax_arr, vmaxp_arr, xmin_arr, rh_arr, rp_arr, nzeros)
       
        A, hj, V = self._A, self._hj, self._V
        A_arr = np.full_like(nzeros, A)
        hj_arr = np.full_like(nzeros, hj)
        J = np.zeros(n-1)
        X_t_arr = np.full(jf, self._X_t)
        Q_in_arr = np.zeros(m)
        V_arr = np.full(m, V)

        def dy_dt(t, QC_ins, QC, dQC_ins):
            dQC[-(n+1)] = dQC_ins[0,-1]
            Q_in = QC_ins[0,-1]
            Q_e = Q_in - Q_s
            C_in = QC_ins[0,:-1]
            dC_in = dQC_ins[0,:-1]
            Z_in = C_in*(1-x)
            X_in = sum(C_in*imass*x)           # influent TSS
            dX_in = sum(dC_in*imass*x)
            X_min_arr[:] = X_in * fns
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
            condition = (X_rolled[:jf]<X_t_arr)
            J[:jf][condition] = VX[:jf][condition]
            settle_out[:-1] = J
            settle_in[1:] = J
            dQC[-n:] = ((flow_in - flow_out)/A_arr + settle_in - settle_out)/hj_arr       # (n,)
            #*********solubles**********
            Q_in_arr[:] = Q_in
            dQC[:m] = Q_in_arr*(Z_in - Z)/V_arr
            # instrumental variables
            dX_comp[:] = (dC_in * X_in - dX_in * C_in) * x / X_in**2
            _update_dstate()

        self._ODE = dy_dt
    
    _units = {
        'Number of clarifiers': 'ea',
        'Volumetric flow': 'm3/day',
        'Clarifier depth': 'm',
        'Surface area': 'm2',
        'Clarifier diameter': 'm',
        'Clarifier volume': 'm3',
        'Design solids loading rate': 'kg/m2/hr',
        'Surface overflow rate': 'm3/day/m2',
        'Hydraulic Retention Time': 'hr', 
        'Center feed depth': 'm',
        'Downward flow velocity': 'm/hr',
        'Center feed diameter': 'm',
        'Volume of concrete wall': 'm3',
        'Stainless steel': 'kg',
        'Pump pipe stainless steel' : 'kg',
        'Pump stainless steel': 'kg',
        'Number of pumps': 'ea'
    }
     
    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
    
        self._ras.copy_like(self.outs[1])
        self._was.copy_like(self.outs[2])
        
        ins_dct = {
            'ras': self._ras,
            'was': self._was,
            }
        
        D = self.design_results
        
        ras_flow = self._ras.get_total_flow('m3/hr')
        was_flow = self._was.get_total_flow('m3/hr')
        
        ras_flow_u = ras_flow/D['Number of clarifiers']*0.00634
        was_flow_u = was_flow/D['Number of clarifiers']*0.00634
        
        Q_mgd = {
            'ras': ras_flow_u,
            'was': was_flow_u,
            }
        
        type_dct = dict.fromkeys(pumps, 'sludge')
        inputs_dct = dict.fromkeys(pumps, (1,))
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                pump = WWTpump(
                    ID=ID, ins=ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=Q_mgd[i], add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
        return pipe_ss, pump_ss
     
    def _design(self):
        
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
        D = self.design_results
        
        # Number of clarifiers based on tentative suggestions by Jeremy 
        # (would be verified through collaboration with industry)
        total_flow = (mixed.get_total_flow('m3/hr')*24)/3785 # in MGD
        if total_flow <= 3:
            D['Number of clarifiers'] = 2
        elif total_flow > 3 and total_flow <= 8:
            D['Number of clarifiers'] = 3
        elif total_flow > 8 and total_flow <=20:
            D['Number of clarifiers'] = 4
        else:
            D['Number of clarifiers'] = 4
            total_flow -= 20
            D['Number of clarifiers'] += np.ceil(total_flow/20)
                
        D['Volumetric flow'] =  (mixed.get_total_flow('m3/hr')*24)/D['Number of clarifiers'] #m3/day
        
        # Sidewater depth of a cylindrical clarifier lies between 4-5 m (MOP 8)
        D['Clarifier depth'] = self._h # in m
        
        # Area of clarifier 
        # D['Surface area'] = solids_clarifier/D['Solids loading rate'] #m2
        D['Surface area'] = self._A/D['Number of clarifiers']
        D['Clarifier diameter'] = np.sqrt(4*D['Surface area']/np.pi) # in m
        D['Clarifier volume'] = D['Surface area']*D['Clarifier depth'] # in m3
        
        # Checks on SLR,, SOR, and HRT 
        
        D['Design solids loading rate'] = self._slr # kg/(m2*hr)
        
        total_solids = mixed.get_TSS()*mixed.get_total_flow('m3/hr')/1000 # in kg/hr (mg/l * m3/hr)
        solids_clarifier = total_solids/D['Number of clarifiers'] # in kg/hr
        simulated_slr = solids_clarifier/D['Surface area'] # in kg/(m2*hr)
        
        # Consult Joy on the margin or error
        if simulated_slr < 0.8*D['Design solids loading rate'] or simulated_slr > 1.2*D['Design solids loading rate']:
            design_slr = D['Design solids loading rate']
            warn(f'Solids loading rate = {simulated_slr} is not within 20% of the recommended design level of {design_slr} kg/hr/m2')
        
        # Check on SLR [3, 4, 5] 
        if simulated_slr > 14:
            warn(f'Solids loading rate = {simulated_slr} is above recommended level of 14 kg/hr/m2')
        
        # Check on SOR [3, 4, 5]
        D['Surface overflow rate'] = D['Volumetric flow']/D['Surface area']  # in m3/m2/hr
        if D['Surface overflow rate'] > 49:
            sor = D['Surface overflow rate']
            warn(f'Surface overflow rate = {sor} is above recommended level of 49 m3/day/m2')
        
        # HRT
        D['Hydraulic Retention Time'] = D['Clarifier volume']*24/D['Volumetric flow'] # in hr
        
        # Clarifiers can be center feed or peripheral feed. The design here is for the more commonly deployed center feed.
        # Depth of the center feed lies between 30-75% of sidewater depth. [2]
        D['Center feed depth'] = 0.5*D['Clarifier depth']
        # Criteria for downward velocity of flow determine 
        D['Downward flow velocity'] = self._downward_flow_velocity # in m/hr
        Center_feed_area = (D['Volumetric flow']/24)/D['Downward flow velocity'] # in m2
        D['Center feed diameter'] = np.sqrt(4*Center_feed_area/np.pi) 

        #Sanity check: Diameter of the center feed lies between 20-25% of tank diameter [2]
        if D['Center feed diameter'] < 0.20*D['Clarifier diameter'] or D['Center feed diameter']  > 0.25*D['Clarifier diameter']:
            cf_dia = D['Center feed diameter'] 
            tank_dia = D['Clarifier diameter']
            warn(f'Diameter of the center feed does not lie between 20-25% of tank diameter. It is {cf_dia*100/tank_dia} % of tank diameter')
            
        # Amount of concrete required
        D_tank = D['Clarifier depth']*39.37 # m to inches 
        # Thickness of the wall concrete, [m]. Default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft. (Brian's code)
        thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
        inner_diameter = D['Clarifier diameter']
        outer_diameter = inner_diameter + 2*thickness_concrete_wall
        D['Volume of concrete wall']  = (np.pi*D['Clarifier depth']/4)*(outer_diameter**2 - inner_diameter**2)
        
        # Concrete slab thickness, [ft], default to be 2 in thicker than the wall thickness. (Brian's code)
        thickness_concrete_slab = thickness_concrete_wall + (2/12)*0.3048 # from inch to m
        # From Brian's code
        D['Volume of concrete slab']  = (thickness_concrete_slab + thickness_concrete_wall)*D['Surface area']
        
        # Amount of metal required for center feed
        thickness_metal_wall = 0.3048 # equal to 1 feet, in m (!! NEED A RELIABLE SOURCE !!)
        inner_diameter_center_feed = D['Center feed diameter']
        outer_diameter_center_feed = inner_diameter_center_feed + 2*thickness_metal_wall
        volume_center_feed = (np.pi*D['Center feed depth']/4)*(outer_diameter_center_feed**2 - inner_diameter_center_feed **2)
        density_ss = 7930 # kg/m3, 18/8 Chromium
        D['Stainless steel'] = volume_center_feed*density_ss # in kg
       
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        # For secondary clarifier
        D['Number of pumps'] = 2*D['Number of clarifiers']
        
    def _cost(self):
       
        D = self.design_results
        C = self.baseline_purchase_costs
       
        # Construction of concrete and stainless steel walls
        C['Wall concrete'] = D['Number of clarifiers']*D['Volume of concrete wall']*self.wall_concrete_unit_cost
        
        C['Slab concrete'] = D['Number of clarifiers']*D['Volume of concrete slab']*self.slab_concrete_unit_cost
        
        C['Wall stainless steel'] = D['Number of clarifiers']*D['Stainless steel']*self.stainless_steel_unit_cost
        
        # Cost of equipment 
        
        # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
        # Scraper 
        # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
        base_cost_scraper = 2500
        base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
        clarifier_flow =  D['Volumetric flow']/24 # in m3/hr
        C['Scraper'] = D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
        base_power_scraper = 2.75 # in kW
        scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # v notch weir
        # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
        base_cost_v_notch_weir = 6888
        base_flow_v_notch_weir = 10 # in m3/hr
        C['v notch weir'] = D['Number of clarifiers']*base_cost_v_notch_weir*(clarifier_flow/base_flow_v_notch_weir)**0.6
       
        # Pump (construction and maintainance)
        pumps = self.pumps
        add_OPEX = self.add_OPEX
        pump_cost = 0.
        building_cost = 0.
        opex_o = 0.
        opex_m = 0.
       
        # i would be 0 and 1 for RAS and WAS respectively 
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost = p.baseline_purchase_costs
            p_add_opex = p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']
            
        # All costs associated with pumping need to be multiplied by number of clarifiers 
        C['Pumps'] = pump_cost*D['Number of clarifiers']
        C['Pump building'] = building_cost*D['Number of clarifiers']
        add_OPEX['Pump operating'] = opex_o*D['Number of clarifiers']
        add_OPEX['Pump maintenance'] = opex_m*D['Number of clarifiers']
       
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
        
        pumping = pumping*D['Number of clarifiers']
        
        self.power_utility.rate += pumping
        self.power_utility.consumption += scraper_power
        
# %% 
   
class IdealClarifier(SanUnit):

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=2000, solids_removal_efficiency=.995,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream',
                 F_BM_default=None, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.sludge_flow_rate = sludge_flow_rate
        self.solids_removal_efficiency = solids_removal_efficiency
        self.sludge_MLSS = sludge_MLSS

    @property
    def sludge_flow_rate(self):
        '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
        return self._Qs

    @sludge_flow_rate.setter
    def sludge_flow_rate(self, Qs):
        if Qs is not None: self._Qs = Qs
        elif self.ins[0].isempty(): self._Qs = None
        else: self._Qs = self._calc_Qs()

    @property
    def solids_removal_efficiency(self):
        return self._e_rmv

    @solids_removal_efficiency.setter
    def solids_removal_efficiency(self, f):
        if f is not None:
            if f > 1 or f < 0:
                raise ValueError(f'solids removal efficiency must be within [0, 1], not {f}')
            self._e_rmv = f
        elif self.ins[0].isempty(): self._e_rmv = None
        else: self._e_rmv = self._calc_ermv()

    @property
    def sludge_MLSS(self):
        return self._MLSS

    @sludge_MLSS.setter
    def sludge_MLSS(self, MLSS):
        if MLSS is not None: self._MLSS = MLSS
        elif self.ins[0].isempty(): self._MLSS = None
        else: self._MLSS = self._calc_SS()[1]

    def _calc_Qs(self, TSS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if TSS_in is None: TSS_in = self.ins[0].get_TSS()
        return Q_in*TSS_in*self._e_rmv/(self._MLSS-TSS_in)

    def _calc_ermv(self, TSS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if TSS_in is None: TSS_in = self.ins[0].get_TSS()
        return self._Qs*(self._MLSS-TSS_in)/TSS_in/(Q_in-self._Qs)

    def _calc_SS(self, SS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if SS_in is None: SS_in = self.ins[0].get_TSS()
        SS_e = (1-self._e_rmv)*SS_in
        Qs = self._Qs
        Qe = Q_in - Qs
        return SS_e, (Q_in*SS_in - Qe*SS_e)/Qs

    def _run(self):
        inf, = self.ins
        eff, sludge = self.outs
        cmps = self.components
        Q_in = inf.get_total_flow('m3/d')
        TSS_in = (inf.conc*cmps.x*cmps.i_mass).sum()
        params = (Qs, e_rmv, MLSS) = self._Qs, self._e_rmv, self._MLSS
        if sum([i is None for i in params]) > 1:
            raise RuntimeError('must specify two of the following parameters: '
                               'sludge_flow_rate, solids_removal_efficiency, sludge_MLSS')
        if Qs is None:
            Qs = self._calc_Qs(TSS_in, Q_in)
            Xs = MLSS / TSS_in * inf.conc * cmps.x
            Xe = (1-e_rmv) * inf.conc * cmps.x
        elif e_rmv is None:
            e_rmv = self._calc_ermv(TSS_in, Q_in)
            Xs = MLSS / TSS_in * inf.conc * cmps.x
            Xe = (1-e_rmv) * inf.conc * cmps.x
        else:
            Xe, Xs = self._calc_SS(inf.conc * cmps.x, Q_in)
        Zs = Ze = inf.conc * (1-cmps.x)
        Ce = dict(zip(cmps.IDs, Ze+Xe))
        Cs = dict(zip(cmps.IDs, Zs+Xs))
        Ce.pop('H2O', None)
        Cs.pop('H2O', None)
        eff.set_flow_by_concentration(Q_in-Qs, Ce, units=('m3/d', 'mg/L'))
        sludge.set_flow_by_concentration(Qs, Cs, units=('m3/d', 'mg/L'))

    def _design(self):
        pass
   
    
# %%    
class PrimaryClarifierBSM2(SanUnit):
   
    """
    A Primary clarifier based on BSM2 Layout. [1]

    Parameters
    ----------
    ID : str
        ID for the clarifier.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 3.
    outs : class:`WasteStream`
        Sludge (uf) and treated effluent (of).
    Hydraulic Retention time : float
        Hydraulic Retention Time in days. The default is 0.04268 days, based on IWA report.[1]
    ratio_uf : float
        The ratio of sludge to primary influent. The default is 0.007, based on IWA report.[1]
    f_corr : float
        Dimensionless correction factor for removal efficiency in the primary clarifier.[1]
    
    # cylindrical_depth : float, optional
    #     The depth of the cylindrical portion of clarifier [in m].  
    # upflow_velocity : float, optional
    #     Speed with which influent enters the center feed of the clarifier [m/hr]. The default is 43.2.
    # F_BM : dict
    #     Equipment bare modules.

    Examples
    --------
    >>> from qsdsan import set_thermo, Components, WasteStream
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.sanunits import PrimaryClarifierBSM2
    >>> PC = PrimaryClarifierBSM2(ID='PC', ins= (ws,), outs=('eff', 'sludge'))
    >>> PC.simulate()
    >>> uf, of = PC.outs
    >>> uf.imass['X_OHO']/ws.imass['X_OHO'] # doctest: +ELLIPSIS
    0.280...
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
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
    outs...
    [0] eff
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    70
                    S_NH4  140
                    X_OHO  4.2e+03
                    H2O    7e+03
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 428873.3 mg/L
         BOD        : 244072.7 mg/L
         TC         : 156644.5 mg/L
         TOC        : 156644.5 mg/L
         TN         : 43073.0 mg/L
         TP         : 8085.4 mg/L
         TK         : 2011.4 mg/L
    [1] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    9.93e+03
                    S_NH4  1.99e+04
                    X_OHO  1.08e+04
                    H2O    9.93e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 19982.3 mg/L
         BOD        : 12762.2 mg/L
         TC         : 6873.2 mg/L
         TOC        : 6873.2 mg/L
         TN         : 20145.0 mg/L
         TP         : 293.5 mg/L
         TK         : 49.6 mg/L
   
    References
    ----------
    [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    [3] Otterpohl R. and Freund M. (1992). Dynamic Models for clarifiers of activated sludge 
    plants with dry and wet weather flows. Water Sci. Technol., 26(5-6), 1391-1400. 
    """
    
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False

    # Costs
    wall_concrete_unit_cost = 650 / 0.765 # $/m3, 0.765 is to convert from $/yd3
    stainless_steel_unit_cost=1.8 # $/kg (Taken from Joy's METAB code) https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
   
    pumps = ('sludge',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  isdynamic=False, init_with='WasteStream', Hydraulic_Retention_Time=0.04268,
                  ratio_uf=0.007, f_corr=0.65, cylindrical_depth = 5, upflow_velocity = 43.2, 
                  F_BM=default_F_BM, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                          init_with=init_with)
        self.Hydraulic_Retention_Time = Hydraulic_Retention_Time #in days
        self.ratio_uf = ratio_uf
        self.f_corr = f_corr
        self.cylindrical_depth = cylindrical_depth # in m 
        self.upflow_velocity = upflow_velocity # in m/hr (converted from 12 mm/sec)
        self.F_BM.update(default_F_BM)
        self._mixed = self.ins[0].copy(f'{ID}_mixed')
        self._sludge = self.outs[1].copy(f'{ID}_sludge')
       
    @property
    def Hydraulic_Retention_Time(self):
        '''The Hydraulic Retention time in days.'''
        return self._HRT

    @Hydraulic_Retention_Time.setter
    def Hydraulic_Retention_Time(self, HRT):
        if HRT is not None:
            self._HRT = HRT
        else:
            raise ValueError('HRT expected from user')

    @property
    def ratio_uf(self):
        return self._r

    @ratio_uf.setter
    def ratio_uf(self, r):
        if r is not None:
            if r > 1 or r < 0:
                raise ValueError(f'Sludge to influent ratio must be within [0, 1], not {r}')
            self._r = r
        else:
            raise ValueError('Sludge to influent ratio expected from user')
           
    @property
    def f_corr(self):
        return self._corr

    @f_corr.setter
    def f_corr(self, corr):
        if corr is not None:
            # if corr > 1 or corr < 0:
            #     raise ValueError(f'correction factor must be within [0, 1], not {corr}')
            self._corr = corr
        else:
            raise ValueError('correction factor expected from user')
   
    def _f_i(self):
        xcod = self._mixed.composite('COD', particle_size='x')
        fx = xcod/self._mixed.COD
        corr = self._corr
        HRT = self._HRT
        n_COD = corr*(2.88*fx - 0.118)*(1.45 + 6.15*np.log(HRT*24*60))
        f_i = 1 - (n_COD/100)
        return f_i
   
    def _run(self):
        uf, of = self.outs
        cmps = self.components
        mixed = self._mixed
        mixed.mix_from(self.ins)
    
        r = self._r
        f_i = self._f_i()
       
        Xs = (1 - f_i)*mixed.mass*cmps.x
        Xe = (f_i)*mixed.mass*cmps.x
       
        Zs = r*mixed.mass*cmps.s
        Ze = (1-r)*mixed.mass*cmps.s
       
        Ce = Ze + Xe
        Cs = Zs + Xs
        of.set_flow(Ce,'kg/hr')
        uf.set_flow(Cs,'kg/hr')
       
    def _init_state(self):
        # if multiple wastestreams exist then concentration and total flow
        # would be calculated assuming perfect mixing
        Qs = self._ins_QC[:,-1]
        Cs = self._ins_QC[:,:-1]
        self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.
       
        uf, of = self.outs
        s_flow = uf.F_vol/(uf.F_vol+of.F_vol)
        denominator = uf.mass + of.mass
        denominator += (denominator == 0)
        s = uf.mass/denominator
        self._sludge = np.append(s/s_flow, s_flow)
        self._effluent = np.append((1-s)/(1-s_flow), 1-s_flow)
       
    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Primary Clarifier'''
        self._outs[0].state = self._sludge * self._state
        self._outs[1].state = self._effluent * self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Primary Clarifier'''
        self._outs[0].dstate = self._sludge * self._dstate
        self._outs[1].dstate = self._effluent * self._dstate
     
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            #Because there are multiple inlets
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            dQ_ins = dQC_ins[:, -1]
            dC_ins = dQC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            Q_dot = dQ_ins.sum()
            C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            _dstate[-1] = Q_dot
            _dstate[:-1] = C_dot
            _update_state()
            _update_dstate()
        self._AE = yt

   
    # _units = {
    #     'Number of clarifiers': '?',
    #     'Cylindrical volume': 'm3',
    #     'Cylindrical depth': 'm',
    #     'Cylindrical diameter': 'm',
       
    #     'Conical radius': 'm',
    #     'Conical depth': 'm',
    #     'Conical volume': 'm3',
       
    #     'Volume': 'm3',
    #     'Center feed depth': 'm',
    #     'Upflow velocity': 'm/hr',
    #     'Center feed diameter': 'm',
    #     'Volume of concrete wall': 'm3',
    #     'Stainless steel': 'kg',
    #     'Pump pipe stainless steel' : 'kg',
    #     'Pump stainless steel': 'kg',
    #     'Number of pumps': '?'
    # }
   
   
    # def _design_pump(self):
    #     ID, pumps = self.ID, self.pumps
    #     self._sludge.copy_like(self.outs[1])
        
    #     ins_dct = {
    #         'sludge': self._sludge,
    #         }
       
    #     type_dct = dict.fromkeys(pumps, 'sludge')
    #     inputs_dct = dict.fromkeys(pumps, (1,),)
        
    #     D = self.design_results
    #     influent_Q = self._sludge.get_total_flow('m3/hr')/D['Number of clarifiers']
    #     influent_Q_mgd = influent_Q*0.00634 # m3/hr to MGD 
       
    #     for i in pumps:
    #         if hasattr(self, f'{i}_pump'):
    #             p = getattr(self, f'{i}_pump')
    #             setattr(p, 'add_inputs', inputs_dct[i])
    #         else:
    #             ID = f'{ID}_{i}'
    #             capacity_factor=1
    #             pump = WWTpump(
    #                 ID=ID, ins=ins_dct[i], pump_type=type_dct[i],
    #                 Q_mgd=influent_Q_mgd, add_inputs=inputs_dct[i],
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
       
    #     D = self.design_results
    #     total_flow = self._mixed.get_total_flow('m3/hr')
        
    #     if total_flow <= 1580: # 10 MGD
    #         design_flow = 790  # 5 MGD
    #     elif total_flow >1580 and total_flow <= 4730: # Between 10 and 30 MGD
    #         design_flow = 2365 # 15 MGD
    #     elif total_flow > 4730 and total_flow <= 15770: # Between 30 and 100 MGD
    #         design_flow = 3940 # 25 MGD
    #     else:
    #         design_flow = 5520 # 35 MGD 
        
    #     D['Number of clarifiers'] = np.ceil(total_flow/design_flow)
       
    #     total_volume = 24*self._HRT*design_flow #in m3
    #     working_volume = total_volume/0.8 # Assume 80% working volume
       
    #     D['Cylindrical volume'] = working_volume
    #     # Sidewater depth of a cylindrical clarifier lies between 2.5-5m
    #     D['Cylindrical depth'] = self.cylindrical_depth # in m
    #     # The tank diameter can lie anywhere between 3 m to 100 m
    #     D['Cylindrical diameter'] = (4*working_volume/(3.14*D['Cylindrical depth']))**(1/2) # in m
       
    #     D['Conical radius'] = D['Cylindrical diameter']/2
    #     # The slope of the bottom conical floor lies between 1:10 to 1:12
    #     D['Conical depth'] = D['Conical radius']/10
    #     D['Conical volume'] = (3.14/3)*(D['Conical radius']**2)*D['Conical depth']
       
    #     D['Volume'] = D['Cylindrical volume'] + D['Conical volume']
       
    #     # Primary clarifiers can be center feed or peripheral feed. The design here is for the more
    #     # commonly deployed center feed.
       
    #     # Depth of the center feed lies between 30-75% of sidewater depth
    #     D['Center feed depth'] = 0.5*D['Cylindrical depth']
    #     # Typical conventional feed wells are designed for an average downflow velocity
    #     # of 10-13 mm/s and maximum velocity of 25-30 mm/s
    #     peak_flow_safety_factor = 2.5 # assumed based on average and maximum velocities
    #     upflow_velocity = self.upflow_velocity # in m/hr (converted from 12 mm/sec)
    #     D['Upflow velocity'] = upflow_velocity*peak_flow_safety_factor # in m/hr
    #     Center_feed_area = design_flow/D['Upflow velocity'] # in m2
    #     D['Center feed diameter'] = ((4*Center_feed_area)/3.14)**(1/2) # Sanity check: Diameter of the center feed lies between 15-25% of tank diameter

    #     # Amount of concrete required
    #     D_tank = D['Cylindrical depth']*39.37 # m to inches 
    #     # Thickness of the wall concrete, [m]. Default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
    #     thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
    #     inner_diameter = D['Cylindrical diameter']
    #     outer_diameter = inner_diameter + 2*thickness_concrete_wall
    #     volume_cylindrical_wall = (3.14*D['Cylindrical depth']/4)*(outer_diameter**2 - inner_diameter**2)
    #     volume_conical_wall = (3.14/3)*(D['Conical depth']/4)*(outer_diameter**2 - inner_diameter**2)
    #     D['Volume of concrete wall'] = volume_cylindrical_wall + volume_conical_wall # in m3
       
    #     # Amount of metal required for center feed
    #     thickness_metal_wall = 0.5 # in m (!! NEED A RELIABLE SOURCE !!)
    #     inner_diameter_center_feed = D['Center feed diameter']
    #     outer_diameter_center_feed = inner_diameter_center_feed + 2*thickness_metal_wall
    #     volume_center_feed = (3.14*D['Center feed depth']/4)*(outer_diameter_center_feed**2 - inner_diameter_center_feed **2)
    #     density_ss = 7930 # kg/m3, 18/8 Chromium
    #     D['Stainless steel'] = volume_center_feed*density_ss # in kg
       
    #     # Pumps
    #     pipe, pumps = self._design_pump()
    #     D['Pump pipe stainless steel'] = pipe
    #     D['Pump stainless steel'] = pumps
    #     # For primary clarifier 
    #     D['Number of pumps'] = D['Number of clarifiers']
       
    # def _cost(self):
       
    #     self._mixed.mix_from(self.ins)
    #     D = self.design_results
    #     C = self.baseline_purchase_costs
       
    #     # Construction of concrete and stainless steel walls
    #     C['Wall concrete'] = D['Number of clarifiers']*D['Volume of concrete wall']*self.wall_concrete_unit_cost
    #     C['Wall stainless steel'] = D['Number of clarifiers']*D['Stainless steel']*self.stainless_steel_unit_cost
        
    #     # Cost of equipment 
        
    #     # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
    #     # Scraper 
    #     # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
    #     base_cost_scraper = 2500
    #     base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
    #     clarifier_flow = self._mixed.get_total_flow('m3/hr')/D['Number of clarifiers']
    #     C['Scraper'] = D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
    #     base_power_scraper = 2.75 # in kW
    #     scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
    #     # v notch weir
    #     # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
    #     base_cost_v_notch_weir = 6888
    #     base_flow_v_notch_weir = 10 # in m3/hr
    #     C['v notch weir'] = D['Number of clarifiers']*base_cost_v_notch_weir*(clarifier_flow/base_flow_v_notch_weir)**0.6
        
    #     # Pump (construction and maintenance)
    #     pumps = self.pumps
    #     add_OPEX = self.add_OPEX
    #     pump_cost = 0.
    #     building_cost = 0.
    #     opex_o = 0.
    #     opex_m = 0.
       
    #     for i in pumps:
    #         p = getattr(self, f'{i}_pump')
    #         p_cost = p.baseline_purchase_costs
    #         p_add_opex = p.add_OPEX
    #         pump_cost += p_cost['Pump']
    #         building_cost += p_cost['Pump building']
    #         opex_o += p_add_opex['Pump operating']
    #         opex_m += p_add_opex['Pump maintenance']

    #     C['Pumps'] = pump_cost*D['Number of pumps']
    #     C['Pump building'] = building_cost*D['Number of pumps']
    #     add_OPEX['Pump operating'] = opex_o*D['Number of pumps']
    #     add_OPEX['Pump maintenance'] = opex_m*D['Number of pumps']
       
    #     # Power
    #     pumping = 0.
    #     for ID in self.pumps:
    #         p = getattr(self, f'{ID}_pump')
    #         if p is None:
    #             continue
    #         pumping += p.power_utility.rate
            
    #     pumping = pumping*D['Number of pumps']
        
    #     self.power_utility.consumption += pumping
    #     self.power_utility.consumption += scraper_power
        
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

class PrimaryClarifier(SanUnit):
    
    """
    Primary clarifier adapted from the design of thickener as defined in BSM-2. [1]
    ----------
    ID : str
        ID for the Primary Clarifier. The default is ''.
    ins : class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1. 
    outs : class:`WasteStream`
        Sludge and treated effluent.
    thickener_perc : float
        The percentage of solids in the underflow of the clarifier.[1]
    TSS_removal_perc : float
        The percentage of suspended solids removed in the clarifier.[1]
    surface_overflow_rate : float
        Surface overflow rate in the clarifier in [(m3/day)/m2]. [3]
        Design SOR value for clarifier is 41 (m3/day)/m2 if it does not receive WAS.
        Design SOR value for clarifier is 29 (m3/day)/m2 if it receives WAS.
        Typically SOR lies between 30-50 (m3/day)/m2. 
        Here default value of 41 (m3/day)/m2 is used.
    depth_clarifier : float
        Depth of clarifier. Typical depths range from 3 m to 4.9 m [2, 3]. 
        Default value of 4.5 m would be used here. 
    downward_flow_velocity : float, optional
        Speed on the basis of which center feed diameter is designed [m/hr]. [4]
        The default is 36 m/hr. (10 mm/sec)
    F_BM : dict
        Equipment bare modules.
        
    Examples
    --------
    >>> from qsdsan import set_thermo, Components, WasteStream
    >>> cmps = Components.load_default()
    >>> cmps_test = cmps.subgroup(['S_F', 'S_NH4', 'X_OHO', 'H2O'])
    >>> set_thermo(cmps_test)
    >>> ws = WasteStream('ws', S_F = 10, S_NH4 = 20, X_OHO = 15, H2O=1000)
    >>> from qsdsan.sanunits import Thickener
    >>> TC = Thickener(ID='TC', ins= (ws), outs=('sludge', 'effluent'))
    >>> TC.simulate()
    >>> sludge, effluent = TC.outs
    >>> sludge.imass['X_OHO']/ws.imass['X_OHO']
    0.98
    >>> TC.show() # doctest: +ELLIPSIS
    Thickener: TC
    ins...
    [0] ws
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1e+04
                    S_NH4  2e+04
                    X_OHO  1.5e+04
                    H2O    1e+06
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 23873.0 mg/L
         BOD        : 14963.2 mg/L
         TC         : 8298.3 mg/L
         TOC        : 8298.3 mg/L
         TN         : 20363.2 mg/L
         TP         : 367.6 mg/L
         TK         : 68.3 mg/L
    outs...
    [0] sludge
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    1.56e+03
                    S_NH4  3.11e+03
                    X_OHO  1.47e+04
                    H2O    1.56e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 95050.4 mg/L
         BOD        : 55228.4 mg/L
         TC         : 34369.6 mg/L
         TOC        : 34369.6 mg/L
         TN         : 24354.4 mg/L
         TP         : 1724.0 mg/L
         TK         : 409.8 mg/L
    [1] effluent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): S_F    8.44e+03
                    S_NH4  1.69e+04
                    X_OHO  300
                    H2O    8.44e+05
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 9978.2 mg/L
         BOD        : 7102.9 mg/L
         TC         : 3208.8 mg/L
         TOC        : 3208.8 mg/L
         TN         : 19584.1 mg/L
         TP         : 102.9 mg/L
         TK         : 1.6 mg/L

    References
    ----------
    .. [1] Gernaey, Krist V., Ulf Jeppsson, Peter A. Vanrolleghem, and John B. Copp.
    Benchmarking of control strategies for wastewater treatment plants. IWA publishing, 2014.
    [2] Metcalf, Leonard, Harrison P. Eddy, and Georg Tchobanoglous. Wastewater 
    engineering: treatment, disposal, and reuse. Vol. 4. New York: McGraw-Hill, 1991.
    [3] Chapter-10: Primary Treatment. Design of water resource recovery facilities. 
    WEF Manual of Practice No. 8. 6th Edition. Virginia: McGraw-Hill, 2018. 
    [4] Introduction to Wastewater Clarifier Design by Nikolay Voutchkov, PE, BCEE.
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    # Costs
    wall_concrete_unit_cost = 1081.73 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    slab_concrete_unit_cost = 582.48 # $/m3 (Hydromantis. CapdetWorks 4.0. https://www.hydromantis.com/CapdetWorks.html)
    stainless_steel_unit_cost=1.8 # Alibaba. Brushed Stainless Steel Plate 304. https://www.alibaba.com/product-detail/brushed-stainless-steel-plate-304l-stainless_1600391656401.html?spm=a2700.details.0.0.230e67e6IKwwFd
    
    pumps = ('sludge',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', thickener_perc=7, 
                  TSS_removal_perc=98, surface_overflow_rate = 41, depth_clarifier=4.5, 
                  downward_flow_velocity=36, F_BM=default_F_BM, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with)
        self.thickener_perc = thickener_perc 
        self.TSS_removal_perc = TSS_removal_perc
        self.surface_overflow_rate = surface_overflow_rate
        self.depth_clarifier = depth_clarifier
        self.downward_flow_velocity = downward_flow_velocity
        self.F_BM.update(F_BM)
        self._mixed = WasteStream(f'{ID}_mixed')        
        self._sludge = self.outs[0].copy(f'{ID}_sludge')
        
    @property
    def thickener_perc(self):
        '''tp is the percentage of Suspended Sludge in the underflow of the clarifier'''
        return self._tp

    @thickener_perc.setter
    def thickener_perc(self, tp):
        if tp is not None:
            if tp>=100 or tp<=0:
                raise ValueError(f'should be between 0 and 100 not {tp}')
            self._tp = tp
        else: 
            raise ValueError('percentage of SS in the underflow of the thickener expected from user')
            
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
            
    @property
    def TSS_removal_perc(self):
        '''The percentage of suspended solids removed in the clarifier'''
        return self._TSS_rmv

    @TSS_removal_perc.setter
    def TSS_removal_perc(self, TSS_rmv):
        if TSS_rmv is not None:
            if TSS_rmv>=100 or TSS_rmv<=0:
                raise ValueError(f'should be between 0 and 100 not {TSS_rmv}')
            self._TSS_rmv = TSS_rmv
        else: 
            raise ValueError('percentage of suspended solids removed in the clarifier expected from user')
            
    @property
    def thickener_factor(self):
        self._mixed.mix_from(self.ins)
        inf = self._mixed
        _cal_thickener_factor = self._cal_thickener_factor
        if not self.ins: return
        elif inf.isempty(): return
        else: 
            TSS_in = inf.get_TSS()
            thickener_factor = _cal_thickener_factor(TSS_in)
        return thickener_factor
    
    @property
    def thinning_factor(self):
        self._mixed.mix_from(self.ins)
        inf = self._mixed
        TSS_in = inf.get_TSS()
        _cal_thickener_factor = self._cal_thickener_factor
        thickener_factor = _cal_thickener_factor(TSS_in)
        _cal_parameters = self._cal_parameters
        Qu_factor, thinning_factor = _cal_parameters(thickener_factor)
        return thinning_factor
    
    def _cal_thickener_factor(self, TSS_in):
        if TSS_in > 0:
            thickener_factor = self._tp*10000/TSS_in
            if thickener_factor<1:
                thickener_factor=1
            return thickener_factor
        else: return None
            
    def _cal_parameters(self, thickener_factor):
        if thickener_factor<1:
            Qu_factor = 1
            thinning_factor=0
        else:
            Qu_factor = self._TSS_rmv/(100*thickener_factor)
            thinning_factor = (1 - (self._TSS_rmv/100))/(1 - Qu_factor)
        return Qu_factor, thinning_factor
    
    def _update_parameters(self):
        
        # Thickener_factor, Thinning_factor, and Qu_factor need to be 
        # updated again and again. while dynamic simulations 
        
        cmps = self.components 
    
        TSS_in = np.sum(self._state[:-1]*cmps.i_mass*cmps.x)
        _cal_thickener_factor = self._cal_thickener_factor
        self.updated_thickener_factor = _cal_thickener_factor(TSS_in)
        _cal_parameters = self._cal_parameters
        
        updated_thickener_factor = self.updated_thickener_factor
        self.updated_Qu_factor, self.updated_thinning_factor = _cal_parameters(updated_thickener_factor)
        
    def _run(self):
        self._mixed.mix_from(self.ins)
        inf = self._mixed
        sludge, eff = self.outs
        cmps = self.components
        
        TSS_rmv = self._TSS_rmv
        thinning_factor = self.thinning_factor
        thickener_factor = self.thickener_factor
        
        # The following are splits by mass of particulates and solubles 
        
        # Note: (1 - thinning_factor)/(thickener_factor - thinning_factor) = Qu_factor
        Zs = (1 - thinning_factor)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        Ze = (thickener_factor - 1)/(thickener_factor - thinning_factor)*inf.mass*cmps.s
        
        Xe = (1 - TSS_rmv/100)*inf.mass*cmps.x
        Xs = (TSS_rmv/100)*inf.mass*cmps.x
        
        # e stands for effluent, s stands for sludge 
        Ce = Ze + Xe 
        Cs = Zs + Xs
    
        eff.set_flow(Ce,'kg/hr')
        sludge.set_flow(Cs,'kg/hr')
       
    def _init_state(self):
       
        # This function is run only once during dynamic simulations 
    
        # Since there could be multiple influents, the state of the unit is 
        # obtained assuming perfect mixing 
        Qs = self._ins_QC[:,-1]
        Cs = self._ins_QC[:,:-1]
        self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.
        
        # To initialize the updated_thickener_factor, updated_thinning_factor
        # and updated_Qu_factor for dynamic simulation 
        self._update_parameters()
        
    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Thickener''' 
        
        # This function is run multiple times during dynamic simulation 
        
        # Remember that here we are updating the state array of size n, which is made up 
        # of component concentrations in the first (n-1) cells and the last cell is flowrate. 
        
        # So, while in the run function the effluent and sludge are split by mass, 
        # here they are split by concentration. Therefore, the split factors are different. 
        
        # Updated intrinsic modelling parameters are used for dynamic simulation 
        thickener_factor = self.updated_thickener_factor
        thinning_factor = self.updated_thinning_factor
        Qu_factor = self.updated_Qu_factor
        cmps = self.components
        
        # For sludge, the particulate concentrations are multiplied by thickener factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        uf, of = self.outs
        if uf.state is None: uf.state = np.zeros(len(cmps)+1)
        uf.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thickener_factor
        uf.state[-1] = self._state[-1]*Qu_factor
        
        # For effluent, the particulate concentrations are multiplied by thinning factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        if of.state is None: of.state = np.zeros(len(cmps)+1)
        of.state[:-1] = self._state[:-1]*cmps.s*1 + self._state[:-1]*cmps.x*thinning_factor
        of.state[-1] = self._state[-1]*(1 - Qu_factor)

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Thickener'''
        
        # This function is run multiple times during dynamic simulation 
        
        # Remember that here we are updating the state array of size n, which is made up 
        # of component concentrations in the first (n-1) cells and the last cell is flowrate. 
        
        # So, while in the run function the effluent and sludge are split by mass, 
        # here they are split by concentration. Therefore, the split factors are different. 
        
        # Updated intrinsic modelling parameters are used for dynamic simulation
        thickener_factor = self.updated_thickener_factor
        thinning_factor = self.updated_thinning_factor
        Qu_factor = self.updated_Qu_factor
        cmps = self.components
        
        # For sludge, the particulate concentrations are multiplied by thickener factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same. 
        uf, of = self.outs
        if uf.dstate is None: uf.dstate = np.zeros(len(cmps)+1)
        uf.dstate[:-1] = self._dstate[:-1]*cmps.s*1 + self._dstate[:-1]*cmps.x*thickener_factor
        uf.dstate[-1] = self._dstate[-1]*Qu_factor
        
        # For effluent, the particulate concentrations are multiplied by thinning factor, and
        # flowrate is multiplied by Qu_factor. The soluble concentrations remains same.
        if of.dstate is None: of.dstate = np.zeros(len(cmps)+1)
        of.dstate[:-1] = self._dstate[:-1]*cmps.s*1 + self._dstate[:-1]*cmps.x*thinning_factor
        of.dstate[-1] = self._dstate[-1]*(1 - Qu_factor)
     
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        
        # This function is run multiple times during dynamic simulation 
        
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        _update_parameters = self._update_parameters
        def yt(t, QC_ins, dQC_ins):
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            dQ_ins = dQC_ins[:, -1]
            dC_ins = dQC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            Q_dot = dQ_ins.sum()
            C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            _dstate[-1] = Q_dot
            _dstate[:-1] = C_dot
    
            _update_parameters()
            _update_state()
            _update_dstate()
        self._AE = yt
    
    def _design_pump(self):
        ID, pumps = self.ID, self.pumps
        self._sludge.copy_like(self.outs[0])
        sludge = self._sludge
        
        ins_dct = {
            'sludge': sludge,
            }
        
        type_dct = dict.fromkeys(pumps, 'sludge')
        inputs_dct = dict.fromkeys(pumps, (1,))
        
        D = self.design_results
        influent_Q = sludge.get_total_flow('m3/hr')/D['Number of clarifiers']
        influent_Q_mgd = influent_Q*0.00634 # m3/hr to MGD
       
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                ID = f'{ID}_{i}'
                capacity_factor=1
                pump = WWTpump(
                    ID=ID, ins= ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=influent_Q_mgd, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)

        pipe_ss, pump_ss = 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p.simulate()
            p_design = p.design_results
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
        return pipe_ss, pump_ss
    
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
        'Pump pipe stainless steel' : 'kg',
        'Pump stainless steel': 'kg',
        'Number of pumps': 'ea'
    }
    
    def _design(self):
        
        self._mixed.mix_from(self.ins)
        mixed = self._mixed
        D = self.design_results
        
        # Number of clarifiers based on tentative suggestions by Jeremy 
        # (would be verified through collaboration with industry)
        total_flow = (mixed.get_total_flow('m3/hr')*24)/3785 # in MGD
        if total_flow <= 3:
            D['Number of clarifiers'] = 2
        elif total_flow > 3 and total_flow <= 8:
            D['Number of clarifiers'] = 3
        elif total_flow > 8 and total_flow <=20:
            D['Number of clarifiers'] = 4
        else:
            D['Number of clarifiers'] = 4
            total_flow -= 20
            D['Number of clarifiers'] += np.ceil(total_flow/20)
            
        D['SOR'] = self.surface_overflow_rate # in (m3/day)/m2
        D['Volumetric flow'] =  (mixed.get_total_flow('m3/hr')*24)/D['Number of clarifiers'] # m3/day
        D['Surface area'] = D['Volumetric flow']/D['SOR'] # in m2
        D['Cylindrical diameter'] = np.sqrt(4*D['Surface area']/np.pi) #in m
        
        #Check on cylindrical diameter [2, 3]
        if D['Cylindrical diameter'] < 3 or D['Cylindrical diameter'] > 60:
            Cylindrical_dia = D['Cylindrical diameter']
            warn(f'Cylindrical diameter = {Cylindrical_dia} is not between 3 m and 60 m')
        
        D['Conical radius'] = D['Cylindrical diameter']/2
        # The slope of the bottom conical floor lies between 1:10 to 1:12 [3, 4]
        D['Conical depth'] = D['Conical radius']/12 
        D['Clarifier depth'] = self.depth_clarifier #in m 
        D['Cylindrical depth'] = D['Clarifier depth'] -  D['Conical depth']
        
        # Check on cylindrical and conical depths 
        if D['Cylindrical depth'] < D['Conical depth']:
            Cylindrical_depth = D['Cylindrical depth']
            Conical_depth = D['Conical depth']
            warn(f'Cylindrical depth = {Cylindrical_depth} is lower than Conical depth = {Conical_depth}')
        
        D['Cylindrical volume'] = np.pi*np.square(D['Cylindrical diameter']/2)*D['Cylindrical depth'] #in m3
        D['Conical volume'] = (3.14/3)*(D['Conical radius']**2)*D['Conical depth'] #in m3
        D['Volume'] = D['Cylindrical volume'] + D['Conical volume'] #in m3
        
        D['Hydraulic Retention Time'] = D['Volume']/(D['Volumetric flow']/24) #in hrs
        
        # Check on cylinderical HRT [3]
        if D['Hydraulic Retention Time'] < 1.5 or D['Hydraulic Retention Time'] > 2.5:
            HRT = D['Hydraulic Retention Time']
            warn(f'HRT = {HRT} is not between 1.5 and 2.5 hrs')
        
        # The design here is for center feed of clarifier.
        
        # Depth of the center feed lies between 30-75% of sidewater depth. [3, 4]
        D['Center feed depth'] = 0.5*D['Cylindrical depth']
        # Typical conventional feed wells are designed for an average downflow velocity
        # of 10-13 mm/s and maximum velocity of 25-30 mm/s. [4]
        peak_flow_safety_factor = 2.5 # assumed based on average and maximum velocities
        D['Downward flow velocity'] = self.downward_flow_velocity*peak_flow_safety_factor # in m/hr
        
        Center_feed_area = (D['Volumetric flow']/24)/D['Downward flow velocity'] # in m2
        
        D['Center feed diameter'] = np.sqrt(4*Center_feed_area/np.pi) 

        #Sanity check: Diameter of the center feed lies between 15-25% of tank diameter [4]
        #The lower limit of this check has been modified to 10% based on advised range of down flow velocity in [4]. 
        if D['Center feed diameter'] < 0.10*D['Cylindrical diameter'] or D['Center feed diameter']  > 0.25*D['Cylindrical diameter']:
            cf_dia = D['Center feed diameter'] 
            tank_dia = D['Cylindrical diameter']
            warn(f'Diameter of the center feed does not lie between 15-25% of tank diameter. It is {cf_dia*100/tank_dia}% of tank diameter')

        # Amount of concrete required
        D_tank = D['Cylindrical depth']*39.37 # m to inches 
        # Thickness of the wall concrete [m]. Default to be minimum of 1 feet with 1 inch added for every feet of depth over 12 feet.
        thickness_concrete_wall = (1 + max(D_tank-12, 0)/12)*0.3048 # from feet to m
        inner_diameter = D['Cylindrical diameter']
        outer_diameter = inner_diameter + 2*thickness_concrete_wall
        volume_cylindercal_wall = (np.pi*D['Cylindrical depth']/4)*(outer_diameter**2 - inner_diameter**2)
        D['Volume of concrete wall'] = volume_cylindercal_wall # in m3
        
        # Concrete slab thickness, [ft], default to be 2 in thicker than the wall thickness. (Brian's code)
        thickness_concrete_slab = thickness_concrete_wall + (2/12)*0.3048 # from inch to m
        outer_diameter_cone = inner_diameter + 2*(thickness_concrete_wall + thickness_concrete_slab)
        volume_conical_wall = (np.pi/(3*4))*(((D['Conical depth'] + thickness_concrete_wall + thickness_concrete_slab)*(outer_diameter_cone**2)) - (D['Conical depth']*(inner_diameter)**2))
        D['Volume of concrete slab'] = volume_conical_wall
        
        # Amount of metal required for center feed
        thickness_metal_wall = 0.3048 # equal to 1 feet, in m (!! NEED A RELIABLE SOURCE !!)
        inner_diameter_center_feed = D['Center feed diameter']
        outer_diameter_center_feed = inner_diameter_center_feed + 2*thickness_metal_wall
        volume_center_feed = (3.14*D['Center feed depth']/4)*(outer_diameter_center_feed**2 - inner_diameter_center_feed**2)
        density_ss = 7930 # kg/m3, 18/8 Chromium
        D['Stainless steel'] = volume_center_feed*density_ss # in kg
       
        # Pumps
        pipe, pumps = self._design_pump()
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        
        #For primary clarifier 
        D['Number of pumps'] = D['Number of clarifiers']
        
    def _cost(self):
        
        self._mixed.mix_from(self.ins)
       
        D = self.design_results
        C = self.baseline_purchase_costs
       
        # Construction of concrete and stainless steel walls
        C['Wall concrete'] = D['Number of clarifiers']*D['Volume of concrete wall']*self.wall_concrete_unit_cost
        
        C['Slab concrete'] = D['Number of clarifiers']*D['Volume of concrete slab']*self.slab_concrete_unit_cost
        
        C['Wall stainless steel'] = D['Number of clarifiers']*D['Stainless steel']*self.stainless_steel_unit_cost
       
        # Cost of equipment 
        
        # Source of scaling exponents: Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol by NREL.
        
        # Scraper 
        # Source: https://www.alibaba.com/product-detail/Peripheral-driving-clarifier-mud-scraper-waste_1600891102019.html?spm=a2700.details.0.0.47ab45a4TP0DLb
        base_cost_scraper = 2500
        base_flow_scraper = 1 # in m3/hr (!!! Need to know whether this is for solids or influent !!!)
        clarifier_flow = D['Volumetric flow']/24
        C['Scraper'] =  D['Number of clarifiers']*base_cost_scraper*(clarifier_flow/base_flow_scraper)**0.6
        base_power_scraper = 2.75 # in kW
        scraper_power = D['Number of clarifiers']*base_power_scraper*(clarifier_flow/base_flow_scraper)**0.6
        
        # v notch weir
        # Source: https://www.alibaba.com/product-detail/50mm-Tube-Settler-Media-Modules-Inclined_1600835845218.html?spm=a2700.galleryofferlist.normal_offer.d_title.69135ff6o4kFPb
        base_cost_v_notch_weir = 6888
        base_flow_v_notch_weir = 10 # in m3/hr
        C['v notch weir'] = D['Number of clarifiers']*base_cost_v_notch_weir*(clarifier_flow/base_flow_v_notch_weir)**0.6
       
        # Pump (construction and maintainance)
        pumps = self.pumps
        add_OPEX = self.add_OPEX
        pump_cost = 0.
        building_cost = 0.
        opex_o = 0.
        opex_m = 0.
       
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            p_cost = p.baseline_purchase_costs
            p_add_opex = p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost*D['Number of clarifiers']
        C['Pump building'] = building_cost*D['Number of clarifiers']
        add_OPEX['Pump operating'] = opex_o*D['Number of clarifiers']
        add_OPEX['Pump maintenance'] = opex_m*D['Number of clarifiers']
       
        # Power
        pumping = 0.
        for ID in self.pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate
        
        pumping = pumping*D['Number of clarifiers']
        
        self.power_utility.rate += pumping
        self.power_utility.rate += scraper_power
