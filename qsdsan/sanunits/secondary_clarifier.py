from numpy import maximum as npmax, minimum as npmin, exp as npexp
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
                    ID=ID, ins=ins_dct[i], thermo = self.thermo, pump_type=type_dct[i],
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