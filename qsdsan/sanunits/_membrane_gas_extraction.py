# -*- coding: utf-8 -*-
"""
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Ian Song <song0231@umn.edu>
    
    Saumitra Rai <raisaumitra9@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
"""

from qsdsan import SanUnit
import numpy as np

__all__ = ('GasExtractionMembrane', 'MembraneGasExtraction',)

class GasExtractionMembrane(SanUnit):
    
    """
    Gas Extraction Membrane 

    Parameters
    ----------
    ID : str
        ID for the Gas Extraction Membrane. The default is 'GEM'.
    ins : class:`WasteStream`
        Influent to the Gas Extraction Membrane. Expected number of influent is 3. 
    outs : class:`WasteStream`
        Gas and Liquid streams are expected effluents.
    FiberID : float
        Inner diameter of the membrane [m].
    FiberOD : float
        Outer diameter of the membrane [m].
    NumTubes : float
        The number of fibers in the membrane.  
    ShellDia : float
        The diameter of the shell [m]. 
    SurfArea : float
        Surface area of membrane [m^2]. 
    GasID : array
        Array containing IDs of gases to be extracted. 
    PVac : float
        Operating vaccum pressure in the membrane [Pa].
    segs : float
        Number of segments considered in the membrane. 
    GasPerm : dict
        Dictionary of permeability of gases.
    HenryPreFac : dict
        Dictionary of Henry's Law Factor for gases. 
    HenrySlope : dict
        Dictionary of Henry's Slope for gases. 
    WilkeChang : dict
        Dictionary of Wilke Chang correlation for gases. 
    """
    
    _N_ins = 1
    _N_outs = 2
    
    # All gas properties are in form of dictionaries
    
    _GasPerm = {
        'H2': 650*(3.35e-16), 
        'O2': 600*(3.35e-16),
        'N2': 280*(3.35e-16), 
        'CO2': 3250*(3.35e-16), 
        'CH4': 950*(3.35e-16),
        'H2O': 36000*(3.35e-16)
        }
    
    _HenryPreFac = {
        'H2': 7.8e-6, 
        'O2': 1.2e-5, 
        'N2': 6e-6, 
        'CO2' : 3.5e-4, 
        'CH4': 1.3e-5,
        'H2O': 1
        }
    
    _HenrySlope = {
        'H2': 640, 
        'O2': 1800, 
        'N2': 1300, 
        'CO2' : 2600, 
        'CH4' : 1900,
        'H2O': 1
        }
    
    _WilkeChang = {
        'H2': 9.84, 
        'O2': 1.90,
        'N2': 1.77, 
        'CO2': 2.6, 
        'CH4': 2.2, 
        'H2O': 1
        }
    
    # Constructor: Initialize the instance variables
    def __init__(self, ID='GEM', ins=None, outs=(), thermo=None, isdynamic=True, 
                  init_with='WasteStream', F_BM_default=None,   FiberID=190e-6, 
                  FiberOD=300e-6, NumTubes=1512, ShellDia=1.89e-2, SurfArea=0.1199,   
                  GasID = ['H2', 'O2', 'N2', 'CO2', 'CH4', 'H2O'], PVac = 97.325, 
                  segs = 50, GasPerm = {}, HenryPreFac = {}, HenrySlope = {}, 
                  WilkeChang = {}):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, 
                         init_with=init_with, F_BM_default=F_BM_default)
        self.FiberID = FiberID        # Fiber Inner Diameter [m]
        self.FiberOD = FiberOD        # Fiber Outer Diameter [m]
        self.MemThick = (FiberOD - FiberID)/2 # Membrane Thickness [m]
        self.NumTubes = NumTubes      # Number of Tubes []
        self.ShellDia = ShellDia      # Shell Diameter [m]
        self.SurfArea = SurfArea      # Surface Area [m^2]
        self.GasID = GasID            # IDs of gas used in the process
        self.PVac = PVac              # Operating Vacuum Pressure [-kPa]
        self.segs = segs              # Number of segments 
        #self.Volume = VolBatchTank   # Volume of the bioreactor (Don't think this is needed)
                
        self.indexer = GasID.index
        dct_gas_perm = GasPerm or self._GasPerm
        self.set_GasPerm(**dct_gas_perm)
        dct_gas_hpf =  HenryPreFac or self._HenryPreFac
        self.set_HenryPreFac(**dct_gas_hpf)
        dct_gas_hs = HenrySlope or self._HenrySlope
        self.set_HenrySlope(**dct_gas_hs)
        dct_gas_wc = WilkeChang or self._WilkeChang
        self.set_WilkeChang(**dct_gas_wc)
        
        cmps = self.thermo.chemicals
        # self.indexer = cmps.index
        # self.idx ensures that the indexing in further code is only for gases 
        # and not all components in the influent
        self.idx = cmps.indices(self.GasID) 
        # to save the index of water in the array of gases, to be used later 
        # for i, ID in enumerate(GasID):
        #     if ID == 'H2O': 
        #         self.h2o_j = i
        #         break
        #!!! alternatively    
        self.h2o_j = GasID.index('H2O')
        self.gas_mass2mol = (cmps.i_mass/cmps.chem_MW)[self.idx]
    
    @property 
    def FiberOD(self):
        return self._FiberOD
    
    @FiberOD.setter
    def FiberOD(self, FiberOD):
        if FiberOD is not None:
            self._FiberOD = FiberOD
        else:
            raise ValueError('FiberOD expected from user')
            
    @property 
    def FiberID(self):
        return self._FiberID
        
    @FiberID.setter
    def FiberID(self, FiberID):
        if FiberID is not None:
            self._FiberID = FiberID
        else:
            raise ValueError('Inner Diameter of fiber expected from user')
            
    @property 
    def NumTubes(self):
        return self._NumTubes
        
    @NumTubes.setter
    def NumTubes(self, NumTubes):
        if NumTubes is not None:
            self._NumTubes = NumTubes
        else:
            raise ValueError('Number of tubes expected from user')
    
    @property 
    def ShellDia(self):
        return self._ShellDia
        
    @ShellDia.setter
    def ShellDia(self, ShellDia):
        if ShellDia is not None:
            self._ShellDia = ShellDia
        else:
            raise ValueError('Diameter of the shell expected from user')
            
    @property 
    def SurfArea(self):
        return self._SurfArea
        
    @SurfArea.setter
    def SurfArea(self, SurfArea):
        if SurfArea is not None:
            self._SurfArea = SurfArea
        else:
            raise ValueError('Surface Area of Membrane expected from user')
        
    # Calculate the volume fraction of the lumen to the shell. 
    @property
    def VolFrac(self):
        lumenVol = self.NumTubes*np.pi*((self.FiberID/2)**2)
        shellVol = (np.pi*((self.ShellDia/2)**2)) - self.NumTubes*np.pi*((self.FiberOD/2)**2)
        return lumenVol/shellVol

    # Calculate the effective length of the membrane by taking the ratio of the 
    # declared surface area, and dividing it by the area per length of the tube. 
    @property
    def Length(self):
        memSurf = self.NumTubes*np.pi*self.FiberOD
        return self.SurfArea/memSurf

    # Determine the Shell cross-sectional area
    @property
    def ShellAc(self):
        return np.pi*((self.ShellDia/2)**2) - self.NumTubes*np.pi*((self.FiberOD/2)**2)
    
    @property
    def HeL(self):
        # Define the constant properties of gas
        TRefH = 298.15 # Reference T for Henry's Law [K]
        NIST_HeL = self._hpf*(np.exp(self._hs*(1/self.ins[0].T - 1/TRefH)))
        return 1/NIST_HeL
    
    @property
    def Diff(self):
        inf, = self.ins
        cmps = inf.components
        self._Vc = np.array([cmp.Vc for cmp in cmps])
        Phi = self._wc # Wilke Chang Correlation Factor 
        
        # Define the constant properties of gas
        TRefMu = 300                # Reference T for H2O Viscosity [K]
        MWH2O = cmps.H2O.MW         # Molar Weight of the Solvent [Da]

        # Reduced Temp for Water Viscosity
        Tb = self.ins[0].T/TRefMu 

        # Temperature Dependent Water Viscosity [cP]
        mu = (1*10**(-3))*(280.68*(Tb**(-1.9)) + 511.45*(Tb**(-7.7)) + 61.131*(Tb**(-19.6)) + 0.45903*(Tb**(-40)) ) 

        # Molar Volume at Normal Boiling Point [cm^3/mol]
        # Yoel: Critical models are not very reliable
        # Yoel: Eqns of State or other models already available in QSDsan can be looked into
        V1 = 0.285*(self._Vc*1000000)**(1.048) # unit conversion from m^3/mol (QSDsan) to cm^3/mol (here)
        # Diffusion Coefficient [m^2/s]
        D = 0.0001*(7.4*10**(-8))*np.sqrt(MWH2O*Phi)*self.ins[0].T/(mu*V1**(0.6))
        return D    
    
    def set_GasPerm(self, **kwargs):
        self.set_prop('_gasp', **kwargs)
            
    def set_WilkeChang(self, **kwargs):
        self.set_prop('_wc', **kwargs)
            
    def set_HenryPreFac(self, **kwargs):
        self.set_prop('_hpf', **kwargs)
            
    def set_HenrySlope(self, **kwargs):
        self.set_prop('_hs', **kwargs)
            
    def set_prop(self, attr_name, **kwargs):
        idxr = self.indexer
        try: attr = getattr(self, attr_name)
        except: attr = self.__dict__[attr_name] = np.zeros(len(self.chemicals))
        for k, v in kwargs.items():
            attr[idxr(k)] = v
            
    def _init_state(self):
        # inf, = self.ins
        # cmps = inf.components
        # C = self._ins_QC[0,:-1]/cmps.chem_MW*cmps.i_mass
        # Cs = C[self.idx] #idx selects only gases 
        # Seg = self.segs
        numGas = len(self.GasID)
        # self._state = np.zeros(2*Seg*numGas)
        # for i in range(0, 2*Seg*numGas, 2*numGas):
        #     for j in range(numGas):
        #         self._state[j+i] = Cs[j]
        seg_i = np.zeros(2*numGas)
        seg_i[:numGas] = self._ins_QC[0,self.idx]*self.gas_mass2mol
        self._state = np.tile(seg_i, self.segs)
        self._dstate = self._state*0

    def _update_state(self):
        inf, = self.ins
        gas, liq = self.outs
        idx = self.idx  
        numGas = len(self.GasID)
        mass2mol = self.gas_mass2mol
        y = self._state
        # Need to add effluent streams for liquid (this includes all cmps of influent ws) and gas.
        # The multiplication of any of the first n-1 array element with last element should give out g/day values.
        
        if gas.state is None:
            cmps = inf.components
            gas.state = np.zeros(len(cmps) + 1)
            liq.state = np.zeros(len(cmps) + 1)
            
        liq.state[:] = inf.state
        liq.state[idx] = y[-2*numGas: -numGas]/mass2mol
        
        gas.state[-1] = 1
        gas.state[idx] = (y[:numGas] - y[-2*numGas: -numGas])/mass2mol * liq.state[-1]
        
        # The of the effluent gas in extraction membrane is the difference 
        # between lumen concentration in the last and first segment
        #!!! why? It seems this only holds when the unit is at steady state
        # gas_state_in_unit = y[:numGas] - y[-2*numGas: -numGas]  # in mol/m3
        # Molar_flow_gases = self._ins_QC[0,-1]*gas_state_in_unit # (m3/day)*(mol/m3) = mol/day
        # Mass_flow_gases = Molar_flow_gases*cmps.chem_MW[idx] #(mol/day)*(g/mol) = (g/day)
        
        # self._outs[0].state[idx] =  Mass_flow_gases # (g/day)
        # self._outs[0].state[-1] = 1 #(So the mutiplication with Q would give out g/day values)
        
        # # The state of effluent Liquid stream is simply the concentration of 
        # # the last lumen segment in the extraction membrane 
        # liquid_state_in_unit = y[-2*numGas: -numGas]  # in mol/m3
        # liquid_state_in_unit = (liquid_state_in_unit*cmps.chem_MW[idx])/cmps.i_mass[idx] # (mol/m3)*(g/mol) = g/m3 = mg/l
        
        # self._outs[1].state[:] = self._ins_QC[0]
        # self._outs[1].state[idx] = liquid_state_in_unit
        
        
    def _update_dstate(self):        
        inf, = self.ins
        gas, liq = self.outs
        numGas = len(self.GasID)
        mass2mol = self.gas_mass2mol
        idx = self.idx
        dy = self._dstate
        
        if gas.dstate is None:
            cmps = inf.components
            gas.dstate = np.zeros(len(cmps) + 1)
            liq.dstate = np.zeros(len(cmps) + 1)
            
        liq.dstate[:] = inf.dstate
        liq.dstate[idx] = dy[-2*numGas: -numGas]/mass2mol
        
        #!!! this is probably wrong
        gas.dstate[idx] = (dy[:numGas] - dy[-2*numGas: -numGas])/mass2mol * liq.dstate[-1]
        
        # self._outs[0].dstate = np.zeros(len(cmps) + 1)
        # # The of the effluent gas in extraction membrane is the difference 
        # # between lumen concentration in the last and first segment
        # gas_dstate_in_unit =  self._dstate[ :numGas] - self._dstate[ -2*numGas: -numGas]# in mol/m3
        # Molar_dflow_gases = self._ins_dQC[0,-1]*gas_dstate_in_unit # (m3/day)*(mol/m3) = mol/day
        # Mass_dflow_gases = Molar_dflow_gases*cmps.chem_MW[self.idx] #(mol/day)*(g/mol) = (g/day)
        
        # self._outs[0].dstate[idx] =  Mass_dflow_gases # (g/day)
        # self._outs[0].dstate[-1] = 0 # Just differentiating constant 1 to 0 
        
        # self._outs[1].dstate = np.zeros(len(cmps) + 1)
        # # The state of effluent Liquid stream is simply the concentration of 
        # # the last lumen segment in the extraction membrane 
        # liquid_dstate_in_unit = self._dstate[-2*numGas: -numGas]  # in mol/m3
        # liquid_dstate_in_unit = (liquid_dstate_in_unit*cmps.chem_MW[self.idx])/cmps.i_mass[self.idx] # (mol/m3)*(g/mol) = g/m3 = mg/l
        
        # self._outs[1].dstate = self._ins_dQC[0]
        # self._outs[1].dstate[idx] = liquid_dstate_in_unit

    def _run(self):
        s_in, = self.ins
        gas, liq = self.outs
        gas.phase = 'g'
        liq.copy_like(s_in)
        
    
    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE    
    
    def _compile_ODE(self):
        # Synthesizes the ODEs to simulate a batch reactor with side-flow gas extraction. The code takes in an object of class Membrane (Mem) and an array of objects of class Gas (GasVec). It also takes in an array of experimental conditions ExpCond. 
        
        # Extract Operating Parameters from ExpCond
        # Q = self.ins[0].F_vol  # Volumetric Flowrate [m3/sec]
        T = self.ins[0].T  # Temperature [K]
        P = self.PVac*1000 # Vacuum Pressure [Pa]
        #V = self.Volume  # Volume of the Batch Tank [L]

        idx = self.idx
        h2o_j = self.h2o_j
        mass2mol = self.gas_mass2mol
        
        # Calculate vapor pressure of water at operating temperature
        TCel = T-273.15 # Temperature [C]
        PVapH2O = np.exp(34.494-(4924.99/(TCel+237.1)))/((TCel+105)**1.57)   # Saturated Vapor Pressure of Water [Pa]

        # Define Constants
        R = 8.314       # Universal Gas Constant [J/K mol]
        
        # Extract Membrane Properties
        D = self.FiberID            # Membrane Fiber ID [m]
        l = self.MemThick           # Membrane Thickness [m]
        num_tubes = self.NumTubes   # Number of Tubes in the Module []
        L = self.Length             # Membrane Length [m]
        Segs = self.segs        # Number of segments? Ask Ian
        vFrac = self.VolFrac        # Lumen/Shell Volume Fraction [m^3/m^3]

        # Pre-allocate vectors for gas thermophysical properties
        numGas = len(self.GasID)
        # numVec = 2*numGas

        inf, = self.ins
        # cmps = inf.components
        # Extract Gas Properties
        #for i in range(0,numGas):
        
        #Diff[i] = GasVec[i].Diff()
        Diff = self.Diff
            
        #Perm_SI[i] = Mem.PermDict[GasVec[i].Name]
        Perm_SI = self._gasp 
            
        #H[i] = GasVec[i].HeL()
        H = self.HeL

        # Diff = np.array([0.0265e-7, 0.0296e-7, 0.0253e-7, 0.3199e-7])
        # Calculate dx and u
        dx = L/Segs # Length of Segments [m]
        # u = Q/((np.pi*D**2/4)*num_tubes)    # Linear Flow Velocity [m/s]
        cross_section_A = ((np.pi*D**2/4)*num_tubes)

        # Calculate the Kinematic Viscosity of Water
        Tb = T/300  # Reduced Temperature []
        mu = (1e-6)*(280.68*(Tb**(-1.9)) + 511.45*(Tb**(-7.7)) + 61.131*(Tb**(-19.6)) + 0.45903*(Tb**(-40))) # Absolute Viscosity of Water [Pa s]
        rho = (-13.851 + 0.64038*T - 1.9124e-3*T**2 + 1.8211e-6*T**3)*18   # Liquid Density of Water [kg/m^3] from DIPPR
        nu = mu/rho # Kinematic Viscosity of Water [m^2/s]

        #!!! These variables should be calculated within ODE because it is 
        #!!! dependent on the influent Q, which could change during simulation
        #!!! unless we can assume this change is negligible
        
        # # Calculate the dimensionless numbers
        # # Reynolds
        # Re = u*D/nu
        # Schmidt
        Sc = nu/Diff
        # # Sherwood
        # Sh = 1.615*(Re*Sc*D/L)**(1/3)
        
        # # Calculate Mass Transfer Coefficients
        KMem = Perm_SI/l
        # KLiq = Sh*Diff/D
        # KTot = 1/(1/KLiq + 1/(KMem*H))
        # # for j in range(0, numGas):
        # #     #if GasVec[j].Name == 'H2O':
        # #     if j == self.h2o_j:
        # #         KTot[j] = KMem[j]
        # #!!! alternatively
        # KTot[self.h2o_j] = KMem[self.h2o_j]
        
        # Initialize
        # C = self._state
        dC_lumen = np.zeros((Segs, numGas))
        dC_shell = dC_lumen.copy()

        sumCp_init = P/(R*T)
        # sumCp_fin = np.zeros(Segs)
        
        # C = self._ins_QC[0,:-1]/cmps.chem_MW*cmps.i_mass # conc. in mol/m^3 as defined by Ian 
        # Cs = C[self.idx] #self.idx ensures its only for gases 
        
        dC = self._dstate
        _update_dstate = self._update_dstate

            
        def dy_dt(t, QC_ins, QC, dQC_ins):
            # QC is exactly 'the state' as we define in _init_
            # C = QC
            Q = QC_ins[0,-1]/24/3600
            C_in = QC_ins[0, idx] * mass2mol  # mol/m^3
            u = Q/cross_section_A
            
            # Calculate the dimensionless numbers
            # Reynolds
            Re = u*D/nu
            # Sherwood
            Sh = 1.615*(Re*Sc*D/L)**(1/3)
            
            # Calculate Mass Transfer Coefficients
            KLiq = Sh*Diff/D
            KTot = 1/(1/KLiq + 1/(KMem*H))
            KTot[h2o_j] = KMem[h2o_j]
            
            QC = QC.reshape((Segs, numGas*2))
            C_lumen = QC[:,:numGas]
            C_shell = QC[:,numGas:]
            #!!! alternatively
            
            # # For the first segment:
            # for j in range(0, numGas):
            #     #if GasVec[j].Name == 'H2O':
            #     if j == h2o_j:
            #         dC[j+numGas] = (KTot[j]/(D/4))*(PVapH2O- (C[j+numGas]/sumCp_init)*P)*vFrac
            #         dC[j] = 0
            #     else:
            #         # dC[j] = (u/dx)*(Cs[j] - C[j]) - (KTot[j]/(D/4))*(C[j] - (C[j+numGas]/sumCp_init)*P/H[j])
            #         #!!! It seems like Cs should be the dissolved gas concentration in influent
            #         dC[j] = (u/dx)*(C_in[j] - C[j]) - (KTot[j]/(D/4))*(C[j] - (C[j+numGas]/sumCp_init)*P/H[j])
            #         dC[j+numGas] = (KTot[j]/(D/4))*(C[j]-(C[j+numGas]/sumCp_init)*P/H[j])*vFrac

            #     # Calculate the total gas concentration in the shell after the change
            #     sumCp_fin[0] += C[j+numGas] + dC[j+numGas]
            
            transmembrane = (KTot/(D/4))*(C_lumen - C_shell/sumCp_init*P/H)
            dC_lumen[0] = (u/dx)*(C_in - C_lumen[0]) - transmembrane[0]            
            dC_lumen[1:] = (u/dx)*(C_lumen[:-1] - C_lumen[1:]) - transmembrane[1:]
            dC_lumen[:,h2o_j] = 0
            dC_shell[:] = transmembrane*vFrac
            dC_shell[:,h2o_j] = (KTot[h2o_j]/(D/4))*(PVapH2O - (C_shell[:,h2o_j]/sumCp_init)*P)*vFrac
            
            # sumCp_fin = np.sum(C_shell+dC_shell, axis=1)
            

            # for i in range(1,Segs):
            #     # For the remaining segments:
            #     # Calculate the rate of change of the shell and lumen for all remaining segments.
            #     for j in range(0, numGas):
                    
            #         # Lumen
            #         dC[numVec*(i)+j] =  (u/dx)*(C[numVec*(i-1)+j] - C[numVec*(i)+j]) - (KTot[j]/(D/4))*(C[numVec*(i)+j] - (C[numVec*(i)+j+numGas]/sumCp_init)*(P/H[j]))

            #         # Shell
            #         dC[numVec*(i)+j+numGas] = (KTot[j]/(D/4))*(C[numVec*(i)+j] - (C[numVec*(i)+j+numGas]/sumCp_init)*(P/H[j]))*vFrac
                    
            #         # If the gas is H2O, then it follows a different formula:
            #         #if GasVec[j].Name == 'H2O':
            #         if j == self.h2o_j:
            #             dC[numVec*i+j+numGas] = (KTot[j]/(D/4))*(PVapH2O-(C[numVec*i+j+numGas]/sumCp_init)*P)*vFrac
            #             dC[numVec*i+j] = 0

            #         # Calculate the total gas concentration in the shell after the change
            #         sumCp_fin[i] += C[numVec*(i)+j+numGas] + dC[numVec*(i)+j+numGas]
                    
            # Re-scale the shell concentration so that vacuum pressure stays constant during the entire process. Given the change of concentration that we have calculated above for the shell, we can re-scale the shell concentrations with the sumCp at each segment. 
                
            # This FOR LOOP maintains consistent pressure in the shell
            # for i in range(0, Segs):
            #     for j in range(0, numGas):
            #         # Calculate the new concentration of gases in the shell
            #         newCp = (C[numVec*(i)+j+numGas] + dC[numVec*(i)+j+numGas])
                    
            #         # Re-scale the concentration to vacuum pressure
            #         newCp = (newCp/sumCp_fin[i])*P/(R*T)
                    
            #         # Calculate the actual difference of concentration that the function will output
            #         dC[numVec*(i)+j+numGas] = newCp- C[numVec*(i)+j+numGas]
            #         # Return the difference in concentration
                    
            #         #return dC
            
            new_C_shell = C_shell + dC_shell
            sumCp_fin = np.sum(new_C_shell, axis=1)
            dC_shell[:] = np.diag(sumCp_init/sumCp_fin) @ new_C_shell - C_shell
            dC[:] = np.hstack((dC_lumen, dC_shell)).flatten()
            _update_dstate()
        self._ODE = dy_dt


# For naming consistency
MembraneGasExtraction = GasExtractionMembrane