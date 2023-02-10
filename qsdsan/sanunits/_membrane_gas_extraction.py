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

__all__ = ('Membrane')

class GasExtractionMembrane(SanUnit):
    
    _N_ins = 1
    _N_outs = 1
    
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
        'CO2' : 3.5e-4, 
        'CH4': 1.3e-5,
        'O2': 1.2e-5, 
        'N2': 6e-6, 
        'H2O': 1
        }
    
    _HenrySlope = {
        'H2': 640, 
        'CO2' : 2600, 
        'CH4' : 1900,
        'O2': 1800, 
        'N2': 1300, 
        'H2O': 1
        }
    
    _WilkeChang = {
        'H2': 9.84, 
        'CO2': 2.6, 
        'CH4': 2.2, 
        'O2': 1.90,
        'N2': 1.77, 
        'H2O': 1
        }
    
    # Constructor: Initialize the instance variables
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None,   FiberID=190e-6, 
                  FiberOD=300e-6, NumTubes=1512, ShellDia=1.89e-2, SurfArea=0.1199,   
                  GasID = ['H2', 'CO2', 'CH4', 'O2', 'N2', 'H2O'], PVac = 97.325, 
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
        self.segs = segs              # Number of segments ??Ask Ian??
        #self.Volume = VolBatchTank   # Volume of the bioreactor (Don't think this is needed)
        
        dct_gas_perm = GasPerm or self._GasPerm
        self.set_GasPerm(**dct_gas_perm)
        dct_gas_hpf =  HenryPreFac or self._HenryPreFac
        self.set_HenryPreFac(**dct_gas_hpf)
        dct_gas_hs = HenrySlope or self._HenrySlope
        self.set_HenrySlope(**dct_gas_hs)
        dct_gas_wc = WilkeChang or self._WilkeChang
        self.set_WilkeChang(**dct_gas_wc)
        
        inf, = self.ins
        cmps = inf.components
        self.indexer = cmps.index
        # self.idx ensures that the indexing in further code is only for gases 
        # and not all components in the influent
        self.idx = cmps.indices(self.GasID) 
        # to save the index of water in the array of gases, to be used later 
        for i, ID in enumerate(GasID):
            if ID == 'H2O': 
                self.h2o_j = i
                break
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
            raise ValueError('Outer diameter of fiber expected from user')
    
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
        TRefH = 298.15              # Reference T for Henry's Law [K]
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
    
    def _init_state(self):
        inf, = self.ins
        cmps = inf.components
        # ASSUMPTION: Only 1 influent 
        C = self._ins_QC[:-1]/cmps.chem_MW*cmps.i_mass # conc. in mol/m^3 as defined by Ian 
        Cs = C[self.idx]     #self.idx ensures its only for gases 
        #Q = self._ins_QC[-1]
        Seg = self.segments
        numGas = len(self.GasID)
        
        self._state = np.zeros(2*Seg*numGas)
        for i in range(0, 2*Seg*numGas, 2*numGas):
            for j in range(numGas):
                self._state[j+i] = Cs[j] 
                
        self._dstate = self._state*0
        
            
    #def transientGasExtraction(t, C, ExpCond, GasVec, Mem, Segs, SS):
    def _compile_ODE(self):
        # Synthesizes the ODEs to simulate a batch reactor with side-flow gas extraction. The code takes in an object of class Membrane (Mem) and an array of objects of class Gas (GasVec). It also takes in an array of experimental conditions ExpCond. 
        
        # Extract Operating Parameters from ExpCond
        Q = self.ins.F_vol*(1000/60)  # Volumetric Flowrate [L/min]
        T = self.ins[0].T  # Temperature [K]
        P = self.PVac*1000 # Vacuum Pressure [Pa]
        #V = self.Volume  # Volume of the Batch Tank [L]

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
        Segs = self.segments        # Number of segments? Ask Ian
        vFrac = self.VolFrac        # Lumen/Shell Volume Fraction [m^3/m^3]

        # Pre-allocate vectors for gas thermophysical properties
        numGas = len(self.GasID)
        numVec = 2*numGas
# =============================================================================
#         Diff = np.zeros(numGas)
#         Perm_SI = np.zeros(numGas)
#         H = np.zeros(numGas)
#         Cin = np.zeros(numGas)
#         MM = np.zeros(numGas)
# =============================================================================

        inf, = self.ins
        cmps = inf.components
        # Extract Gas Properties
        #for i in range(0,numGas):
        
        #Diff[i] = GasVec[i].Diff()
        Diff = self.Diff
            
        #Perm_SI[i] = Mem.PermDict[GasVec[i].Name]
        Perm_SI = self._gasp 
            
        #H[i] = GasVec[i].HeL()
        H = self.HeL

        #Diff = np.array([0.0265e-7, 0.0296e-7, 0.0253e-7, 0.3199e-7])
        # Calculate dx and u
        dx = L/Segs # Length of Segments [m]
        u = Q/((np.pi*D**2/4)*num_tubes*1000*60)    # Linear Flow Velocity [m/s]

        # Calculate the Kinematic Viscosity of Water
        Tb = T/300  # Reduced Temperature []
        mu = (1e-6)*(280.68*(Tb**(-1.9)) + 511.45*(Tb**(-7.7)) + 61.131*(Tb**(-19.6)) + 0.45903*(Tb**(-40))) # Absolute Viscosity of Water [Pa s]
        rho = (-13.851 + 0.64038*T - 1.9124e-3*T**2 + 1.8211e-6*T**3)*18   # Liquid Density of Water [kg/m^3] from DIPPR
        nu = mu/rho # Kinematic Viscosity of Water [m^2/s]

        # Calculate the dimensionless numbers
        # Reynolds
        Re = u*D/nu
        # Schmidt
        Sc = nu/Diff
        # Sherwood
        Sh = 1.615*(Re*Sc*D/L)**(1/3)
        
        
        # Calculate Mass Transfer Coefficients
        KMem = Perm_SI/l
        KLiq = Sh*Diff/D
        KTot = 1/(1/KLiq + 1/(KMem*H))
        for j in range(0, numGas):
            #if GasVec[j].Name == 'H2O':
            if j == self.h2o_j:
                KTot[j] = KMem[j]
      
        # Initialize
        C = self._state
        dC = self._dstate

        sumCp_init = P/(R*T)
        sumCp_fin = np.zeros(Segs)
        
        C = self._ins_QC[:-1]/cmps.chem_MW*cmps.i_mass # conc. in mol/m^3 as defined by Ian 
        Cs = C[self.idx] #self.idx ensures its only for gases 
        
        # For the first segment:
        for j in range(0, numGas):
            #if GasVec[j].Name == 'H2O':
            if j == self.h2o_j:
                dC[j+numGas] = (KTot[j]/(D/4))*(PVapH2O- (C[j+numGas]/sumCp_init)*P)*vFrac
                dC[j] = 0
            else:
                dC[j] = (u/dx)*(Cs[j] - C[j]) - (KTot[j]/(D/4))*(C[j] - (C[j+numGas]/sumCp_init)*P/H[j])
                dC[j+numGas] = (KTot[j]/(D/4))*(C[j]-(C[j+numGas]/sumCp_init)*P/H[j])*vFrac

            # Calculate the total gas concentration in the shell after the change
            sumCp_fin[0] += C[j+numGas] + dC[j+numGas]

        for i in range(1,Segs):
            # For the remaining segments:
            # Calculate the rate of change of the shell and lumen for all remaining segments.
            for j in range(0, numGas):
                
                # Lumen
                dC[numVec*(i)+j] =  (u/dx)*(C[numVec*(i-1)+j] - C[numVec*(i)+j]) - (KTot[j]/(D/4))*(C[numVec*(i)+j] - (C[numVec*(i)+j+numGas]/sumCp_init)*(P/H[j]))

                # Shell
                dC[numVec*(i)+j+numGas] = (KTot[j]/(D/4))*(C[numVec*(i)+j] - (C[numVec*(i)+j+numGas]/sumCp_init)*(P/H[j]))*vFrac
                
                # If the gas is H2O, then it follows a different formula:
                #if GasVec[j].Name == 'H2O':
                if j == self.h2o_j:
                    dC[numVec*i+j+numGas] = (KTot[j]/(D/4))*(PVapH2O-(C[numVec*i+j+numGas]/sumCp_init)*P)*vFrac
                    dC[numVec*i+j] = 0

                # Calculate the total gas concentration in the shell after the change
                sumCp_fin[i] += C[numVec*(i)+j+numGas] + dC[numVec*(i)+j+numGas]

        # Re-scale the shell concentration so that vacuum pressure stays constant during the entire process. Given the change of concentration that we have calculated above for the shell, we can re-scale the shell concentrations with the sumCp at each segment. 
            
        for i in range(0, Segs):
            for j in range(0, numGas):
                # Calculate the new concentration of gases in the shell
                newCp = (C[numVec*(i)+j+numGas] + dC[numVec*(i)+j+numGas])
                
                # Re-scale the concentration to vacuum pressure
                newCp = (newCp/sumCp_fin[i])*P/(R*T)
                
                # Calculate the actual difference of concentration that the function will output
                dC[numVec*(i)+j+numGas] = newCp- C[numVec*(i)+j+numGas]

        # Return the difference in concentration
        return dC