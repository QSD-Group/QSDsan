#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from math import log, exp
from thermosteam import indexer, equilibrium
from qsdsan import SanUnit, Construction

__all__ = ('MembraneDistillation',)


# =============================================================================
# Membrane Distillation
# =============================================================================

class MembraneDistillation(SanUnit):
    '''
    Membrane distillation recovers nitrogen as ammonia sulfate based on vapor
    pressure difference across the hydrophobic membrane. Ignore water flux across
    membrane since it will not affect system performance (either TEA or LCA).
    
    Parameters
    ----------
    ins : Iterable(stream)
        influent, acid, base, mem_in.
    outs : Iterable(stream)
        ammonium sulfate, ww, mem_out.
    influent_pH: float
        Influent pH.
    target_pH: float
        Target pH for membrane distillation.
    N_S_ratio: float
        mol(N) to mol(S) ratio.
    m2_2_m3: float
        m2 to m3 factor, 1/specific surface area, [m3/m2].
    Dm: float
        NH3 molecular diffusity in air, [m2/s]. 
    porosity: float
        Membrane porosity.
    thickness: float
        Membrane thickness, [m].
    tortuosity: float
        Membrane tortuosity.
    Henry: float
        NH3 Henry constant, [atm*m3/mol].
    Ka: float
        Overall mass transfer coefficient, [m/s].
    capacity: float
        Membrane treatement capacity (permeate flux), [kg/m2/h].
    membrane_price: float
        Membrane price, [$/kg] ([$/m2]).
        
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    .. [2] Doran, P. M. Chapter 11 - Unit Operations. In Bioprocess Engineering
        Principles (Second Edition); Doran, P. M., Ed.; Academic Press: London,
        2013; pp 445–595. https://doi.org/10.1016/B978-0-12-220851-5.00011-3.
    .. [3] Spiller, L. L. Determination of Ammonia/Air Diffusion Coefficient Using
        Nafion Lined Tube. Analytical Letters 1989, 22 (11–12), 2561–2573.
        https://doi.org/10.1080/00032718908052375.
    .. [4] Scheepers, D. M.; Tahir, A. J.; Brunner, C.; Guillen-Burrieza, E.
        Vacuum Membrane Distillation Multi-Component Numerical Model for Ammonia
        Recovery from Liquid Streams. Journal of Membrane Science
        2020, 614, 118399. https://doi.org/10.1016/j.memsci.2020.118399.
    .. [5] Ding, Z.; Liu, L.; Li, Z.; Ma, R.; Yang, Z. Experimental Study of Ammonia
        Removal from Water by Membrane Distillation (MD): The Comparison of Three
        Configurations. Journal of Membrane Science 2006, 286 (1), 93–103.
        https://doi.org/10.1016/j.memsci.2006.09.015.
    .. [6] Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
        Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
        Thermal Efficiency, Sensitivity Study and Cost Estimation.
        Journal of Membrane Science 2008, 323 (1), 85–98.
        https://doi.org/10.1016/j.memsci.2008.06.006.
    .. [7] Kogler, A.; Farmer, M.; Simon, J. A.; Tilmans, S.; Wells, G. F.;
        Tarpeh, W. A. Systematic Evaluation of Emerging Wastewater Nutrient Removal
        and Recovery Technologies to Inform Practice and Advance Resource
        Efficiency. ACS EST Eng. 2021, 1 (4), 662–684.
        https://doi.org/10.1021/acsestengg.0c00253.
    .. [8] Pikaar, I.; Guest, J.; Ganigue, R.; Jensen, P.; Rabaey, K.; Seviour, T.;
        Trimmer, J.; van der Kolk, O.; Vaneeckhaute, C.; Verstraete, W.; Resource
        Recovery from Water: Principles and Applicaiton. IWA 2022.
    '''
    _N_ins = 4
    _N_outs = 4
    _F_BM_default = {'Membrane': 1}
    
    _units = {'Area': 'm2',
              'Total volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 influent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 [1]
                 target_pH=10,
                 N_S_ratio=2,
                 # S is excess since not all N can be transfered to form ammonia sulfate
                 # for now, assume N_S_ratio = 2 is ok
                 m2_2_m3=1/1200, # specific surface area, for hollow fiber membrane [2]
                 Dm=2.28*10**(-5), # (2.28 ± 0.12)*10^-5 m^2/s NH3 molecular diffusity in air [3]
                 # (underestimate, this value may be at 15 or 25 C, our feed is 60 C, should be higher)
                 porosity=0.9, # [4]
                 thickness=7*10**(-5), # m [4]
                 tortuosity=1.2, # [4]
                 Henry=1.61*10**(-5), # atm*m3/mol
                 # https://webwiser.nlm.nih.gov/substance?substanceId=315&identifier=\
                 # Ammonia&identifierType=name&menuItemId=81&catId=120#:~:text=The%20\
                 # Henry's%20Law%20constant%20for,m%2Fmole(2). (accessed 11-11-2022)
                 Ka=1.75*10**(-5), # overall mass transfer coefficient 1.2~2.3*10^-5 m/s [5]
                 capacity=6.01, # kg/m2/h [6]
                 # permeate flux, similar values can be found in many other papers ([176], [222] ,[223] in [7])
                 # [177], [223] in [7] show high nitrogen recovery ratio (>85% under optimal conditions)
                 membrane_price=93.29,
                 yearly_operation_hour=7920,
                 ): # $90/m2 2008 [6]
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.influent_pH = influent_pH
        self.target_pH = target_pH
        self.N_S_ratio = N_S_ratio
        self.m2_2_m3 = m2_2_m3
        self.Dm = Dm
        self.porosity = porosity
        self.thickness = thickness
        self.tortuosity = tortuosity
        self.Henry = Henry
        self.Ka = Ka
        self.capacity = capacity
        self.membrane_price = membrane_price
        self.construction = [
            Construction('membrane', linked_unit=self, item='RO', quantity_unit='m2'),
            ]
    
    def _run(self):
        
        influent, acid, base, mem_in = self.ins
        ammoniumsulfate, ww, mem_out, ammoniumsulfatesolution = self.outs
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        if self.CHG.CHGout_N == 0:
            ww.copy_like(influent)
            self.membrane_area = self.F_vol_in*1000/self.capacity
        else:
            NaOH_conc = 10**(self.target_pH - 14) - 10**(self.influent_pH - 14)
            NaOH_mol = NaOH_conc*self.ins[0].F_mass
            base.imass['NaOH'] = NaOH_mol*39.997/1000
            
            acid.imass['H2SO4'] = self.CHG.CHGout_N/14.0067/self.N_S_ratio*98.079
            acid.imass['H2O'] = acid.imass['H2SO4']*1000/98.079/0.5*1.05 -\
                                acid.imass['H2SO4']
            
            self.pKa = pKa = -log(exp(52.22/8.3145*1000*(1/298 - 1/self.ins[0].T))*10**(-9.252), 10)
            # calculation of pKa under different T follows van't Hoff relationship [8] page 292

            ammonia_to_ammonium = 10**(-pKa)/10**(-self.target_pH)
            ammonia = self.CHG.CHGout_N*ammonia_to_ammonium/(1 +\
                       ammonia_to_ammonium)*17.031/14.0067
            others = influent.F_mass - ammonia
            
            self.membrane_area = self.F_vol_in*1000/self.capacity
            N2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.79*1.204
            O2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.21*1.204
            # N2:O2 = 0.79:0.21 in the air, air density is 1.204 kg/m3
            # https://en.wikipedia.org/wiki/Density_of_air#:~:text=It%20also%20\
            # changes%20with%20variation,International%20Standard%20Atmosphere%2\
            # 0(ISA). (accessed 11-14-2022)
            
            imass = indexer.MassFlowIndexer(l=[('H2O', others),
                                               ('NH3', ammonia),
                                               ('N2', 0),
                                               ('O2', 0)],
                                            g=[('H2O', 0),
                                               ('NH3', 0), 
                                               ('N2', N2_in_air),
                                               ('O2', O2_in_air)])
            # N2 amount will be changed based on design, maybe also add O2 (N2:O2 = 4:1)
            
            vle = equilibrium.VLE(imass)
            vle(T=influent.T, P=influent.P)
            X_NH3_f_m = vle.imol['g','NH3']/(vle.imol['g','H2O'] + vle.imol['g','NH3'])
            X_NH3_f = vle.imol['l','NH3']/(vle.imol['l','H2O'] + vle.imol['l','NH3'])
    
            km = self.Dm*self.porosity/self.tortuosity/self.thickness
            
            dimensionless_Henry = self.Henry/8.20575/(10**(-5))/influent.T # H' = H/RT
            # https://www.sciencedirect.com/topics/chemistry/henrys-law#:~:text=\
            # Values%20for%20Henry's%20law%20constants,gas%20constant%20(8.20575%\
            # 20%C3%97%2010 (accessed 11-11-2022)
            
            kf = 1/(1/self.Ka - 1/dimensionless_Henry/km*(1 + 10**(-pKa)*10**\
                 (-self.target_pH)/(10**(-14))))
            
            J = kf*ammonia/influent.F_mass*1000*log(X_NH3_f_m/X_NH3_f)*3600 # in kg/m2/h
            
            NH3_mass_flow = J*self.membrane_area
            
            ammonia_transfer_ratio = min(1, NH3_mass_flow/ammonia)
            
            ammoniumsulfate.imass['NH42SO4'] = ammonia*ammonia_transfer_ratio/34.062*132.14
            ammoniumsulfatesolution.imass['H2O'] = acid.imass['H2O']
            ammoniumsulfatesolution.imass['H2SO4'] = acid.imass['H2SO4'] -\
                                             ammoniumsulfate.imass['NH42SO4']/\
                                             132.14*98.079
                                            
            ww.copy_like(influent) # ww has the same T and P as influent
            
            ww.imass['N'] = self.CHG.CHGout_N*(1 - ammonia_to_ammonium/(1 +\
                             ammonia_to_ammonium)*ammonia_transfer_ratio)
                                                              
            ww.imass['C'] = self.CHG.CHGout_C
                             
            ww.imass['P'] = self.CHG.CHGout_P
                       
            ww.imass['H2O'] -= (ww.imass['C'] + ww.imass['N'] + ww.imass['P'])
            
            ww.imass['H2O'] += self.ins[2].F_mass
            
            ammoniumsulfate.T = ammoniumsulfatesolution.T = acid.T
            ammoniumsulfate.P = ammoniumsulfatesolution.P = acid.P
            # ammoniumsulfate has the same T and P as acid
            
            #!!! Use the lifetime attr
            mem_in.imass['Membrane'] = 0.15*self.membrane_area/7920 # kg/hr (m2/hr)
            mem_out.copy_like(mem_in)
            # add membrane as streams to include 15% membrane replacement per year [6]
    
    @property
    def N_recovery_ratio(self):
        return 1 - self.outs[1].imass['N']/self.CHG.CHGout_N
        
    def _design(self):
        Design = self.design_results
        Design['Area'] = self.membrane_area
        Design['Total volume'] = Design['Area']*self.m2_2_m3
        
        self.construction[0].quantity = Design['Area']
    
    def _cost(self):
        self.ins[3].price = self.membrane_price
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs
        purchase_costs['Membrane'] = Design['Area']*self.membrane_price