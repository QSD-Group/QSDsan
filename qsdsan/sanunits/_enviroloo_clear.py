#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created by Yuyao Huang and Siqi Tang for Enviroloo (EL) Clear Toilet system
'''
# %%
from qsdsan import WasteStream
from qsdsan import SanUnit, Construction
from qsdsan.equipments import Blower
from qsdsan.sanunits import IdealClarifier
from qsdsan.sanunits._tank import StorageTank
from qsdsan.processes._decay import Decay
from qsdsan.utils import ospath, load_data, data_path, price_ratio
# %% This callable file will be reposited to qsdsan.SanUnit subbranch with the name of _enviroloo
__all__ = (
    'EL_CT', # Collection tank
    'EL_PC', # Primary clarifier
    'EL_Anoxic', # Anoxic tank
    'EL_Aerobic', # Aerobic tank
    'EL_MBR', # Membrane filter
    'EL_CWT', # Clear water tank
    'EL_PT', # Pressure tank
    'EL_blower', # blower
    'EL_System', # System-level summary
    'EL_Housing', # Housing of EL_System, such as equipment's armor
    )

EL_su_data_path = ospath.join(data_path, 'sanunit_data/el')  # need change

# %%
CollectionTank_path = ospath.join(EL_su_data_path, '_EL_CT.tsv')
@price_ratio()
class EL_CT(StorageTank):

    '''
    Name
    ----
    Collection tank in the Enviroloo (EL) Clear Toilet system.
    
    Parameters
    ----------
    Ins: 
    (1) Mixed wastewater
    (2) Primary clarifier effluent spill
    (3) Clear water tank effluent spill 
    (4) Primary clarifier sludge return

    Outs:
    (1) Treated water
    (2) Methane (CH4)
    (3) Nitrous oxide (N2O)

    Attributes
    ----------
    length_to_diameter : float
        Ratio of the tank length to diameter.
    vessel_material : str
        Material used for constructing the vessel.
    sludge_moisture_content : float
        Moisture content of the sludge (mass of water/total mass).
    COD_removal : float
        Fraction of COD removed in the collection tank.

    References
    ----------
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

    See Also
    ----------
    class:`biosteam.units.StorageTank`
    '''
    _N_ins = 4 
    _N_outs = 1
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, ppl=None, baseline_ppl=None,
                 vessel_type= 'Horizontal', tau=24, V_wf=None, vessel_material='Stainless steel', kW_per_m3=0.1,
                 init_with='WasteStream', F_BM_default=1,
                 include_construction=True, **kwargs):
        StorageTank.__init__(self, 
                      # Basic parameters
                      ID=ID, # The unique identifier of the tank
                      ins=ins, # The input stream to the tank
                      outs=outs, # The output streams from the tank
                      thermo=thermo, # The thermodynamic property package for simulating physical or chemical processes
                      # Other control parameters
                      init_with=init_with, # The method to initialize the tank contents.
                      F_BM_default=F_BM_default, # The default bare module factor for the tank cost estimation
                      include_construction=include_construction, # A Boolean value indicating whether the tank's construction material is considered in the cost analysis or life cycle assessment.
                      # Design parameters
                      vessel_type=vessel_type, # The type of the tank
                      tau=tau, # The retention time of the tank contents, the important parameters involved in mixing, reaction or separation processes.
                      V_wf=V_wf, # The volume working fraction of the tank, the ratio of the volume of the tank contents to the total volume of the tank.
                      vessel_material=vessel_material, # The material of the tank
                      kW_per_m3=kW_per_m3, # The power consumption per unit volume of the tank
                      )
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
    
        self._mixed = WasteStream()
        
        data = load_data(path=CollectionTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item = 'Steel', linked_unit=self, quantity_unit='kg'),]
      
        
    def _run(self):

        # Input stream
        WasteWater = self.ins[0]
        sludge_return = self.ins[1]
        PC_spill_return = self.ins[2]
        MT_spill_return = self.ins[3]

        # Output stream
        TreatedWater = self.outs[0]

        # Ensure existing input streams
        input_streams = [WasteWater, sludge_return]
        if PC_spill_return:  # If there is PC_spill_return
            input_streams.append(PC_spill_return)
        if MT_spill_return:  # If there is MT_spill_return
            input_streams.append(MT_spill_return)

        # Mix all inputs into a single stream
        self._mixed.mix_from([input_streams])

        # Copy the mixed result to the outflow
        TreatedWater.copy_like(self._mixed)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.collection_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
    
        price_ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * price_ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        self.power_utility(self.power_demand_CT)
    
    def _calc_replacement_cost(self):
        scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CT_replacement_cost = (
            self.collection_tank_cost / self.collection_tank_lifetime +               
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        CT_replacement_cost = CT_replacement_cost / (365 * 24); # convert to USD/hr
        return CT_replacement_cost
        



# %%
PrimaryClarifier_path = ospath.join(EL_su_data_path, '_EL_PC.tsv'); # need change

@price_ratio()
class EL_PC(IdealClarifier):
    
    """
    Name
    ----
    Primary clarifier in the Enviroloo (EL) Clear Toilet system.

    Introduction
    ------------
    The primary treatment of the EL system uses anaerobic digestion
    to treat wastes (similar to a septic tank).

    It can be used in conjunction with a membrane bioreactor (MBR)
    to recovery N and P as struvite.

    The following impact items should be pre-constructed for life cycle assessment:
    FRP.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from lift pump
    (2) nitrate return flow from membrane tank

    Outs:
    (1) effluent of treated wastewater
    (2) sludge return flow to collection tank
    (3) spill flow to collection tank
    (4) fugitive CH4 emission
    (5) fugitive N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.PrimaryClarifier module

     """
    
    _N_ins = 2
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, max_overflow=None,
                 ppl = None, baseline_ppl = None, sludge_flow_rate=None,
                 solids_removal_efficiency=None, 
                 F_BM_default=1, init_with='WasteStream', **kwargs):
        """
        Initialize the primary clarifier with default parameters:
        - sludge_flow_rate: Default to  m3/d
        - solids_removal_efficiency: Default to 50% (0.5)
        """
        IdealClarifier().__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, sludge_flow_rate=sludge_flow_rate,
                                  solids_removal_efficiency=solids_removal_efficiency,
                                  isdynamic=isdynamic, init_with=init_with, F_BM_default=F_BM_default
                                  )
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        self._mixed = WasteStream()  # Create a temporary mixed stream and its properties and actions are like 'WasteStream'
        self._f_spill = None  # Spill return
        self._f_overflow = None  # Overflow
        self.max_overflow = max_overflow
        # self.if_with_MBR = if_with_MBR

        data = load_data(path=PrimaryClarifier_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item = 'Steel', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
        """
        Run the enhanced primary clarifier process with spill return handling.
        """
        # Input streams
        WasteWater = self.ins[0]  # Incoming wastewater
        MT_sludge_return = self.ins[1]  # Returned sludge from membrane tank
        
        # Outputs
        TreatedWater = self.outs[0]  # Normal overflow to anoxic tank
        spill_return = self.outs[1]  # Spill return to the collection tank
        PC_sludg_return = self.outs[2]  # Settled sludge to upstream
        
        # Step 1: Mix influent and returned sludge
        self._mixed.mix_from([WasteWater, MT_sludge_return])
        
        # Step 2: Calculate normal overflow and underflow
        TreatedWater.copy_like(self._mixed)
        PC_sludg_return.copy_like(self._mixed)
        
        # Apply solids removal efficiency
        PC_sludg_return.F_mass[:] *= self.solids_removal_efficiency
        TreatedWater.F_mass[:] *= (1 - self.solids_removal_efficiency)
        
        # Step 3: Handle spill return based on overflow volume
        if self.max_overflow is not None:
            if TreatedWater.F_vol > self.max_overflow:
                # Calculate excess volume
                spill_vol = TreatedWater.F_vol - self.max_overflow
                
                # Calculate the fraction of spill
                self._f_spill = spill_vol / TreatedWater.F_vol
                
                # Update the normal overflow fraction
                self._f_overflow = 1 - self._f_spill
                
                # Assign spill fraction to spill return
                spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill
                
                # Adjust the normal overflow to within max capacity
                TreatedWater.F_mass[:] *= self._f_overflow
            else:
                # No spill return if overflow is within capacity
                spill_return.empty()
                self._f_spill = 0.0
                self._f_overflow = 1.0
        else:
            # No max_overflow defined, no spill handling
            spill_return.empty()
            self._f_spill = 0.0
            self._f_overflow = 1.0

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.PC_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
    
        ratio = self.price_ratio 
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        self.power_utility(self.power_demand_PC)
    
    def _calc_replacement_cost(self):
        scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        PC_replacement_cost = (
            self.PC_tank_cost / self.PC_tank_lifetime +               
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        PC_replacement_cost = PC_replacement_cost / (365 * 24); # convert to USD/hr
        return PC_replacement_cost

# %%
Anoxic_path = ospath.join(EL_su_data_path, '_EL_Anoxic.tsv')

@price_ratio()
class EL_Anoxic(SanUnit, Decay):
    '''
    Anoxic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from primary clarifier
    (2) nitrate return flow from membrane tank
    (3) glucose addition
    (4) agitation pump

    Outs:
    (1) effluent of treated wastewater
    (2) fugitive CH4 emission
    (3) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.PrimaryClarifier module
    
    
    '''
    _N_ins = 4
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    baseline_ppl = 30

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',),  F_BM_default=1, ppl = None, baseline_ppl = None,
                 if_capture_biogas=False, if_N2O_emission=True, **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path=Anoxic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),]      
        
    def _run(self):
        # Input stream
        WasteWater = self.ins[0]
        sludge_return = self.ins[1]
        glucose = self.ins[2]
        agitation = self.ins[3]
        
        # Output stream
        TreatedWater = self.outs[0]
        CH4_emission = self.outs[1]
        N2O_emission = self.outs[2]

        # Mix all inputs into a single stream
        agitation.F_mass = 0  # Attend nothing
        WasteWater.F_mass += sludge_return.F_mass + agitation.F_mass
        for component in WasteWater.components:
            WasteWater.imass[component] += sludge_return.imass[component]

        # Add glucose to the waste water
        WasteWater.imass['Glucose'] += glucose.imass['Glucose']
        WasteWater.F_mass += glucose.F_mass

        # Simulate complete glucose consumption
        consumed_glucose = WasteWater.imass['Glucose']
        WasteWater.imass['Glucose'] -= consumed_glucose
        WasteWater.F_mass -= consumed_glucose

        # Adjust N2O emission factor based on reduction ratio
        original_N2O_EF = self.N2O_EF_decay 
        self.N2O_EF_decay *= 0.25

        # Call Decay._first_order_run with the updated N2O_EF_decay
        super()._first_order_run(waste=WasteWater,
                                 treated=TreatedWater,
                                 CH4=CH4_emission,
                                 N2O=N2O_emission
                                 )

        # Restore the original N2O emission factor
        self.N2O_EF_decay = original_N2O_EF

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scale
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.anoxic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Chemcial_glucose'] = self.chemical_glucose_dosage * self.ins[0] * self.chemical_glucose_price; # make sense the unit of treated water flow

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()

        self.power_utility(self.power_demand_PC)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Anoxic_tank_replacement_cost = (self.anoxic_tank_cost /self.anoxic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Anoxic_tank_replacement_cost = Anoxic_tank_replacement_cost / (365 * 24); # convert to USD/hr
        return Anoxic_tank_replacement_cost

# %%
Aerobic_path = ospath.join(EL_su_data_path, '_EL_Aerobic.tsv'); # need change

@price_ratio()
class EL_Aerobic(SanUnit, Decay):

    '''
    Aerobic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from anoxic tank
    (2) PAC addition
    (3) blower

    Outs:
    (1) effluent of treated wastewater
    (2) fugitive CH4 emission
    (3) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.PrimaryClarifier module
    '''

    _N_ins = 3; # treated water, PAC, blower
    _N_outs = 3; # treated water, CH4, N2O
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
                 if_capture_biogas=False, if_N2O_emission=True,**kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path=Aerobic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),]  
    
    def _run(self):

        # Input streams
        WasteWater = self.ins[0]
        PAC = self.ins[1]
        air = self.ins[2]

        # Output streams
        TreatedWater = self.outs[0]
        CH4_emission = self.outs[1]
        N2O_emission = self.outs[2]

        # Mix PAC into WasteWater (affects mass balance)
        WasteWater.imass['PAC'] += PAC.F_mass  
        WasteWater.F_mass += PAC.F_mass

        # Air (used for aeration) does not affect mass balance
        # Air is only used for operational purposes, no changes in WasteWater mass

        # Copy updated WasteWater to TreatedWater
        TreatedWater.copy_like(WasteWater)

        # No biogas captured in aerobic tanks, calculate emissions
        super()._first_order_run(waste=WasteWater,
                                 treated=TreatedWater,
                                 CH4=CH4_emission,
                                 N2O=N2O_emission
                                 )

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.arobic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Chemical_PAC'] = self.chemical_PAC_dosage * self.ins[0] * self.chemical_PAC_price

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()

        self.power_utility(self.power_demand_AerobicTank) # kWh, defined in .tsv file
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
        Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24); # convert to USD/hr
        return Aerobic_tank_replacement_cost

# %%

MBR_path = ospath.join(EL_su_data_path, '_EL_MBR.tsv')

@price_ratio()
class EL_MBR(SanUnit, Decay):

    '''
    Aerobic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from aerobic tank
    (2) blower

    Outs:
    (1) effluent of treated wastewater
    (2) nitrate return flow to anoxic tank
    (3) nitrate return flow to primary clarifier
    (4) fugitive CH4 emission
    (5) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the exposan.eco-san.MBR module
    '''
    _N_ins = 2  # treated water, Blower
    _N_outs = 5  # treated water, CH4, N2O, Nitrate return to Primary Clarifier, Nitrate return to Anoxic Tank
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.6

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
                 if_capture_biogas=False, if_N2O_emission=True, **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path=MBR_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
                             Construction(item ='Membrane_Material', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
        
        # Input streams
        WasteWater = self.ins[0]
        air = self.ins[1]

        # Output streams
        TreatedWater = self.outs[0]
        PC_sludge_return = self.outs[1]
        AnoT_sludge_return = self.outs[2]
        CH4_emission = self.outs[3]
        N2O_emission = self.outs[4]

        # Air does not affect mass balance
        # Air is used for aeration, so we do not modify WasteWater mass

        # Step 1: Copy WasteWater to TreatedWater
        TreatedWater.copy_like(WasteWater)

        # Step 2: Split sludge for returns (1% to PC, 99% to AnoT)
        total_f_sludge = 0.05  # Total sludge accounts for 5% of the wastewater mass
        f_PC_sludge = 0.01     # Primary clarifier accounts for 1% of the total sludge
        f_AnoT_sludge = 0.99   # Anoxic tank accounts for 99% of the total sludge

        # Sludge mass
        total_sludge_mass = WasteWater.F_mass * total_f_sludge
        PC_sludge_return.mass = total_sludge_mass * f_PC_sludge
        AnoT_sludge_return.mass = total_sludge_mass * f_AnoT_sludge

        # Update TreatedWater after sludge removal
        TreatedWater.F_mass -= total_sludge_mass
        for component in TreatedWater.components:
            TreatedWater.imass[component] -= (
                PC_sludge_return.imass[component] + AnoT_sludge_return.imass[component]
            )

        # Step 3: Calculate emissions using Decay._first_order_run
        super()._first_order_run(waste=WasteWater,
                                 treated=TreatedWater,
                                 biogas=None,  # Assume no biogas is captured
                                 CH4=CH4_emission,
                                 N2O=N2O_emission
                                 )

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['Membrane_Material'] = constr[1].quantity = self.membrane_material_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['MBR_tank'] = self.MBR_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
        C['Membrane_cleaning'] = self.membrane_cleaning_fee

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        self.power_utility(self.power_demand_MBR)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        MBR_replacement_cost = (
            self.MBR_tank_cost / self.MBR_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime
            ) * scale
        MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
        return MBR_replacement_cost

# %%

ClearWaterTank_path = ospath.join(EL_su_data_path, '_EL_CWT.tsv')

@price_ratio()
class EL_CWT(StorageTank):

    '''
    Introduction
    ------------
    To only collect the treated water

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from self-priming pump
    (2) O3 dosing
    (3) air-dissolve influent

    Outs:
    (1) effluent of treated wastewater (clear water)
    (2) spill flow to collection tank
    (3) influent to air-dssolved pump
    (4) fugitive CH4 emission
    (5) fugitive N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    _N_ins = 3; # number of ins
    _N_outs = 5; # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    def __init__(self, ID = '', ins = None, outs = (), vessel_material=None, V_wf = None, include_construction = True, 
                 length_to_diameter = None, F_BM_default = 1, kw_per_m3 = None, vessel_type=None,
                 tau = None, max_overflow=None, ppl = None, baseline_ppl = None,
                 thermo = None, init_with = 'WasteStream', **kwargs):
        StorageTank.__init__(self, ID, ins, outs, thermo = thermo, init_with = init_with, include_construction= include_construction,
                             kw_per_m3= kw_per_m3, vessel_type= vessel_type, tau= tau,
                             vessel_material= vessel_material, V_wf= V_wf, F_BM_default= F_BM_default)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.length_to_diameter = length_to_diameter
        self.max_overflow = max_overflow

        data = load_data(path = ClearWaterTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
        
        # Input streams
        WasteWater = self.ins[0]
        ozone = self.ins[1]
        air = self.ins[2]

        # Output streams
        TreatedWater = self.outs[0]
        CT_spill_return = self.outs[1]
        PC_spill_return = self.outs[2]
 
        # Constants for overflow conditions
        CT_spill_fraction = 0.6
        PC_spill_fraction = 0.4

        # Step 1: Process inputs
        # Air and ozone are balanced, not affecting mass balance
        TreatedWater.copy_like(WasteWater)

        # Step 2: Handle overflow if necessary
        if TreatedWater.F_vol > self.max_overflow:
            # Calculate excess volume
            spill_vol = TreatedWater.F_vol - self.max_overflow

            # Calculate spill fractions
            self._f_spill = spill_vol / TreatedWater.F_vol
            self._f_treated = 1 - self._f_spill

            # Assign spill fractions to CT_spill_return and PC_spill_return
            CT_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * CT_spill_fraction
            PC_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * PC_spill_fraction

            # Adjust the normal treated water to within max capacity
            TreatedWater.F_mass[:] *= self._f_treated
        else:
            # No spill return if overflow is within capacity
            CT_spill_return.empty()
            PC_spill_return.empty()
            self._f_spill = 0.0
            self._f_treated = 1.0

    def _design(self):
        design = self.design_results;
        constr = self.construction;
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl); # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Clear water tank'] = self.clear_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['O3 generator'] = self.O3_generation_machine_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        self.power_utility(self.power_demand * self.working_time)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CWR_replacement_cost = (
            self.clear_water_tank_cost / self.clear_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.O3_generation_machine_cost / self.O3_generation_machine_lifetime
            ) * scale
        CWR_replacement_cost = CWR_replacement_cost / (365 * 24) # convert to USD/hr
        return CWR_replacement_cost

# %%
PressureTank_path = ospath.join(EL_su_data_path, '_EL_PT.tsv')

@price_ratio()
class EL_PT(StorageTank):
    
    '''
    Introduction
    ------------
    To only collect the clear water within pressure tank

    Parameters
    ----------
    Ins:
    (1) influent of pressurized clear water from clear water tank

    Outs:
    (1) recycle to flush toilet


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    _N_ins = 1; # number of ins
    _N_outs = 1; # number of outs
    _ins_size_is_fixed = True;
    _outs_size_is_fixed = True;


    def __init__(self, ID = '', ins = None, outs = (), vessel_material = None, V_wf = None, include_construction = True,
                 length_to_diameter = None, F_BM_default = 1, kw_per_m3 = None, vessel_type = 'Steel', tau = None,
                 thermo = None, init_with = 'WasteStream', 
                 ppl = None, baseline_ppl = None, **kwargs):
        StorageTank.__init__(self, ID, ins, outs, thermo = thermo, init_with = init_with, include_construction= include_construction,
                             kw_per_m3= kw_per_m3, vessel_type= vessel_type, tau= tau,
                             vessel_material= vessel_material, V_wf= V_wf, F_BM_default= F_BM_default)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.length_to_diameter = length_to_diameter

        data = load_data(path = PressureTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Pressure water tank'] = self.pressure_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        self.power_utility(self.power_demand_PT)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        PW_replacement_cost = (
            self.pressure_water_tank_cost / self.pressure_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        PW_replacement_cost = PW_replacement_cost / (365 * 24) # convert to USD/hr
        return PW_replacement_cost

# %%
blower_path = ospath.join(EL_su_data_path, '_EL_blower.tsv')

@price_ratio()
class EL_blower(Blower):
    '''
    Introduction
    ------------
    To areate air for aerobic tank and membrane tank

    Parameters
    ----------
    Ins:
    (1) air

    Outs:
    (1) air


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.equipments.Blower module

    '''
    _N_ins = 1; # number of ins
    _N_outs = 1; # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True

    def __init__(self, ID = '', ins = None, outs = (), 
                 F_BM={
                     'Blowers': 2.22,
                     'Blower piping': 1,
                     'Blower building': 1.11,
                     },
                 lifetime=15, lifetime_unit='yr',
                 units={
                     'Total gas flow': 'CFM',
                     'Blower capacity': 'CFM',
                     'Number of blowers': '',
                     'Total blower power': 'kW',
                     },
                 N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
                 gas_demand_per_reactor=1, # gas demand per reactor
                 TDH=6, # total dynamic head for rhe blower, in psi
                 eff_blower=0.7, # efficiency of the blower in fraction
                 eff_motor=0.7, # efficiency of the motor in fraction
                 AFF=3.33, # air flow fraction
                 building_unit_cost=9, # unit cost of the building, in USD/ft2
                 thermo = None, ppl = None, baseline_ppl = None, **kwargs):
        Blower.__init__(self, ID=ID, lifetime = lifetime, lifetime_unit = lifetime_unit, F_BM=F_BM,
                        units=units, N_reactor=N_reactor, gas_demand_per_reactor=gas_demand_per_reactor,
                        TDH=TDH, eff_blower=eff_blower, eff_motor=eff_motor, AFF=AFF, building_unit_cost=building_unit_cost,)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.ins = ins
        self.outs = outs
        self.thermo = thermo

        data = load_data(path = blower_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='Steel', linked_unit=self, quantity_unit='kg'),]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.blower_steel_weight; # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Blower'] = self.blower_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        self.power_utility(self.power_demand_blower)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Blower_replacement_cost = (
            self.blower_cost * self.blower_lifetime +
            self.pipeline_connectors / self.pipeline_lifetime +
            self.weld_female_adapter_fittings / self.fittings_lifetime) * scale
        Blower_replacement_cost = Blower_replacement_cost / (365 * 24) # convert to USD/hr
        return Blower_replacement_cost

# %%
housing_path = ospath.join(EL_su_data_path, '_EL_housing.tsv')

@price_ratio()
class EL_Housing(SanUnit):
    '''
     non_reactive unit for the Enviroloo Clear system
    '''

    ppl_per_MURT = 3; # number of people per MURT

    def __init__(self, ID = '', ins = None, outs = (), thermo = None, init_with = 'WasteStream', 
                 ppl = 1, baseline_ppl = None, F_BM_default= 1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo = thermo, init_with = init_with, F_BM_default=F_BM_default)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path = housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self): # replace the actual materials used in the EL
        self.construction = [
            Construction(item = 'Steel', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'Plastic', linked_unit= self, quantity_unit= 'kg'),]

    def _design(self): # replace the actual materials used in the EL
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = (self.steel_weight + self.steel_framework_weight + self.steel_fittings_weight) * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['Plastic'] = constr[1].quantity = (self.LLDPE_weight) * (self.ppl / self.baseline_ppl)   # assume linear scaling
        self.add_construction(add_cost= False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = (self.frame + self.extrusion + 
                        self.angle_frame + self.angle +
                        self.door_sheet + self.plate +
                        self.powder_coating) * (1 + 0.1 * (self.N_EL -1))
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
    
    @property
    def N_EL(self): # determine the number of EL system needed
        return ceil(self.ppl / self.baseline_ppl)
    
    @property
    def N_toilets(self): # determine the number of toilets needed
        return ceil(self.ppl / self.ppl_per_MURT)

# %%
system_path = ospath.join(EL_su_data_path, '_EL_system.tsv')

@price_ratio()
class EL_System(SanUnit):
    '''
    Relate to connection components in the EL system
    '''
    exponent_scale = 0.6

    def __init__(self, ID='', ins=None, outs= (), thermo = None, init_with = 'WasteStream', 
                 if_gridtied = True, ppl = None, baseline_ppl = None, F_BM_default = 1, **kwargs):
        SanUnit.__init__(ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with,  F_BM_default = F_BM_default)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.if_gridtied = if_gridtied

        data = load_data(path = system_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item = 'PVC', linked_unit= self, quantity_unit= 'm'),
            Construction(item = 'HDPE', linked_unit= self, quantity_unit= 'm'),
            ];

    def _design(self):
        design = self.design_results
        design['PVC'] = self.construction[0].quantity = self.PVC_weight * (self.ppl / self.baseline_ppl)
        design['HDPE'] = self.construction[1].quantity = self.HDPE_weight * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost= False)
    
    def _cost(self): # replace these items below by that listed in the _EL_system.tsv
        C = self.baseline_purchase_costs
        C['System'] = (
            self.membrane_filters_M +
            self.membrane_filters_size +
            self.membrane_filters_pause_size +
            self.membrane_filters_chassis_M +
            self.membrane_filters_air_diffuser +
            self.membrane_filters_air_diffuser_chassis +
            self.overflow_membrane2collection +
            self.overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane +
            self.200m_ozone_pipeline +
            self.32mm_pipeline_fittings +
            self.40mm_pipeline_fittings +
            self.60mm_pipeline_fittings +
            self.50mm_ball_valves +
            self.aerobic_air_diffuser)

        ratio = self.price_ratio; # ratio of the price of the new system to the baseline system
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost(); # add the cost of replacement

        if self.if_gridtied:
            power_demand = (self.power_demand_system / 1000) * self.N_EL; # in W/d
        else:
            power_demand = 0
        
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        system_replacement_cost = (
            self.membrane_filters_M / self.lifetime_membrane_filters_M +
            self.membrane_filters_size / self.lifetime_membrane_filters_size +
            self.aerobic_basin / self.lifetime_aerobic_basin +
            self.aerobic_air_diffuser / self.lifetime_aerobic_air_diffuser +
            self.aerobic_air_diffuser_chassis / self.lifetime_aerobic_air_diffuser_chassis +
            self.overflow_membrane2collection / self.lifetime_overflow_membrane2collection +
            self.overflow_clear_water2collection / self.lifetime_overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc / self.lifetime_overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic / self.lifetime_overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane / self.lifetime_overflow_aerobic2membrane +
            self.200m_ozone_pipeline / self.lifetime_200m_ozone_pipeline +
            self.32mm_pipeline_fittings / self.lifetime_32mm_pipeline_fittings +
            self.40mm_pipeline_fittings / self.lifetime_40mm_pipeline_fittings +
            self.60mm_pipeline_fittings / self.lifetime_60mm_pipeline_fittings +
            self.50mm_ball_valves / self.lifetime_50mm_ball_valves +
            self.aerobic_air_diffuser / self.lifetime_aerobic_air_diffuser) * scale
        system_replacement_cost = system_replacement_cost / (365 * 24);  # convert from USD/year to USD/hour
        return system_replacement_cost
    @property
    def N_EL(self): # determine the number of EL system needed
       return ceil(self.ppl / self.baseline_ppl)
