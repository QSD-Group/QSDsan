# -*- coding: utf-8 -*-
"""Steady-state QSDsan unit operations."""

from ._static_abstract import PhaseChanger, ComponentSplitter
from ._static_activated_sludge_process import ActivatedSludgeProcess
from ._static_anaerobic_reactor import (
    AnaerobicBaffledReactor, AnaerobicDigestion, SludgeDigester,
)
from ._static_biogenic_refinery import (
    BiogenicRefineryCarbonizerBase, BiogenicRefineryControls,
    BiogenicRefineryGrinder, BiogenicRefineryHHX, BiogenicRefineryHHXdryer,
    BiogenicRefineryHousing, BiogenicRefineryIonExchange,
    BiogenicRefineryOHX, BiogenicRefineryPollutionControl,
    BiogenicRefineryScrewPress, BiogenicRefineryStruvitePrecipitation,
)
from ._static_clarifier import (
    FlatBottomCircularClarifier, IdealClarifier,
    PrimaryClarifierBSM2, PrimaryClarifier,
)
from ._static_combustion import BiogasCombustion, CombinedHeatPower
from ._static_crop_application import CropApplication
from ._static_eco_san import (
    EcoSanAerobic, EcoSanAnaerobic, EcoSanAnoxic, EcoSanBioCost,
    EcoSanECR, EcoSanMBR, EcoSanPrimary, EcoSanSolar, EcoSanSystem,
)
from ._static_electrochemical_cell import ElectrochemicalCell
from ._static_excretion import Excretion, ExcretionmASM2d
from ._static_hydroprocessing import Hydrocracking, Hydrotreating
from ._static_hydrothermal import CatalyticHydrothermalGasification, HydrothermalLiquefaction
from ._static_internal_circulation_rx import InternalCirculationRx
from ._static_lagoon import Lagoon
from ._static_membrane_distillation import MembraneDistillation
from ._static_membrane_gas_extraction import GasExtractionMembrane, MembraneGasExtraction
from ._static_metal_dosage import MetalDosage
from ._static_non_reactive import Copier, LumpedCost
from ._static_polishing_filter import PolishingFilter
from ._static_pumping import WWTpump, SludgePump, wwtpump
from ._static_reactor import Reactor
from ._static_reclaimer import (
    ReclaimerECR, ReclaimerHousing, ReclaimerIonExchange,
    ReclaimerSolar, ReclaimerSystem, ReclaimerUltrafiltration,
)
from ._static_screening import Screening
from ._static_sedimentation import Sedimentation
from ._static_septic_tank import SepticTank
from ._static_sludge_pasteurization import SludgePasteurization
from ._static_sludge_thickening import (
    SludgeThickening, BeltThickener, SludgeCentrifuge, SludgeSeparator,
)
from ._static_sludge_treatment import Thickener, Centrifuge, Incinerator
from ._static_toilet import Toilet, MURT, PitLatrine, UDDT
from ._static_treatment_bed import DryingBed, LiquidTreatmentBed
from ._static_trucking import Trucking

__all__ = (
    'PhaseChanger', 'ComponentSplitter',
    'ActivatedSludgeProcess',
    'AnaerobicBaffledReactor', 'AnaerobicDigestion', 'SludgeDigester',
    'BiogenicRefineryCarbonizerBase', 'BiogenicRefineryControls',
    'BiogenicRefineryGrinder', 'BiogenicRefineryHHX', 'BiogenicRefineryHHXdryer',
    'BiogenicRefineryHousing', 'BiogenicRefineryIonExchange',
    'BiogenicRefineryOHX', 'BiogenicRefineryPollutionControl',
    'BiogenicRefineryScrewPress', 'BiogenicRefineryStruvitePrecipitation',
    'FlatBottomCircularClarifier', 'IdealClarifier',
    'PrimaryClarifierBSM2', 'PrimaryClarifier',
    'BiogasCombustion', 'CombinedHeatPower',
    'CropApplication',
    'EcoSanAerobic', 'EcoSanAnaerobic', 'EcoSanAnoxic', 'EcoSanBioCost',
    'EcoSanECR', 'EcoSanMBR', 'EcoSanPrimary', 'EcoSanSolar', 'EcoSanSystem',
    'ElectrochemicalCell',
    'Excretion', 'ExcretionmASM2d',
    'Hydrocracking', 'Hydrotreating',
    'CatalyticHydrothermalGasification', 'HydrothermalLiquefaction',
    'InternalCirculationRx',
    'Lagoon',
    'MembraneDistillation',
    'GasExtractionMembrane', 'MembraneGasExtraction',
    'MetalDosage',
    'Copier', 'LumpedCost',
    'PolishingFilter',
    'WWTpump', 'SludgePump', 'wwtpump',
    'Reactor',
    'ReclaimerECR', 'ReclaimerHousing', 'ReclaimerIonExchange',
    'ReclaimerSolar', 'ReclaimerSystem', 'ReclaimerUltrafiltration',
    'Screening',
    'Sedimentation',
    'SepticTank',
    'SludgePasteurization',
    'SludgeThickening', 'BeltThickener', 'SludgeCentrifuge', 'SludgeSeparator',
    'Thickener', 'Centrifuge', 'Incinerator',
    'Toilet', 'MURT', 'PitLatrine', 'UDDT',
    'DryingBed', 'LiquidTreatmentBed',
    'Trucking',
)
