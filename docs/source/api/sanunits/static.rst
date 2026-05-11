Static QSDsan Units
===================

Steady-state QSDsan unit operations — QSDsan-native units without dynamic state equations. They may include sanitation-specific design, costing, construction, TEA, LCA, or ``WasteStream`` behavior.

All classes are also importable directly from :mod:`qsdsan.sanunits`.

Abstract
--------
:class:`~qsdsan.sanunits.PhaseChanger`,
:class:`~qsdsan.sanunits.ComponentSplitter`

Activated Sludge Process
------------------------
:class:`~qsdsan.sanunits.ActivatedSludgeProcess`

Anaerobic Reactors
------------------
:class:`~qsdsan.sanunits.AnaerobicBaffledReactor`,
:class:`~qsdsan.sanunits.AnaerobicDigestion`,
:class:`~qsdsan.sanunits.SludgeDigester`

Biogenic Refinery
-----------------
:class:`~qsdsan.sanunits.BiogenicRefineryCarbonizerBase`,
:class:`~qsdsan.sanunits.BiogenicRefineryControls`,
:class:`~qsdsan.sanunits.BiogenicRefineryGrinder`,
:class:`~qsdsan.sanunits.BiogenicRefineryHHX`,
:class:`~qsdsan.sanunits.BiogenicRefineryHHXdryer`,
:class:`~qsdsan.sanunits.BiogenicRefineryHousing`,
:class:`~qsdsan.sanunits.BiogenicRefineryIonExchange`,
:class:`~qsdsan.sanunits.BiogenicRefineryOHX`,
:class:`~qsdsan.sanunits.BiogenicRefineryPollutionControl`,
:class:`~qsdsan.sanunits.BiogenicRefineryScrewPress`,
:class:`~qsdsan.sanunits.BiogenicRefineryStruvitePrecipitation`

Clarifiers
----------
:class:`~qsdsan.sanunits.FlatBottomCircularClarifier`,
:class:`~qsdsan.sanunits.IdealClarifier`,
:class:`~qsdsan.sanunits.PrimaryClarifier`,
:class:`~qsdsan.sanunits.PrimaryClarifierBSM2`

Combustion
----------
:class:`~qsdsan.sanunits.BiogasCombustion`,
:class:`~qsdsan.sanunits.CombinedHeatPower`

Crop Application
----------------
:class:`~qsdsan.sanunits.CropApplication`

Eco-San
-------
:class:`~qsdsan.sanunits.EcoSanAerobic`,
:class:`~qsdsan.sanunits.EcoSanAnaerobic`,
:class:`~qsdsan.sanunits.EcoSanAnoxic`,
:class:`~qsdsan.sanunits.EcoSanBioCost`,
:class:`~qsdsan.sanunits.EcoSanECR`,
:class:`~qsdsan.sanunits.EcoSanMBR`,
:class:`~qsdsan.sanunits.EcoSanPrimary`,
:class:`~qsdsan.sanunits.EcoSanSolar`,
:class:`~qsdsan.sanunits.EcoSanSystem`

Electrochemical Cell
--------------------
:class:`~qsdsan.sanunits.ElectrochemicalCell`

Excretion
---------
:class:`~qsdsan.sanunits.Excretion`,
:class:`~qsdsan.sanunits.ExcretionmASM2d`

Hydroprocessing
---------------
:class:`~qsdsan.sanunits.Hydrocracking`,
:class:`~qsdsan.sanunits.Hydrotreating`

Hydrothermal
------------
:class:`~qsdsan.sanunits.CatalyticHydrothermalGasification`,
:class:`~qsdsan.sanunits.HydrothermalLiquefaction`

Internal Circulation Reactor
----------------------------
:class:`~qsdsan.sanunits.InternalCirculationRx`

Lagoon
------
:class:`~qsdsan.sanunits.Lagoon`

Membrane
--------
:class:`~qsdsan.sanunits.MembraneDistillation`,
:class:`~qsdsan.sanunits.GasExtractionMembrane`,
:class:`~qsdsan.sanunits.MembraneGasExtraction`

Metal Dosage
------------
:class:`~qsdsan.sanunits.MetalDosage`

Non-Reactive
------------
:class:`~qsdsan.sanunits.Copier`,
:class:`~qsdsan.sanunits.LumpedCost`

Polishing Filter
----------------
:class:`~qsdsan.sanunits.PolishingFilter`

Pumping (WWT-specific)
----------------------
:class:`~qsdsan.sanunits.WWTpump`,
:class:`~qsdsan.sanunits.SludgePump`,
:func:`~qsdsan.sanunits.wwtpump`

Reactor
-------
:class:`~qsdsan.sanunits.Reactor`

Reclaimer
---------
:class:`~qsdsan.sanunits.ReclaimerECR`,
:class:`~qsdsan.sanunits.ReclaimerHousing`,
:class:`~qsdsan.sanunits.ReclaimerIonExchange`,
:class:`~qsdsan.sanunits.ReclaimerSolar`,
:class:`~qsdsan.sanunits.ReclaimerSystem`,
:class:`~qsdsan.sanunits.ReclaimerUltrafiltration`

Screening
---------
:class:`~qsdsan.sanunits.Screening`

Sedimentation
-------------
:class:`~qsdsan.sanunits.Sedimentation`

Septic Tank
-----------
:class:`~qsdsan.sanunits.SepticTank`

Sludge Pasteurization
---------------------
:class:`~qsdsan.sanunits.SludgePasteurization`

Sludge Thickening
-----------------
:class:`~qsdsan.sanunits.SludgeThickening`,
:class:`~qsdsan.sanunits.BeltThickener`,
:class:`~qsdsan.sanunits.SludgeCentrifuge`,
:class:`~qsdsan.sanunits.SludgeSeparator`

Sludge Treatment
----------------
:class:`~qsdsan.sanunits.Thickener`,
:class:`~qsdsan.sanunits.Centrifuge`,
:class:`~qsdsan.sanunits.Incinerator`

Toilets
-------
:class:`~qsdsan.sanunits.Toilet`,
:class:`~qsdsan.sanunits.MURT`,
:class:`~qsdsan.sanunits.PitLatrine`,
:class:`~qsdsan.sanunits.UDDT`

Treatment Bed
-------------
:class:`~qsdsan.sanunits.DryingBed`,
:class:`~qsdsan.sanunits.LiquidTreatmentBed`

Trucking
--------
:class:`~qsdsan.sanunits.Trucking`
