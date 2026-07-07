Static Units
============

Steady-state QSDsan unit operations — QSDsan-native units without dynamic state equations. They may include sanitation-specific design, costing, construction, TEA, LCA, or ``WasteStream`` behavior.

All classes are also importable directly from :mod:`qsdsan.unit_operations`.

.. toctree::
   :maxdepth: 1

   abstract
   ActivatedSludgeProcess
   anaerobic_reactor
   biogenic_refinery
   combustion
   CropApplication
   eco_san
   ElectrochemicalCell
   Excretion
   g2rt
   hydroprocessing
   hydrothermal
   InternalCirculationRx
   Lagoon
   MembraneDistillation
   non_reactive
   PolishingFilter
   pumping
   Reactor
   reclaimer
   Screening
   Sedimentation
   SepticTank
   sludge_thickening
   SludgePasteurization
   toilet
   treatment_bed
   Trucking

Abstract
--------
:class:`~qsdsan.unit_operations.PhaseChanger`,
:class:`~qsdsan.unit_operations.ComponentSplitter`

Activated Sludge Process
------------------------
:class:`~qsdsan.unit_operations.ActivatedSludgeProcess`

Anaerobic Reactors
------------------
:class:`~qsdsan.unit_operations.AnaerobicBaffledReactor`,
:class:`~qsdsan.unit_operations.AnaerobicDigestion`,
:class:`~qsdsan.unit_operations.SludgeDigester`

Biogenic Refinery
-----------------
:class:`~qsdsan.unit_operations.BiogenicRefineryCarbonizerBase`,
:class:`~qsdsan.unit_operations.BiogenicRefineryControls`,
:class:`~qsdsan.unit_operations.BiogenicRefineryGrinder`,
:class:`~qsdsan.unit_operations.BiogenicRefineryHHX`,
:class:`~qsdsan.unit_operations.BiogenicRefineryHHXdryer`,
:class:`~qsdsan.unit_operations.BiogenicRefineryHousing`,
:class:`~qsdsan.unit_operations.BiogenicRefineryIonExchange`,
:class:`~qsdsan.unit_operations.BiogenicRefineryOHX`,
:class:`~qsdsan.unit_operations.BiogenicRefineryPollutionControl`,
:class:`~qsdsan.unit_operations.BiogenicRefineryScrewPress`,
:class:`~qsdsan.unit_operations.BiogenicRefineryStruvitePrecipitation`

Combustion
----------
:class:`~qsdsan.unit_operations.BiogasCombustion`,
:class:`~qsdsan.unit_operations.CombinedHeatPower`

Crop Application
----------------
:class:`~qsdsan.unit_operations.CropApplication`

Eco-San
-------
:class:`~qsdsan.unit_operations.EcoSanAerobic`,
:class:`~qsdsan.unit_operations.EcoSanAnaerobic`,
:class:`~qsdsan.unit_operations.EcoSanAnoxic`,
:class:`~qsdsan.unit_operations.EcoSanBioCost`,
:class:`~qsdsan.unit_operations.EcoSanECR`,
:class:`~qsdsan.unit_operations.EcoSanMBR`,
:class:`~qsdsan.unit_operations.EcoSanPrimary`,
:class:`~qsdsan.unit_operations.EcoSanSolar`,
:class:`~qsdsan.unit_operations.EcoSanSystem`

Electrochemical Cell
--------------------
:class:`~qsdsan.unit_operations.ElectrochemicalCell`,
:class:`~qsdsan.unit_operations.ElectrochemicalStrippingAdsorptionPrecipitation`

Excretion
---------
:class:`~qsdsan.unit_operations.Excretion`

G2RT
----
:class:`~qsdsan.unit_operations.FWMixer`,
:class:`~qsdsan.unit_operations.G2RTBeltSeparation`,
:class:`~qsdsan.unit_operations.G2RTControls`,
:class:`~qsdsan.unit_operations.G2RTExcretion`,
:class:`~qsdsan.unit_operations.G2RThomogenizer`,
:class:`~qsdsan.unit_operations.G2RTHousing`,
:class:`~qsdsan.unit_operations.G2RTLiquidsTank`,
:class:`~qsdsan.unit_operations.G2RTReverseOsmosis`,
:class:`~qsdsan.unit_operations.G2RTSolidsSeparation`,
:class:`~qsdsan.unit_operations.G2RTSolidsTank`,
:class:`~qsdsan.unit_operations.G2RTUltrafiltration`,
:class:`~qsdsan.unit_operations.mSCWOConcentratorModule`,
:class:`~qsdsan.unit_operations.mSCWOGasModule`,
:class:`~qsdsan.unit_operations.mSCWOReactorModule`,
:class:`~qsdsan.unit_operations.UFMixer`,
:class:`~qsdsan.unit_operations.VolumeReductionCombustor`,
:class:`~qsdsan.unit_operations.VolumeReductionFilterPress`,
:class:`~qsdsan.unit_operations.VRConcentrator`,
:class:`~qsdsan.unit_operations.VRdryingtunnel`,
:class:`~qsdsan.unit_operations.VRpasteurization`

Hydroprocessing
---------------
:class:`~qsdsan.unit_operations.Hydrocracking`,
:class:`~qsdsan.unit_operations.Hydrotreating`

Hydrothermal
------------
:class:`~qsdsan.unit_operations.CatalyticHydrothermalGasification`,
:class:`~qsdsan.unit_operations.HydrothermalLiquefaction`

Internal Circulation Reactor
----------------------------
:class:`~qsdsan.unit_operations.InternalCirculationRx`

Lagoon
------
:class:`~qsdsan.unit_operations.Lagoon`

Membrane
--------
:class:`~qsdsan.unit_operations.MembraneDistillation`

Non-Reactive
------------
:class:`~qsdsan.unit_operations.Copier`,
:class:`~qsdsan.unit_operations.LumpedCost`

Polishing Filter
----------------
:class:`~qsdsan.unit_operations.PolishingFilter`

Pumping (WWT-specific)
----------------------
:class:`~qsdsan.unit_operations.WWTpump`,
:class:`~qsdsan.unit_operations.SludgePump`,
:func:`~qsdsan.unit_operations.wwtpump`

Reactor
-------
:class:`~qsdsan.unit_operations.Reactor`

Reclaimer
---------
:class:`~qsdsan.unit_operations.ReclaimerECR`,
:class:`~qsdsan.unit_operations.ReclaimerHousing`,
:class:`~qsdsan.unit_operations.ReclaimerIonExchange`,
:class:`~qsdsan.unit_operations.ReclaimerSolar`,
:class:`~qsdsan.unit_operations.ReclaimerSystem`,
:class:`~qsdsan.unit_operations.ReclaimerUltrafiltration`

Screening
---------
:class:`~qsdsan.unit_operations.Screening`

Sedimentation
-------------
:class:`~qsdsan.unit_operations.Sedimentation`

Septic Tank
-----------
:class:`~qsdsan.unit_operations.SepticTank`

Sludge Pasteurization
---------------------
:class:`~qsdsan.unit_operations.SludgePasteurization`

Sludge Thickening
-----------------
:class:`~qsdsan.unit_operations.SludgeThickening`,
:class:`~qsdsan.unit_operations.BeltThickener`,
:class:`~qsdsan.unit_operations.SludgeCentrifuge`,
:class:`~qsdsan.unit_operations.SludgeSeparator`

Toilets
-------
:class:`~qsdsan.unit_operations.Toilet`,
:class:`~qsdsan.unit_operations.ReinventedToilet`,
:class:`~qsdsan.unit_operations.MURT`,
:class:`~qsdsan.unit_operations.SURT`,
:class:`~qsdsan.unit_operations.PitLatrine`,
:class:`~qsdsan.unit_operations.UDDT`

Treatment Bed
-------------
:class:`~qsdsan.unit_operations.DryingBed`,
:class:`~qsdsan.unit_operations.LiquidTreatmentBed`

Trucking
--------
:class:`~qsdsan.unit_operations.Trucking`
