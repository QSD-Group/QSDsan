Dynamic Units
=============

Dynamic QSDsan unit operations include explicit dynamic-state behavior such as ``state``, ``dstate``, ODE compilation, or algebraic-equation compilation. Use this page when building systems for dynamic simulation.

All classes are also importable directly from :mod:`qsdsan.unit_operations`.

.. note::
   :class:`~qsdsan.unit_operations.Mixer`, :class:`~qsdsan.unit_operations.Splitter`, and :class:`~qsdsan.unit_operations.Pump` below wrap BioSTEAM units (see the :doc:`BioSTEAM-inherited units page <../bst/index>`), but live here rather than there because they support dynamic simulation.

.. toctree::
   :maxdepth: 1

   abstract
   anaerobic_reactor
   clarifier
   DynamicInfluent
   ElectrochemicalCell
   Excretion
   junction
   membrane_bioreactor
   MembraneGasExtraction
   MetalDosage
   pumping
   sludge_treatment
   suspended_growth_bioreactor

Abstract
--------
:class:`~qsdsan.unit_operations.Sampler`,
:class:`~qsdsan.unit_operations.HydraulicDelay`,
:class:`~qsdsan.unit_operations.Mixer`,
:class:`~qsdsan.unit_operations.Splitter`

Anaerobic Reactors
------------------
:class:`~qsdsan.unit_operations.AnaerobicCSTR`

Bioreactors
-----------
:class:`~qsdsan.unit_operations.CSTR`,
:class:`~qsdsan.unit_operations.BatchExperiment`,
:class:`~qsdsan.unit_operations.PFR`,
:class:`~qsdsan.unit_operations.AerobicDigester`

Clarifiers
----------
:class:`~qsdsan.unit_operations.FlatBottomCircularClarifier`,
:class:`~qsdsan.unit_operations.IdealClarifier`,
:class:`~qsdsan.unit_operations.PrimaryClarifier`,
:class:`~qsdsan.unit_operations.PrimaryClarifierBSM2`

Electrochemical Cell
--------------------
:class:`~qsdsan.unit_operations.ESAP`,
:class:`~qsdsan.unit_operations.ESAPRecovery`,
:class:`~qsdsan.unit_operations.ESAPEffluent`

Excretion
---------
:class:`~qsdsan.unit_operations.ExcretionmASM2d`

Influent
--------
:class:`~qsdsan.unit_operations.DynamicInfluent`

Junctions
---------
:class:`~qsdsan.unit_operations.Junction`,
:class:`~qsdsan.unit_operations.ADMjunction`,
:class:`~qsdsan.unit_operations.mADMjunction`,
:class:`~qsdsan.unit_operations.ADMtoASM`,
:class:`~qsdsan.unit_operations.ASMtoADM`,
:class:`~qsdsan.unit_operations.ASM2dtoADM1`,
:class:`~qsdsan.unit_operations.ADM1toASM2d`,
:class:`~qsdsan.unit_operations.ASM2dtomADM1`,
:class:`~qsdsan.unit_operations.mADM1toASM2d`,
:class:`~qsdsan.unit_operations.A1junction`,
:class:`~qsdsan.unit_operations.ADM1ptomASM2d`,
:class:`~qsdsan.unit_operations.mASM2dtoADM1p`

Membrane Bioreactors
--------------------
:class:`~qsdsan.unit_operations.AnMBR`,
:class:`~qsdsan.unit_operations.CompletelyMixedMBR`

Membrane Gas Extraction
-----------------------
:class:`~qsdsan.unit_operations.GasExtractionMembrane`,
:class:`~qsdsan.unit_operations.MembraneGasExtraction`

Metal Dosage
------------
:class:`~qsdsan.unit_operations.MetalDosage`

Pumping
-------
:class:`~qsdsan.unit_operations.Pump`

Sludge Treatment
----------------
:class:`~qsdsan.unit_operations.Thickener`,
:class:`~qsdsan.unit_operations.Centrifuge`,
:class:`~qsdsan.unit_operations.Incinerator`

