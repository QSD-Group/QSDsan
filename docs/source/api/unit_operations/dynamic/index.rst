Dynamic QSDsan Units
====================

Dynamic QSDsan unit operations include explicit dynamic-state behavior such as ``state``, ``dstate``, ODE compilation, or algebraic-equation compilation. Use this page when building systems for dynamic simulation.

All classes are also importable directly from :mod:`qsdsan.sanunits`.

Abstract
--------
:class:`~qsdsan.unit_operations.Sampler`

Anaerobic Reactors
------------------
:class:`~qsdsan.unit_operations.AnaerobicCSTR`

Bioreactors
-----------
:class:`~qsdsan.unit_operations.CSTR`,
:class:`~qsdsan.unit_operations.BatchExperiment`,
:class:`~qsdsan.unit_operations.PFR`,
:class:`~qsdsan.unit_operations.AerobicDigester`

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

Pumping
-------
:class:`~qsdsan.unit_operations.HydraulicDelay`
