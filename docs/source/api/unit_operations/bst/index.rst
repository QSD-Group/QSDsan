BioSTEAM-Inherited Units
========================

BioSTEAM-inherited unit operations are steady-state units that build on BioSTEAM classes and add QSDsan behavior. Use this page when you want familiar BioSTEAM-style unit operations with QSDsan streams, construction, costing, or other QSDsan extensions.

All classes are also importable directly from :mod:`qsdsan.unit_operations`.

.. toctree::
   :maxdepth: 1

   abstract
   compressor
   distillation
   facilities
   Flash
   heat_exchanging
   pumping
   tank

Abstract
--------
:class:`~qsdsan.unit_operations.Mixer`,
:class:`~qsdsan.unit_operations.Splitter`,
:class:`~qsdsan.unit_operations.FakeSplitter`,
:class:`~qsdsan.unit_operations.ReversedSplitter`

Compressor
----------
:class:`~qsdsan.unit_operations.IsothermalCompressor`

Distillation
------------
:class:`~qsdsan.unit_operations.BinaryDistillation`,
:class:`~qsdsan.unit_operations.ShortcutColumn`,
:class:`~qsdsan.unit_operations.MESHDistillation`,
:class:`~qsdsan.unit_operations.AdiabaticMultiStageVLEColumn`

Facilities
----------
:class:`~qsdsan.unit_operations.ProcessWaterCenter`

Flash
-----
:class:`~qsdsan.unit_operations.Flash`

Heat Exchanging
---------------
:class:`~qsdsan.unit_operations.HeatExchangerNetwork`,
:class:`~qsdsan.unit_operations.HXprocess`,
:class:`~qsdsan.unit_operations.HXutility`

Pumping
-------
:class:`~qsdsan.unit_operations.Pump`

Tank
----
:class:`~qsdsan.unit_operations.Tank`,
:class:`~qsdsan.unit_operations.MixTank`,
:class:`~qsdsan.unit_operations.StorageTank`
