BioSTEAM-Inherited Units
========================

BioSTEAM-inherited unit operations are steady-state units that build on BioSTEAM classes and add QSDsan behavior. Use this page when you want familiar BioSTEAM-style unit operations with QSDsan streams, construction, costing, or other QSDsan extensions.

All classes are also importable directly from :mod:`qsdsan.unit_operations`.

.. note::
   Looking for :class:`~qsdsan.unit_operations.Mixer`, :class:`~qsdsan.unit_operations.Splitter`, or :class:`~qsdsan.unit_operations.Pump`? They also wrap BioSTEAM units, but because they support dynamic simulation they're documented on the :doc:`dynamic units page <../dynamic/index>` instead (Abstract and Pumping sections). :class:`~qsdsan.unit_operations.FakeSplitter` and :class:`~qsdsan.unit_operations.ReversedSplitter` below have no dynamic-state support, so they stay here.

.. toctree::
   :maxdepth: 1

   abstract
   compressor
   distillation
   facilities
   Flash
   heat_exchanging
   tank

Abstract
--------
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

Tank
----
:class:`~qsdsan.unit_operations.Tank`,
:class:`~qsdsan.unit_operations.MixTank`,
:class:`~qsdsan.unit_operations.StorageTank`
