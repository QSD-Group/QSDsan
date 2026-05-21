Streams
=======

Which stream class to use?
--------------------------
``QSDsan`` can work with three main stream classes: :class:`biosteam.Stream`, :class:`~.SanStream`, and :class:`~.WasteStream`. They form a capability hierarchy:

* :class:`biosteam.Stream` is the general material stream class from ``BioSTEAM``.
* :class:`~.SanStream` adds stream-level life cycle impact functionality.
* :class:`~.WasteStream` adds wastewater-modeling functionality on top of :class:`~.SanStream`.

Follow this flowchart to decide which class you want to use.

.. figure:: https://lucid.app/publicSegments/view/5aae11af-e8cf-434a-939d-9abe07cb1b82/image.png
   :width: 500
   :align: center


**Essentially,**

+-------------------------------------+-----------------------------+----------------------+----------------------+
|                                     | :class:`thermosteam.Stream` | :class:`SanStream`   | :class:`WasteStream` |
+=====================================+=============================+======================+======================+
| Capable of doing LCA?               | No                          | Yes                  | Yes                  |
+-------------------------------------+-----------------------------+----------------------+----------------------+
| Has wastewater-specific properties? | No                          | No                   | Yes                  |
+-------------------------------------+-----------------------------+----------------------+----------------------+


If you are unsure:

* Use :class:`biosteam.Stream` when you only need flow, composition, thermodynamic, and price information.
* Use :class:`~.SanStream` when you also need stream-level environmental impact accounting.
* Use :class:`~.WasteStream` when you need wastewater-specific quantities such as COD, BOD, TKN, TP, solids, pH, biodegradability fractions, or influent characterization models.

Despite the name, :class:`~.WasteStream` does not require the material to be discarded as waste. It is ``QSDsan``'s wastewater-modeling stream class: use it when the stream is represented using wastewater treatment model states or when you need aggregate wastewater properties and influent characterization methods.


.. note::
	
	Regardless of which stream class you choose, you always use :class:`~.Component` and :class:`~.Components` instead of :class:`biosteam.Chemical` and :class:`biosteam.Chemicals` to indicate which components to include in your system, this is true even if you choose to use :class:`biosteam.Stream` and all of the components are pure chemicals. Although in that case, it might be more straightforward to use ``BioSTEAM`` instead of ``QSDsan`` unless you want to include environmental impacts in your analyses.


SanStream
---------

.. autoclass:: qsdsan.SanStream
   :members:

MissingSanStream
----------------

.. autoclass:: qsdsan.MissingSanStream
   :members:

WasteStream
-----------

.. autoclass:: qsdsan.WasteStream
   :members:

MissingWasteStream
------------------

.. autoclass:: qsdsan.MissingWasteStream
   :members:
