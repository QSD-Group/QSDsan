Streams
=======

Which stream class to use?
--------------------------
``QSDsan`` can work with three main stream classes: :class:`thermosteam.Stream`, :class:`~.SanStream`, and :class:`WasteStream`. Follow this flowchart to decide which class you want to use.

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


.. note::
	
	Regardless of which stream class you choose, you always use :class:`~.Component` and :class:`~.Components` instead of :class:`thermosteam.Chemical` and :class:`thermosteam.Chemicals` to indicate which components to include in your system, this is true even if you choose to use :class:`thermosteam.Stream` and all of the components are pure chemicals. Although in that case, it might be more straightforward to use ``BioSTEAM`` instead of ``QSDsan`` unless you want to include environmental impacts in your analyses.


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