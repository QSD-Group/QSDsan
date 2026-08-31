.. _jupyter_tips:

Jupyter Tips
============

These tutorials are Jupyter notebooks, and a few of Jupyter's built-in features make
them much easier to follow — especially for reading documentation and discovering what a
``QSDsan`` object can do. If you are new to Jupyter itself (cells, kernels, running code),
start with the `Project Jupyter documentation <https://docs.jupyter.org/en/latest/>`_;
the tips below assume you can already run a notebook.


Reading documentation
----------------------

Every ``QSDsan`` class and method carries a docstring. You can read it without leaving
the notebook:

* ``obj?`` — show the docstring (and signature) in a pop-up pager. This is the form used
  throughout the tutorials, e.g. ``qs.WasteStream?`` or ``ws.get_VSS?``.
* ``obj??`` — same as ``?`` but also shows the **source code**, when it is written in
  Python. Handy for seeing exactly what a method does.
* ``help(obj)`` — print the docstring inline in the cell output. This is plain Python, so
  it works outside Jupyter too.
* **Shift+Tab** — with the cursor inside a call's parentheses (e.g.
  ``qs.WasteStream(|)``), pop up a signature/docstring tooltip; press it again to expand.
  (JupyterLab and Jupyter Notebook only.)
* ``obj.__doc__`` — the raw docstring as a string, if you want to handle it
  programmatically.

.. code-block:: ipython3
   :force:

    import qsdsan as qs

    qs.WasteStream?              # docstring of the class
    qs.WasteStream.composite?    # docstring of one method
    qs.WasteStream.composite??   # docstring + source code
    help(qs.WasteStream.composite)

.. note::

   Many ``QSDsan`` classes inherit methods from ``BioSTEAM``/``Thermosteam`` (e.g.
   ``WasteStream`` from ``Stream``). ``?`` still shows the inherited docstring, and
   ``??`` reveals which parent module the method actually lives in.


Discovering what's available
----------------------------

Before looking up *how* a method works, you often want to know *which* attributes and
methods exist:

* **Tab completion** — type ``obj.`` and press **Tab** to list available attributes and
  methods. Works on modules too: ``qs.unit_operations.`` + Tab lists the unit operations.
* ``dir(obj)`` — return the same list as a Python list (useful when completion is slow or
  you want to filter it in code). Names with a leading underscore are considered internal.

.. code-block:: python

    dir(qs.unit_operations)              # all available unit operations
    [m for m in dir(ws) if 'COD' in m]   # methods/attributes mentioning COD


Other handy bits
----------------

* **Restart the kernel** (*Kernel → Restart*) to clear all state and start fresh — useful
  after changing imports, and required once after installing ``QSDsan`` in a fresh session
  (see the Colab instructions in the tutorials index).
* A cell's **last expression** is displayed automatically; you do not need ``print`` to
  see it. Use ``print`` only when you want to show several values from one cell.
