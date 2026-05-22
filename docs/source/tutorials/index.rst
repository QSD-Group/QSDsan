.. _tutorials:

Tutorials
=========

These tutorials teach you how to use ``QSDsan``: defining components, building waste streams,
assembling units into systems, and running techno-economic analysis (TEA), life cycle assessment (LCA), and uncertainty/sensitivity analyses.

All tutorials are Jupyter notebooks. Run them locally in your own Jupyter
environment, or open them in your browser by clicking the badge below:

.. image:: ../images/custom_binder_logo.svg
   :target: https://mybinder.org/v2/gh/QSD-Group/QSDsan-env/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FQSD-group%252FQSDsan%26urlpath%3Dlab%252Ftree%252FQSDsan%252Fdocs%252Fsource%252Ftutorials%26branch%3Dmain
   :alt: Launch Binder


.. _run-in-colab:

Running tutorials in Google Colab
---------------------------------

You can also run any tutorial in `Google Colab <https://colab.research.google.com>`_
without installing anything locally. Colab does not include ``QSDsan`` in its
default environment, so a short one-time setup is needed for each session:

#. Open the notebook in Colab: go to **File → Open notebook → GitHub**, search for
   ``QSD-Group/QSDsan``, and pick the tutorial you want (under
   ``docs/source/tutorials``).
#. Add a new cell at the very top (or use the tutorial's first cell) and install
   ``QSDsan`` and ``EXPOsan``:

   .. code-block:: text

      !pip install qsdsan exposan

#. Run that cell. Colab will warn that some pre-installed packages (such as
   ``numpy`` and ``matplotlib``) were already imported and that you must restart
   the runtime. This is expected, not an error: installing ``QSDsan`` uses different
   versions of those packages, but the versions already loaded into the running session
   cannot be swapped until you restart.
#. Restart the session via **Runtime → Restart session**.
#. Continue running the notebook from the cell *after* the install cell. Do not
   re-run the install cell.

.. note::
   The restart is needed only once per session. Do not import any packages
   before installing, restart once after installing, 
   and then skip the install cell, the warning will not reappear.


.. note::
   **About the YouTube walkthroughs.** Some tutorials have companion videos on
   our `YouTube channel
   <https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA>`_. The videos were
   recorded against earlier versions of ``QSDsan`` and remain useful for the
   concepts and the big picture, but the notebooks here are the authoritative
   reference for syntax and API. Each video description lists the ``QSDsan``
   version it was filmed against.


New to Python or Jupyter?
-------------------------

These tutorials assume you can read Python code, run a Jupyter notebook, and
have ``QSDsan`` installed. If you are new to Python or scientific computing,
work through one of the resources below first;
they cover the fundamentals far better than we could here.

General Python:

- `The official Python tutorial <https://docs.python.org/3/tutorial/>`_, the
  canonical free introduction to the language.
- `Python for Everybody <https://www.py4e.com/>`_ by Charles Severance, a
  textbook and video course aimed at absolute beginners.
- `Real Python <https://realpython.com/>`_, searchable tutorials on individual
  topics once you have the basics.

Python for researchers and engineers:

- `Software Carpentry: Plotting and Programming in Python
  <https://swcarpentry.github.io/python-novice-gapminder/>`_, a short hands-on
  lesson series designed for scientists.
- `Scientific Python Lectures <https://lectures.scientific-python.org/>`_, a
  free course covering NumPy, SciPy, matplotlib, and pandas.

Jupyter:

- `Project Jupyter documentation <https://docs.jupyter.org/en/latest/>`_ and
  `Try Jupyter <https://jupyter.org/try>`_, the official docs and an
  in-browser environment to experiment in.

Throughout the topical tutorials you will also find collapsible **Python Aside**
callouts. These explain a Python concept or idiom that comes up in the surrounding
code, so you can pick it up in context. Expand one if you want the explanation, or
skip it if you already know the concept.

Already comfortable with Python but new to these notebooks? The
:doc:`Jupyter tips <jupyter_tips>` page covers the few Jupyter features that you might
find handy when working with these tutorials.

.. toctree::
   :hidden:

   jupyter_tips

For installing ``QSDsan`` itself, see the installation section on the
:doc:`main documentation page </index>`.


Topical Tutorials
-----------------

The topical tutorials are organized in three parts.

**Part I.** A class-agnostic tour of ``QSDsan``.

.. toctree::
   :maxdepth: 1

   0_Quick_Overview

**Part II.** QSDsan's core classes. Each tutorial focuses on one or a few core classes.

.. toctree::
   :maxdepth: 1

   2_Component
   3_WasteStream
   4_SanUnit_basic
   5_SanUnit_advanced
   6_System
   7_TEA
   8_LCA
   9_Uncertainty_and_Sensitivity_Analyses
   10_Process

**Part III.** Process modeling and dynamic simulation.

.. toctree::
   :maxdepth: 1

   11_Dynamic_Simulation
   12_Anaerobic_Digestion_Model_No_1
   13_Process_Modeling_101


For Contributors
----------------

.. toctree::
   :maxdepth: 1

   14_AI_Assisted_Development


.. 1_Helpful_Basics is kept as an orphan deprecation stub so its URL still
   resolves for anyone with a bookmark; it is excluded from this index and
   from the global navigation via its own ``nbsphinx.orphan`` metadata.


Additional Resources
--------------------

``QSDsan`` is built on `BioSTEAM <https://biosteam.readthedocs.io/en/latest/index.html>`_,
an open-source platform for process modeling and techno-economic analysis focusing on biorefineries.
Some of QSDsan's core classes and methods are inherited from BioSTEAM. 
You may want to refer to BioSTEAM's documentation for the parent classes and methods.