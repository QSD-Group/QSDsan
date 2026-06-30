Process Models
==============

This page lists the process models currently exposed by ``qsdsan.process_models`` for biological or chemical kinetic simulation, along with example implementations in `EXPOsan`_ and the original references.

.. note::

   ``qsdsan.processes`` is a legacy alias for ``qsdsan.process_models``. Both names point to the same package, but ``process_models`` is the current, preferred name.


Biological kinetic models
-------------------------

Each model below is a subclass of :class:`qsdsan.CompiledProcesses` and represents a set of parallel transformation processes (stoichiometry + rate equations) for the listed components.

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Model
     - Implementation
     - Reference
   * - ADM1
     - `adm`_
     - `Batstone`_ et al., 2002; `Rosen and Jeppsson`_, 2006
   * - ADM1p (alias of ``ADM1_p_extension``)
     - `bsm2`_
     - `Alex`_ et al., 2008
   * - ASM1
     - `asm`_, `bsm1`_
     - `Henze`_ et al., 2006
   * - ASM2d
     - `asm`_, `bsm1`_
     - `Henze`_ et al., 2006
   * - mASM2d
     - `bsm2`_, `werf`_
     - `Alex`_ et al., 2008
   * - PM2
     - `pm2_batch`_, `pm2_ecorecover`_
     - --
   * - PM2ASM2d
     - --
     - --
   * - PM2ABACO2
     - --
     - --


Other kinetic modules
---------------------

These are standalone helpers (single :class:`qsdsan.Process` subclasses or rate-function utilities) that can be combined with the models above or used on their own.

.. list-table::
   :header-rows: 1
   :widths: 22 28 50

   * - Module
     - Implementation
     - Reference
   * - Aeration (``DiffusedAeration``)
     - `bsm1`_
     - `EPA design manual`_, 1989; `Mueller`_ et al., 2002
   * - Decay
     - `bwaise`_
     - `Trimmer`_ et al., 2020
   * - KineticReaction
     - --
     - --
   * - ASM_AeDigAddOn
     - --
     - --


Helper functions
----------------

Each model module also exposes a ``create_<modelname>_cmps`` factory that returns a :class:`qsdsan.CompiledComponents` set matching the model's state variables (e.g., :func:`~qsdsan.process_models.create_asm1_cmps`, :func:`~qsdsan.process_models.create_adm1_cmps`, :func:`~qsdsan.process_models.create_masm2d_cmps`). Several modules also expose rate-function utilities used internally by the compiled models — for example :func:`~qsdsan.process_models.pH_inhibit`, :func:`~qsdsan.process_models.Hill_inhibit`, and :func:`~qsdsan.process_models.rhos_adm1` from the ADM1 module.

The full inventory is available as ``qsdsan.process_models.__all__``.


.. Links
.. _EXPOsan: https://github.com/QSD-Group/EXPOsan
.. _adm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm
.. _asm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/asm
.. _bsm1: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1
.. _bsm2: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm2
.. _bwaise: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise
.. _pm2_batch: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_batch
.. _pm2_ecorecover: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_ecorecover
.. _werf: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/werf

.. _Alex: https://iwa-mia.org/wp-content/uploads/2022/09/TR3_BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf
.. _Batstone: https://iwaponline.com/ebooks/book/152/Anaerobic-Digestion-Model-No-1-ADM1
.. _EPA design manual: https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=3000464S.TXT
.. _Henze: https://iwaponline.com/ebooks/book/96/
.. _Mueller: https://www.taylorfrancis.com/chapters/mono/10.1201/9781420010343-5/diffused-aeration-james-mueller-william-boyle-ing-johannes-popel?context=ubx&refId=aa2e0148-4fb1-404b-a213-b46f92fc5c99
.. _Rosen and Jeppsson: https://iwa-mia.org/wp-content/uploads/2019/04/bmadm1_report061127.pdf
.. _Trimmer: https://doi.org/10.1021/acs.est.0c03296


Class references
----------------
.. toctree::
   :maxdepth: 1

   ADM1
   ADM1p
   ASM1
   ASM2d
   ASM_AeDigAddOn
   Aeration
   Decay
   KineticReaction
   PM2
   PM2ASM2d
   PM2ABACO2
