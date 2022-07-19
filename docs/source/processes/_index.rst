processes
=========

This module contains the existing process models within ``QSDsan`` for dynamic simulation, the following table provides the reference and example implementation.

List of Process Models
----------------------
+----------------+------------------+-----------------------+
| Process Models | Implementation   | Reference             |
+================+==================+=======================+
| ADM1           | `adm`_           | `Henze`_ et al., 2006 |
+----------------+------------------+-----------------------+
| Aeration       | `bsm1`_          | `Alex`_ et al., 2008  |
+----------------+------------------+-----------------------+
| ASM1           | `asm`_ & `bsm1`_ | `Henze`_ et al., 2006 |
+----------------+------------------+-----------------------+
| ASM2d          | `asm`_ & `bsm1`_ | `Henze`_ et al., 2006 |
+----------------+------------------+-----------------------+

.. Links
.. _adm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm
.. _Alex: http://iwa-mia.org/wp-content/uploads/2019/04/BSM_TG_Tech_Report_no_1_BSM1_General_Description.pdf
.. _asm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/asm
.. _bsm1: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1
.. _Henze: https://iwaponline.com/ebooks/book/96/


Links to docs
-------------
.. toctree::
   :maxdepth: 1

   ADM1
   Aeration
   ASM1
   ASM2d