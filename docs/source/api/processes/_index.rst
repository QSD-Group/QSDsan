processes
=========

This module contains the existing process models within ``QSDsan`` for dynamic simulation, the following table provides the reference and example implementation.

List of Process Models
----------------------
+----------------+------------------+-----------------------------+
| Process Models | Implementation   | Reference                   |
+================+==================+=============================+
| ADM1           | `adm`_           | `Batstone`_ et al., 2002    |
|                |                  | `Rosen and Jeppsson`_, 2006 |
+----------------+------------------+-----------------------------+
| Aeration       | `bsm1`_          | `EPA design manual`_, 1989  |
|                |                  | `Mueller`_ et al., 2002     |
+----------------+------------------+-----------------------------+
| ASM1           | `asm`_ & `bsm1`_ | `Henze`_ et al., 2006       |
+----------------+------------------+-----------------------------+
| ASM2d          | `asm`_ & `bsm1`_ | `Henze`_ et al., 2006       |
+----------------+------------------+-----------------------------+

.. Links
.. _adm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm
.. _asm: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/asm
.. _bsm1: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1

.. _Batstone: https://iwaponline.com/ebooks/book/152/Anaerobic-Digestion-Model-No-1-ADM1
.. _EPA design manual: https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=3000464S.TXT
.. _Henze: https://iwaponline.com/ebooks/book/96/
.. _Mueller: https://www.taylorfrancis.com/chapters/mono/10.1201/9781420010343-5/diffused-aeration-james-mueller-william-boyle-ing-johannes-popel?context=ubx&refId=aa2e0148-4fb1-404b-a213-b46f92fc5c99
.. _Rosen and Jeppsson: http://iwa-mia.org/wp-content/uploads/2019/04/bmadm1_report061127.pdf


Links to docs
-------------
.. toctree::
   :maxdepth: 1

   ADM1
   Aeration
   ASM1
   ASM2d