.. _developed_systems:

Developed Systems
=================

This page documents the systems that have been/is being developed using ``QSDsan`` with links to the source codes in GitHub and publications.

+--------------------------+--------------------+--------------------------------+-------------------+
| System                   | Source Codes       | Publication                    | Module Status     |
+==========================+====================+================================+===================+
| Biogenic Refinery        | - `br_archived`_   | - `Rowles`_ et al., 2022       | Completed         |
|                          | - `br_current`_    |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| BSM1 (benchmark          | - `bsm1_archived`_ | - `Alex`_ et al., 2008         | Completed         |
| simulation model no. 1)  | - `bsm1_current`_  | - `Li and Zhang`_ et al., 2022 |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Bwaise                   | - `bw_archived`_   | - `Trimmer`_ et al., 2020      | Completed         |
|                          | - `bw_current`_    | - `Li and Zhang`_ et al., 2022 |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| CAS (conventional        | - `cas`_           | - `Shoener`_ et al., 2016      | Completed         |
| activated sludge)        |                    |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Eco-San                  | - `es`_            | NA                             | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+
| NEWgenerator             | - `ng_archived`_   | - Watabe et al., *In Prep.*    | Completed         |
| (under NDA)              | - `ng_current`_    |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Reclaimer                | - `re`_            | - `Trotochaud`_ et al., 2020   | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+
| SCG Zyclonic (under NDA) | - `sz`_            | NA                             | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+

Notes:
    - "Under NDA" indicates that the system is under non-disclosure agreement with the technology design team and unfortunately we are not able to share the codes in full at this stage. The source link code will lead to a private repository that only individuals who have signed the NDA can access.
    - "Archived" is the version used when the linked literature is published.
    - "Current" is the version that has been updated to be compatible with the most up-to-date version of ``QSDsan``.


.. Links
.. _br_archived: https://github.com/QSD-Group/EXPOsan/releases/tag/archive%2FBR_OmniProcessor
.. _br_current: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biogenic_refinery
.. _Rowles: https://doi.org/10.1021/acsenvironau.2c00022

.. _bsm1_archived: https://pypi.org/project/exposan/1.1.4
.. _bsm1_current: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1
.. _Alex: http://iwa-mia.org/wp-content/uploads/2019/04/BSM_TG_Tech_Report_no_1_BSM1_General_Description.pdf
.. _Li and Zhang: https://doi.org/10.1039/d2ew00455k

.. _bw_archived: https://pypi.org/project/exposan/1.1.4
.. _bw_current: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise
.. _Trimmer: https://doi.org/10.1021/acs.est.0c03296

.. _cas: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/cas
.. _Shoener: https://pubs.rsc.org/en/content/articlelanding/2016/ee/c5ee03715h

.. _es: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/eco_san

.. _ng_current: https://github.com/QSD-Group/EXPOsan-private/tree/newgen/exposan/newgen
.. _ng_archived: https://github.com/QSD-Group/EXPOsan-private/tree/main/exposan/new_generator

.. _re: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/reclaimer
.. _Trotochaud: https://doi.org/10.1021/acs.est.0c02755

.. _sz: https://github.com/QSD-Group/EXPOsan-private/tree/main/exposan/scg_zyclonic
