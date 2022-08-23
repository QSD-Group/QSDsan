.. _developed_systems:

Developed Systems
=================

This page documents the systems that have been/is being developed using ``QSDsan`` with links to the source codes in GitHub and publications.

+--------------------------+--------------------+--------------------------------+-------------------+
| System                   | Source Codes       | Publication                    | Module Status     |
+==========================+====================+================================+===================+
| Biogenic Refinery        | - Archived `br1`_  | - `Rowles`_ et al., 2022       | Completed         |
|                          | - Current  `br2`_  |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| BSM1 (benchmark          | - Archived `bsm1`_ | - `Alex`_ et al., 2008         | Completed         |
| simulation model no. 1)  | - Current  `bsm2`_ | - `Li and Zhang`_ et al., 2022 |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Bwaise                   | - Archived `bw1`_  | - `Trimmer`_ et al., 2020      | Completed         |
|                          | - Current  `bw2`_  | - `Li and Zhang`_ et al., 2022 |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| CAS (conventional        | - Current  `cas`_  | - `Shoener`_ et al., 2016      | Completed         |
| activated sludge)        |                    |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Eco-San                  | - Current   `es`_  | NA                             | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+
| NEWgenerator             | - Archived `ng1`_  | - Watabe et al., *In Prep.*    | Completed         |
| (under NDA)              | - Current  `ng2`_  |                                |                   |
+--------------------------+--------------------+--------------------------------+-------------------+
| Reclaimer                | - Current   `re`_  | - `Trotochaud`_ et al., 2020   | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+
| SCG Zyclonic (under NDA) | - Current   `sz`_  | NA                             | Completed         |
+--------------------------+--------------------+--------------------------------+-------------------+

Notes:
    - "Under NDA" indicates that the system is under non-disclosure agreement with the technology design team and unfortunately we are not able to share the codes in full at this stage. The source link code will lead to a private repository that only individuals who have signed the NDA can access.
    - "Archived" is the version used when the linked literature is published.
    - "Current" is the version that has been updated to be compatible with the most up-to-date version of ``QSDsan``.


.. Links
.. _br1: https://github.com/QSD-Group/EXPOsan/releases/tag/archive%2FBR_OmniProcessor
.. _br2: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biogenic_refinery
.. _Rowles: https://doi.org/10.1021/acsenvironau.2c00022

.. _bsm1: https://pypi.org/project/exposan/1.1.4/
.. _bsm2: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1
.. _Alex: http://iwa-mia.org/wp-content/uploads/2019/04/BSM_TG_Tech_Report_no_1_BSM1_General_Description.pdf
.. _Li and Zhang: https://arxiv.org/abs/2203.06243

.. _bw1: https://pypi.org/project/exposan/1.1.4/
.. _bw2: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise
.. _Trimmer: https://doi.org/10.1021/acs.est.0c03296

.. _cas: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/cas
.. _Shoener: https://pubs.rsc.org/en/content/articlelanding/2016/ee/c5ee03715h

.. _es: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/eco_san

.. _ng1: https://github.com/QSD-Group/EXPOsan-private/tree/newgen/exposan/newgen
.. _ng2: https://github.com/QSD-Group/EXPOsan-private/tree/main/exposan/new_generator

.. _re: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/reclaimer
.. _Trotochaud: https://doi.org/10.1021/acs.est.0c02755

.. _sz: https://github.com/QSD-Group/EXPOsan-private/tree/main/exposan/scg_zyclonic
