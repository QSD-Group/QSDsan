====================================================================================
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
====================================================================================

.. image:: https://img.shields.io/pypi/l/qsdsan?color=blue&logo=UIUC&style=flat
   :target: https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
.. image:: https://img.shields.io/pypi/pyversions/qsdsan?style=flat
   :target: https://pypi.python.org/pypi/qsdsan
.. image:: https://img.shields.io/pypi/v/qsdsan?style=flat&color=blue
   :target: https://pypi.org/project/qsdsan/
.. image:: https://readthedocs.org/projects/qsdsan/badge/?version=latest
   :target: https://qsdsan.readthedocs.io/en/latest/
.. image:: https://img.shields.io/appveyor/build/yalinli2/QSDsan/main?label=build-main&logo=appveyor
   :target: https://github.com/QSD-Group/QSDsan/tree/main
.. image:: https://img.shields.io/appveyor/build/yalinli2/QSDsan/beta?label=build-beta&logo=appveyor
   :target: https://github.com/QSD-Group/QSDsan/tree/beta


.. contents::

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) of sanitation and resource recovery systems. It is one of a series of platforms that are being developed for the execution of QSD - a methodology for the research, design, and deployment of technologies and inform decision-making. [1]_ It leverages the structure and modules developed in the `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_ platform [2]_ with additional functions tailored to sanitation processes.

As an open-source and impact-driven platform, QSDsan aims to identify configuration combinations, systematically probe interdependencies across technologies, and identify key sensitivities to contextual assumptions through the use of quantitative sustainable design methods (techno-economic analysis and life cycle assessment and under uncertainty). 

All systems developed with ``QSDsan`` are included in the package `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ - exposition of sanitation and resource recovery systems.


Installation
------------
You can download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

If you use pip:

.. code:: bash

    pip install qsdsan

Note that development of this package is currently under initial stage with limited backward compatibility, please feel free to `submit an issue <https://github.com/QSD-Group/QSDsan/issues>`_ for any questions regarding package upgrading.

If you prefer the most recent version on GitHub, please follow the steps in the `Contributing to QSDsan <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation.


Documentation
-------------
You can find tutorials and documents at:

   https://qsdsan.readthedocs.io


About the authors
-----------------
Development and maintenance of the package is supported by the Quantitative Sustainable Design Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ at the `University of Illinois Urbana-Champaign (UIUC) <https://illinois.edu/>`_. Core contributors are listed below, please refer to the `author page <https://qsdsan.readthedocs.io/en/latest/AUTHORS.html>`_ for the full list of authors.

**Lead developers:**
   - `Yalin Li <zoe.yalin.li@gmail.com>`_ (current maintainer)
   - `Joy Cheung <joycheung1994@gmail.com>`_
   - `Stetson Rowles <lsr@illinois.edu>`_

**Project conception & funding support:**
   - `Jeremy Guest <jsguest@illinois.edu>`_

**Special acknowledgement:**
   - Yoel Cortés-Peña for helping many of the ``QSDsan`` members get started on Python and package development.


Contributing
------------
Please refer to the `Contributing to QSDsan <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation for instructions and guidelines.


License information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.

References
----------
.. [1] Li, Y.; Hand, S.; Trimmer, J. T.; Byrne, D. M.; Chambers, K. G.; Lohman, H. A. C.; Shi, R.; Zhang, X.; Cook, S. M.; Guest, J. S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. In Prep. 2021.

.. [2] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.