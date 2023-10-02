====================================================================================
QSDsan: Quantitative Sustainable Design for Sanitation and Resource Recovery Systems
====================================================================================

.. License
.. image:: https://img.shields.io/pypi/l/qsdsan?color=blue&logo=UIUC&style=flat
   :target: https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt

.. Tested Python version
.. image:: https://img.shields.io/pypi/pyversions/qsdsan?style=flat
   :target: https://pypi.python.org/pypi/qsdsan

.. PyPI version
.. image:: https://img.shields.io/pypi/v/qsdsan?style=flat&color=blue
   :target: https://pypi.org/project/qsdsan

.. DOI
.. image:: https://img.shields.io/badge/DOI-10.1039%2Fd2ew00455k-blue?style=flat
   :target: https://doi.org/10.1039/d2ew00455k

.. Documentation build
.. image:: https://readthedocs.org/projects/qsdsan/badge/?version=latest
   :target: https://qsdsan.readthedocs.io/en/latest

.. GitHub test and coverage of the main branch
.. image:: https://github.com/QSD-Group/QSDsan/actions/workflows/build-coverage.yml/badge.svg?branch=main
   :target: https://github.com/QSD-Group/QSDsan/actions/workflows/build-coverage.yml

.. Codecov
.. image:: https://codecov.io/gh/QSD-Group/QSDsan/branch/main/graph/badge.svg?token=Z1CASBXEOE
   :target: https://codecov.io/gh/QSD-Group/QSDsan

.. Binder launch of tutorials
.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/QSD-Group/QSDsan/main?filepath=%2Fdocs%2Fsource%2Ftutorials

.. Email subscription form
.. image:: https://img.shields.io/badge/news-subscribe-F3A93C?style=flat&logo=rss
   :target: https://groups.webservices.illinois.edu/subscribe/154591

.. Event calendar
.. image:: https://img.shields.io/badge/events-calendar-F3A93C?style=flat&logo=google%20calendar
   :target: https://qsdsan.readthedocs.io/en/latest/Events.html

.. YouTube video
.. image:: https://img.shields.io/endpoint?color=%23ff0000&label=YouTube%20 @qsd-group&url=https%3A%2F%2Fyoutube-channel-badge-blond.vercel.app%2Fapi%2Fvideos
   :target: https://www.youtube.com/@qsd-group

.. Code of Conduct
.. image:: https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg
   :target: https://qsdsan.readthedocs.io/en/latest/CODE_OF_CONDUCT.html

.. AppVeyor test of the stable branch, not in active use
..
   .. image:: https://img.shields.io/appveyor/build/yalinli2/QSDsan/main?label=build-stable&logo=appveyor
   :target: https://github.com/QSD-Group/QSDsan/tree/stable

|

.. contents::

|

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) of sanitation and resource recovery systems [1]_. It is one of a series of platforms that are being developed for the execution of QSD - a methodology for the research, design, and deployment of technologies and inform decision-making [2]_. It leverages the structure and modules developed in the `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_ platform [3]_ with additional functions tailored to sanitation processes.

As an open-source and impact-driven platform, QSDsan aims to identify configuration combinations, systematically probe interdependencies across technologies, and identify key sensitivities to contextual assumptions through the use of quantitative sustainable design methods (techno-economic analysis and life cycle assessment and under uncertainty). 

All systems developed with ``QSDsan`` are included in the package `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ - exposition of sanitation and resource recovery systems.

Additionally, another package, `DMsan <https://github.com/QSD-Group/DMsan>`_ (decision-making for sanitation and resource recovery systems), is being developed for decision-making among multiple dimensions of sustainability with consideration of location-specific contextual parameters.


Installation
------------
The easiest way is through ``pip``, in command-line interface (e.g., Anaconda prompt, terminal):

.. code::

    pip install qsdsan

If you need to upgrade:

.. code::

    pip install -U qsdsan

or for a specific version (replace X.X.X with the version number):

.. code::

    pip install qsdsan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/qsd-group/qsdsan>`_ (note that you can still use the ``-U`` flag for upgrading):

.. code::

    pip install git+https://github.com/QSD-Group/QSDsan.git


.. note::

   If this doesn't give you the newest ``qsdsan``, try ``pip uninstall qsdsan`` first.

   Also, you may need to update some ``qsdsan``'s dependency package (e.g., ' ``biosteam`` and ``thermosteam``) versions in order for the new ``qsdsan`` to run.


or other fork and/or branch (replace ``<USERNAME_OF_THE_FORK>`` and ``<BRANCH_NAME>`` with the desired fork and branch names)

.. code::

    pip install git+https://github.com/<USERNAME_OF_THE_FORK>/QSDsan.git@<BRANCH_NAME>


You can also download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

Note that development of this package is currently under initial stage with limited backward compatibility, please feel free to `submit an issue <https://github.com/QSD-Group/QSDsan/issues>`_ for any questions regarding package upgrading.

If you want to contribute to ``QSDsan``, please follow the steps in the `Contributing Guidelines <CONTRIBUTING#contributing-guidelines>`_ section of the documentation to clone the repository. If you find yourself struggle with the installation of QSDsan/setting up the environment, this extended version of `installation instructions <https://qsdsan.readthedocs.io/en/latest/tutorials/_installation.html>`_ might be helpful to you.


Documentation
-------------
You can find tutorials and documents at:

   https://qsdsan.readthedocs.io

All tutorials are written using Jupyter Notebook, you can run your own Jupyter environment, or you can click the ``launch binder`` badge on the top to launch the environment in your browser.

For each of these tutorials, we are also recording videos where one of the QSD group members will go through the tutorial step-by-step. We are gradually releasing these videos on our `YouTube channel <https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA>`_ so subscribe to receive updates!


About the Authors
-----------------
Please refer to `Contributors <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html#contributors>`_ section for a list of contributors.


Contributing
------------
Please refer to the `Contributing Guidelines <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html#contributing-guidelines>`_ section of the documentation for instructions and guidelines.


Stay Connected
--------------
If you would like to receive news related to the QSDsan platform, you can subscribe to email updates using `this form <https://groups.webservices.illinois.edu/subscribe/154591>`_ (don't worry, you will be able to unsubscribe :)). Thank you in advance for your interest!


QSDsan Events
-------------
We will keep this `calendar <https://calendar.google.com/calendar/embed?src=ep1au561lj8knfumpcd2a7ml08%40group.calendar.google.com&ctz=America%2FChicago>`_ up-to-date as we organize more events (office hours, workshops, etc.), click on the events in the calendar to see the details (including meeting links).


License Information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.


References
----------
.. [1] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289-2303. https://doi.org/10.1039/d2ew00455k.

.. [2] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. (Tutorial Review) Environ. Sci.: Water Res. Technol. 2022, 8 (11), 2439–2465. https://doi.org/10.1039/D2EW00431C.

.. [3] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J.S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.