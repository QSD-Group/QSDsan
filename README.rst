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
.. image:: docs/source/images/custom_binder_logo.svg
   .. :target: https://mybinder.org/v2/gh/QSD-Group/QSDsan/main?filepath=%2Fdocs%2Fsource%2Ftutorials
   :target: https://mybinder.org/v2/gh/QSD-Group/QSDsan-env/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FQSD-group%252FQSDsan%26urlpath%3Dtree%252FQSDsan%252Fdocs%252Fsource%252Ftutorials%26branch%3Dmain

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

If you want to contribute to ``QSDsan``, please follow the steps in the `Contributing Guidelines <CONTRIBUTING.rst#contributing-guidelines>`_ section of the documentation to clone the repository. If you find yourself struggle with the installation of QSDsan/setting up the environment, this extended version of `installation instructions <docs/source/tutorials/_installation.rst>`_ might be helpful to you.


Documentation
-------------
You can find tutorials and documents at:

   https://qsdsan.readthedocs.io

All tutorials are written using Jupyter Notebook, you can run your own Jupyter environment, or you can click the ``launch binder`` badge on the top to launch the environment in your browser.

For each of these tutorials, we are also recording videos where one of the QSD group members will go through the tutorial step-by-step. We are gradually releasing these videos on our `YouTube channel <https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA>`_ so subscribe to receive updates!


About the Authors
-----------------
Please refer to `Contributors <CONTRIBUTING.rst#contributors>`_ section for a list of contributors.


Contributing
------------
Please refer to the `Contributing Guidelines <CONTRIBUTING.rst#contributing-guidelines>`_ section of the documentation for instructions and guidelines.


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


.. Custom launch badges: https://mybinder.readthedocs.io/en/latest/howto/badges.html
.. https://img.shields.io/badge/launch-binder%20%7C%20tutorial-579ACA.svg?logo=data:image/png;base64,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