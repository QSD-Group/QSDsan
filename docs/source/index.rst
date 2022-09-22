QSDsan: Quantitative Sustainable Design for Sanitation and Resource Recovery Systems
====================================================================================

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) [1]_ of sanitation and resource recovery systems [2]_. It leverages existing platforms such as `BioSTEAM <https://biosteam.readthedocs.io>`_ [3]_  with enhanced features tailored to sanitation an resource recovery technologies. Through the integration with `DMsan <https://github.com/QSD-Group/DMsan>`_ (decision-making for sanitation and resource recovery systems), this platform can be used to guide the research, development, and deployment (RD&D) of early-stage technologies considering location-specific parameters and stakeholder priorities.


.. grid:: 2 4 4 4

   .. grid-item-card::  Tutorials
      :text-align: center
      :link: tutorials
      :link-type: ref

      .. figure:: images/tutorials_icon.svg
         :height: 100
         :align: center
          
   .. grid-item-card::  API
      :text-align: center
      :link: api
      :link-type: ref

      .. figure:: images/api_icon.svg
         :height: 100
         :align: center

   .. grid-item-card::  Systems
      :text-align: center
      :link: developed_systems
      :link-type: ref

      .. figure:: images/systems_icon.svg
         :height: 100
         :align: center
      
   .. grid-item-card::  FAQ
      :text-align: center
      :link: faq
      :link-type: ref

      .. figure:: images/faq_icon.svg
         :height: 100
         :align: center


Installation
------------
The easiest way is through ``pip``, in command-line interface (Anaconda prompt, terminal):

.. code::

    pip install qsdsan

If you need to update:

.. code::

    pip install -U qsdsan

Or for a specific version (replace X.X.X with the version number):

.. code::

    pip install qsdsan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/qsd-group/qsdsan>`_:

.. code::

    pip install git+https://github.com/QSD-Group/QSDsan.git


You can also download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

Note that development of this package is currently under initial stage with limited backward compatibility, please feel free to `submit an issue <https://github.com/QSD-Group/QSDsan/issues>`_ for any questions regarding package upgrading.

If you are a developer and want to contribute to ``QSDsan``, please follow the steps in the `contributing`_ section of the documentation to clone the repository.


Join the Community
------------------
We would like to build an open and welcoming community, you can always post issues on our `GitHub homepage <https://github.com/QSD-Group/QSDsan/issues>`_ or contact any of the Quantitative Sustainable Design Group members. We are always excited to have new members in our team.

If you would like to contribute, please follow our `contributing`_ guidelines and the `code of conduct <https://qsdsan.readthedocs.io/en/latest/CODE_OF_CONDUCT.html>`_ (a bonus if you use our `templates <https://github.com/QSD-Group/QSDsan/tree/main/docs/source/community/templates>`_), thank you for making ``QSDsan`` better!


``QSDsan`` is and will stay open source under University of Illinois/NCSA Open Source License. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). Please refer to the `license <https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt>`_ page for details.

If you would like to receive news related to the QSDsan platform, you can subscribe to email updates using `this form <https://groups.webservices.illinois.edu/subscribe/154591>`_ (don't worry, you will be able to unsubscribe :)). Thank you in advance for your interest!


QSDsan Events
-------------

We will keep the calendar up-to-date as we organize more events (office hours, workshops, etc.), click on the events in the calendar to see the details (including meeting links).

.. raw:: html

    <embed>
        <iframe src="https://calendar.google.com/calendar/embed?src=ep1au561lj8knfumpcd2a7ml08%40group.calendar.google.com&ctz=America%2FChicago" style="border: 0" width="100%" height="600" scrolling="no"></iframe>
    </embed>


.. Hidden TOCs for navigation bar
.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials/_index


.. toctree::
   :maxdepth: 1
   :hidden:

   api/_index


.. toctree::
   :maxdepth: 1
   :hidden:

   Developed_Systems


.. toctree::
   :maxdepth: 2
   :hidden:

   FAQ


.. toctree::
   :maxdepth: 1
   :hidden:

   core_developers/_index


.. toctree::
   :maxdepth: 1
   :hidden:

   calendar


.. toctree::
   :hidden:
   :maxdepth: 1

   CHANGELOG


References
----------
.. [1] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. Available on `ChemRxiv <https://chemrxiv.org/engage/chemrxiv/article-details/629df71e97e76a377cc7f06e>`_.

.. [2] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. Accepted, 2022. https://doi.org/10.1039/d2ew00455k.

.. [3] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J.S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.


.. Links
.. _contributing: https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html