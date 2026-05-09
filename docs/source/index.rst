QSDsan: Quantitative Sustainable Design for Sanitation and Resource Recovery Systems
====================================================================================

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) [1]_ of sanitation and resource recovery systems [2]_. Built in Python, it integrates process modeling, system simulation, techno-economic analysis (TEA), and life cycle assessment (LCA) to support transparent, reproducible, and comprehensive evaluation of emerging technologies. The platform provides modular and extensible tools to compare treatment configurations, evaluate resource recovery opportunities, and assess energy, cost, and sustainability tradeoffs. The platform has a range of built-in :ref:`unit operations <unit_operations>` and :ref:`systems <systems>`, including commonly used treatment processes and biokinetic models. These capabilities support the research, development, and deployment (RD&D) of early-stage technologies for advancing sustainable water management and circular resource recovery.

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
      :link: systems
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
``QSDsan`` requires Python 3.12 or newer. The easiest way to install ``QSDsan`` is through ``pip`` in a command-line interface (e.g., Anaconda Prompt, terminal):

.. code::

    pip install qsdsan

To upgrade an existing installation:

.. code::

    pip install -U qsdsan

To install a specific version, replace ``X.X.X`` with the version number:

.. code::

    pip install qsdsan==X.X.X

To install the latest GitHub version from the `main branch <https://github.com/QSD-Group/QSDsan>`_:

.. code::

    pip install git+https://github.com/QSD-Group/QSDsan.git

To install from another fork and/or branch, replace ``<USERNAME_OF_THE_FORK>`` and ``<BRANCH_NAME>``:

.. code::

    pip install git+https://github.com/<USERNAME_OF_THE_FORK>/QSDsan.git@<BRANCH_NAME>

You can also download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

For diagram generation, ``QSDsan`` uses Graphviz. If diagrams fail to render, install Graphviz following the `official Graphviz download instructions <https://graphviz.org/download/>`_ and see :ref:`graphviz-installation` for a quick check.


Join the Community
------------------
We would like to build an open and welcoming community, you can always post issues on our `GitHub homepage <https://github.com/QSD-Group/QSDsan/issues>`_ or contact any of the Quantitative Sustainable Design Group members. We are always excited to have new members in our team.

If you would like to contribute, please follow our `Contributing Guidelines`_ and the `Code of Conduct <CODE_OF_CONDUCT.html>`_, thank you for making ``QSDsan`` better!


``QSDsan`` is and will stay open source under University of Illinois/NCSA Open Source License. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). Please refer to the `license <https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt>`_ page for details.

.. If you would like to receive news related to the QSDsan platform, you can subscribe to email updates using `this form <https://groups.webservices.illinois.edu/subscribe/154591>`_ (don't worry, you will be able to unsubscribe :)). Thank you in advance for your interest!


.. QSDsan Events
.. -------------

.. We will keep the calendar up-to-date as we organize more events (office hours, workshops, etc.), click on the events in the calendar to see the details (including meeting links).

.. .. raw:: html
..
..     <embed>
..         <iframe src="https://calendar.google.com/calendar/embed?src=ep1au561lj8knfumpcd2a7ml08%40group.calendar.google.com&ctz=America%2FChicago" style="border: 0" width="100%" height="600" scrolling="no"></iframe>
..     </embed>


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

   Systems


.. toctree::
   :maxdepth: 2
   :hidden:

   FAQ


.. toctree::
   :maxdepth: 1
   :hidden:

   CONTRIBUTING


.. toctree::
   :maxdepth: 1
   :hidden:

   CODE_OF_CONDUCT


.. .. toctree::
..    :maxdepth: 1
..    :hidden:
..
..    calendar


.. toctree::
   :hidden:
   :maxdepth: 1

   CHANGELOG


References
----------
.. [1] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. (Tutorial Review) Environ. Sci.: Water Res. Technol. 2022, 8 (11), 2439–2465. https://doi.org/10.1039/D2EW00431C.

.. [2] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289-2303. https://doi.org/10.1039/d2ew00455k.

.. Links
.. _Contributing Guidelines: CONTRIBUTING.html#contributing-guidelines
