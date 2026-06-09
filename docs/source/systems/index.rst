.. _systems:

Systems & Publications
======================

``QSDsan`` is built upon the quantitative sustainable design (QSD) methodology, which provides a structured approach to prioritize the research, development, and deployment of early-stage technologies. Leveraging the QSD methodology, ``QSDsan`` powers a growing library of systems for wastewater treatment, sanitation, and resource recovery technologies. This page introduces the foundation metholodologies and includes a non-exhaustive list of the published systems with links to their source code in ``EXPOsan`` and their publications.

Platform and methodology
-------------------------

.. grid:: 1
   :gutter: 4

   .. grid-item-card:: Structured methodology: (QSD)

      .. image:: ../images/systems/qsd.png
         :width: 400px
         :align: left

      Research, development, and deployment (RD&D) of innovative technologies are often impeded by the lack of transparent, systematic, and agile approaches to prioritize investment across the expansive landscape of technologies and design/operational decisions. This tutorial review synthesizes research on sustainability analyses to present Quantitative Sustainable Design (QSD) – a structured methodology to expedite the RD&D of water, sanitation, and resource recovery technologies.

      *Li, Y.; Trimmer, J. T.; Hand, S.; Zhang, X.; Chambers, K. G.; Lohman, H. A. C.;
      Shi, R.; Byrne, D. M.; Cook, S. M.; Guest, J. S. Quantitative Sustainable Design
      (QSD): A Methodology for the Prioritization of Research, Development, and Deployment
      of Technologies. Environ. Sci.: Water Res. Technol. 2022, 8 (11), 2439–2465.*

      +++

      .. button-link:: https://doi.org/10.1039/D2EW00431C
         :color: primary

         Read Paper

   .. grid-item-card:: Integrated platform: QSDsan

      .. image:: ../images/systems/qsdsan.png
         :width: 400px
         :align: left

      An open-source Python tool that integrates system design, simulation, and sustainability characterization (techno-economic analysis and life cycle assessment) to quickly identify critical barriers, prioritize research opportunities, and navigate multi-dimensional sustainability tradeoffs for technology RD&D.

      *Li, Y.; Zhang, X.; Morgan, V. L.; Lohman, H. A. C.; Rowles, L. S.; Mittal, S.;
      Kogler, A.; Cusick, R. D.; Tarpeh, W. A.; Guest, J. S. QSDsan: An Integrated
      Platform for Quantitative Sustainable Design of Sanitation and Resource Recovery
      Systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289–2303.*

      +++

      .. button-link:: https://doi.org/10.1039/D2EW00455K
         :color: primary

         Read Paper

   .. grid-item-card:: Decision-making: DMsan

      .. image:: ../images/systems/dmsan.png
         :width: 400px
         :align: left

      A multi-criteria decision analysis package that integrates with ``QSDsan`` to
      compare alternatives across technical, resource-recovery, economic, environmental,
      and social criteria.

      *Lohman, H. A. C.; Morgan, V. L.; Li, Y.; Zhang, X.; Rowles, L. S.; Cook, S. M.;
      Guest, J. S. DMsan: A Multi-Criteria Decision Analysis Framework and Package to
      Characterize Contextualized Sustainability of Sanitation and Resource Recovery
      Technologies. ACS Environ. Au 2023, 3 (3), 179–192.*

      +++

      .. button-link:: https://doi.org/10.1021/acsenvironau.2c00067
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/DMsan/
         :color: secondary

         Source Code


Water Resource Recovery Facilities
----------------------------------

Benchmark Simulation Models
***************************
The Modelling and Integrated Assessment (MIA) Specialist Group of the International Water Association has established benchmark simulation models (BSMs) to provide a consistent environment for wastewater treatment plant (WWTP)/water resource recovery facility (WRRF) evaluation (see `BSM webpage <https://iwa-mia.org/benchmarking>`_ and `MATLAB implementation and report <https://github.com/wwtmodels/Benchmark-Simulation-Models>`_).

When publishing the paper that introduces QSDsan, we validated the process modeling and dynamic simulation capacities of QSDsan through BSM1. BSM2 is also implemented in ``EXPOsan``.

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: Benchmark Simulation Model No. 1 (BSM1)

      Validated Python implementation of the `BSM1 system <https://iwa-mia.org/benchmarking/#BSM1>`_ by the International Water Association (IWA).

      +++

      .. button-link:: https://doi.org/10.1039/D2EW00455K
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1
         :color: secondary

         Source Code

   .. grid-item-card:: Benchmark Simulation Model No. 2 (BSM2)

      Validated Python implementation of the IWA `BSM2 system <https://iwa-mia.org/benchmarking/#BSM2>`_.

      +++

      .. button-link:: https://doi.org/10.1039/D2EW00455K
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm2
         :color: secondary

         Source Code


WERF Treatment Trains
*********************
In Zhang et al., 2026, we developed 18 benchmark combinations of liquid and solid treatment trains, which cover over 70% of the total treatment capacity of publicly owned treatment works (POTWs) in the Contiguous United States. These configurations were based on the Water Environment Research Foundation (WERF, now a part of the Water Research Foundation, WRF), report on net-zero energy solutions for WRRFs.

These simulation models have been implemented in ``EXPOsan``. More details can be found in :ref:`the interactive page <wrrf_interactive>`.

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: WERF Treatment Trains

      *Zhang, X.; Rai, S.; Wang, Z.; Li, Y.; Guest, J. S. An Agile Benchmarking Framework for Wastewater Resource Recovery Technologies. npj Clean Water 2025, 9 (1), 4.*

      +++

      .. button-ref:: wrrf_interactive
         :ref-type: ref
         :color: primary

         Interactive Page

      .. button-link:: https://doi.org/10.1038/s41545-025-00537-4
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/werf
         :color: secondary

         Source Code


Other WRRFs
***********

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: Anaerobic Digestion Model No. 1 (ADM1)

      Validated Python implementation of IWA `ADM1 <https://doi.org/10.2166/9781780403052>`_.

      +++

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm
         :color: secondary

         Source Code

   .. grid-item-card:: Activated Sludge Models (ASMs)

      Validated Python implementation of IWA `ASM <https://doi.org/10.2166/9781780402369>`_.

      +++

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/asm
         :color: secondary

         Source Code

   .. grid-item-card:: Conventional Activated Sludge

      Python implementation of the CAS system as described in `Shoerner et al. <https://doi.org/10.1039/C5EE03715H>`_

      +++

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/cas
         :color: secondary

         Source Code

      .. button-link:: https://github.com/QSD-Group/AnMBR
         :color: secondary

         Source Code | MATLAB

   .. grid-item-card::  Modular Encapsulated Two-stage Anaerobic Biological system (METAB)

      *Zhang, X.; Arnold, W. A.; Wright, N.; Novak, P. J.; Guest, J. S. Prioritization of Early-Stage Research and Development of a Hydrogel-Encapsulated Anaerobic Technology for Distributed Treatment of High Strength Organic Wastewater. Environ. Sci. Technol. 2024, 58 (44), 19651–19665.*

      +++

      .. button-link:: https://doi.org/10.1021/acs.est.4c05389
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/metab
         :color: secondary

         Source Code

   .. grid-item-card:: EcoRecover

      *Kim, G.-Y.; Molitor, H. R.; Zhang, X.; Li, Y.; Shoener, B. D.; Schramm, S. M.;
      Morgenroth, E.; Snowling, S. D.; Hartnett, E.; Bradley, I. M.; Pinto, A. J.; Guest,
      J. S. Development of an Open-Source Process Simulator for Microalgae-Based Tertiary
      Phosphorus Recovery. npj Clean Water 2025, 9 (1), 13.*

      +++

      .. button-link:: https://doi.org/10.1038/s41545-025-00545-4
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_batch
         :color: secondary

         Source Code | Batch

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_ecorecover
         :color: secondary

         Source Code | EcoRecover


----------

Non-sewered sanitation systems (NSSSs)
--------------------------------------

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: NSS Typologies across 77 Countries

      *Lohman, H. A. C.; Li, Y.; Zhang, X.; Morgan, V. L.; Watabe, S.; Rowles, L. S.; Cusick, R. D.; Guest, J. S. Defining Economic and Environmental Typologies across 77 Countries to Prioritize Opportunities for Nonsewered Sanitation. Environ. Sci. Technol. 2025, 59 (29), 15101–15114.*

      +++

      .. button-link:: https://doi.org/10.1021/acs.est.5c02064
         :color: primary

         Read Paper

   .. grid-item-card:: Biogenic Refinery

      *Rowles, L. S.; Morgan, V. L.; Li, Y.; Zhang, X.; Watabe, S.; Stephen, T.; Lohman,
      H. A. C.; DeSouza, D.; Hallowell, J.; Cusick, R. D.; Guest, J. S. Financial
      Viability and Environmental Sustainability of Fecal Sludge Treatment with Pyrolysis
      Omni Processors. ACS Environ. Au 2022, 2 (5), 455–466.*

      +++

      .. button-link:: https://doi.org/10.1021/acsenvironau.2c00022
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biogenic_refinery
         :color: secondary

         Source Code


   .. grid-item-card:: Sanitation Alternatives in Bwaise, Uganda

      *Trimmer, J. T.; Lohman, H. A. C.; Byrne, D. M.; Houser, S. A.; Jjuuko, F.; Katende,
      D.; Banadda, N.; Zerai, A.; Miller, D. C.; Guest, J. S. Navigating Multidimensional
      Social–Ecological System Trade-Offs across Sanitation Alternatives in an Urban
      Informal Settlement. Environ. Sci. Technol. 2020, 54 (19), 12641–12653.*

      +++

      .. button-link:: https://doi.org/10.1021/acs.est.0c03296
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise
         :color: secondary

         Source Code

      .. button-link:: https://github.com/QSD-Group/Bwaise-sanitation-alternatives
         :color: secondary

         Source Code | Original


   .. grid-item-card:: Eco-San

      Based on the Eco-San system developed by Yixing Eco-sanitary Manufacture Co., Ltd.

      +++

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/eco_san
         :color: secondary

         Source Code

   .. grid-item-card:: NEWgenerator

      *Watabe, S.; Lohman, H. A. C.; Li, Y.; Morgan, V. L.; Rowles, L. S.; Stephen, T.;
      Shyu, H.-Y.; Bair, R. A.; Castro, C. J.; Cusick, R. D.; Yeh, D. H.; Guest, J. S.
      Advancing the Economic and Environmental Sustainability of the NEWgenerator
      Nonsewered Sanitation System. ACS Environ. Au 2023, 3 (4), 209–222.*

      **Note:** the NEWgenerator system is under non-disclosure agreement (NDA), thus unit design is not publicly available, but the system design is implemented in ``EXPOsan``.

      +++

      .. button-link:: https://doi.org/10.1021/acsenvironau.3c00001
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/new_generator
         :color: secondary

         Source Code

   .. grid-item-card:: Pasteurization Mechanical Dewatering (PMD) & Supercritical Water Oxidation (SCWO)

      *Wang, Z.; Feng, J.; Shi, B.; Mendoza, J. A.; Zhang, X.; Trousdale, N.; Cusick, R. D.; Yee, S.; Guest, J. S. The Potential of Thermomechanical and Thermochemical Processes to Enable Sustainable Household Sanitation. Environ. Sci. Technol. 2026, 60 (8), 6227–6238.*

      +++

      .. button-link:: https://doi.org/10.1021/acs.est.5c15639
         :color: primary

         Read Paper


   .. grid-item-card:: Reclaimer

      Based on the work described in `Trotochaud et al. <https://doi.org/10.1021/acs.est.0c02755>`_ for the Reclaimer 2.0 system designed by researchers at Duke University.

      +++

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/reclaimer
         :color: secondary

         Source Code


----------

Other Systems
-------------

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: Hydrothermal systems for biobinder and biofuels from food waste

      *Ahmad, A.; Kawale, H. D.; Summers, S.; Bogarin Cantero, B. C.; Allen, C. M.; Hajj,
      R. M.; Davidson, P. C.; Zhang, Y.; Li, Y. Financial Viability and Carbon Intensity
      of Hydrothermal Waste Valorization Systems for Bio-Based Asphalt Binder. Chemical
      Engineering Journal 2026, 528, 172283.*

      +++

      .. button-link:: https://doi.org/10.1016/j.cej.2025.172283
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biobinder
         :color: secondary

         Source Code

   .. grid-item-card:: Hydroxyapatite (HAp) Synthesis from Urine

      *Müller, I. E.; Lin, A. Y. W.; Otani, Y.; Zhang, X.; Wu, Z.-Y.; Kisailus, D.;
      Mouncey, N. J.; Guest, J. S.; Rad, B.; Ercius, P.; Yoshikuni, Y. Cost-Effective
      Urine Recycling Enabled by a Synthetic Osteoyeast Platform for Production of
      Hydroxyapatite. Nat Commun 2025, 16 (1), 4216.*

      +++

      .. button-link:: https://doi.org/10.1038/s41467-025-59416-8
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/hap
         :color: secondary

         Source Code

   .. grid-item-card:: Hydrothermal systems for fuel and fertilizer production from wet organic wastes

      *Feng, J.; Strathmann, T. J.; Guest, J. S. Hydrothermal-Based Wastewater Solids
      Management for Targeted Resource Recovery and Decarbonization in the Contiguous U.S.
      Environ. Sci. Technol. 2025.*

      +++

      .. button-link:: https://doi.org/10.1021/acs.est.5c07190
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/htl
         :color: secondary

         Source Code

   .. grid-item-card:: Point-of-use disinfection technologies

      *Elijah, B. C.; Ahmad, A.; Li, Y.; Plazas-Tuttle, J.; Rowles, L. S. Assessing the
      Relative Sustainability of Point-of-Use Water Disinfection Technologies for Off-Grid
      Communities. ACS Environ. Au 2024, 4 (5), 248–259.*

      +++

      .. button-link:: https://doi.org/10.1021/acsenvironau.4c00017
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pou_disinfection
         :color: secondary

         Source Code

   .. grid-item-card:: Hydrothermal systems for sustainable aviation fuel from food waste

      *Si, B.; Wang, Z.; Watson, J.; Summers, S.; Li, Y.; Yu, S.; Yang, H.; Yang, Z.;
      Heyne, J. S.; Jiang, J.; Ren, Z. J.; Ma, H.; Wang, C.; Wang, P.; Zhang, Y. A
      Circular Hydrothermal Refinery for Sustainable Aviation Fuel from Food Waste. Nat
      Sustain 2026, 1–11.*

      +++

      .. button-link:: https://doi.org/10.1038/s41893-026-01848-1
         :color: primary

         Read Paper

      .. button-link:: https://github.com/QSD-Group/EXPOsan/tree/main/exposan/saf
         :color: secondary

         Source Code


----------

Additional Publications
-----------------------

Papers published by external users.

#. `Advancing the Economic and Environmental Sustainability of Rare Earth Element Recovery from Phosphogypsum <https://doi.org/10.1021/acs.est.5c04952>`_
#. `N2O as reactant rather than pollutant at wastewater treatment plants: Life Cycle Assessment and Techno-Economic Analysis of N2O-to-phenol <https://chemrxiv.org/doi/10.26434/chemrxiv-2026-m4vwz>`_
#. `Cost and Carbon Implications of Industrial Organic Load Reduction across Water Resource Recovery Facility Typologies <https://doi.org/10.1021/acs.est.5c11447>`_


.. toctree::
   :hidden:

   wrrf_interactive
