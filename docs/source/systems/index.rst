.. _systems:

Published Systems
=================

This page documents the systems that have been/is being developed using ``QSDsan`` with links to the source codes in GitHub and publications.

Water Resource Recovery Facilities
-----------------------------------

Benchmark Simulation Models
***************************
The Modelling and Integrated Assessment (MIA) Specialist Group of the International Water Association has established benchmark simulation models (BSMs) to provide a consistent environment for wastewater treatment plant (WWTP)/water resource recovery facility (WRRF) evaluation (see `BSM webpage <https://iwa-mia.org/benchmarking>`_ and `MATLAB implementation and report <https://github.com/wwtmodels/Benchmark-Simulation-Models>`_).

When publishing the paper that introduces QSDsan [1]_, we validated the process modeling and dynamic simulation capacities of QSDsan through BSM1 (`bsm1 EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1>`_, `bsm1 archived codes <https://pypi.org/project/exposan/1.1.4>`_). BSM2 is also implemented in `EXPOsan <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm2>`_.


WERF Treatment Trains
*********************
In Zhang et al., 2026 [2]_, we developed 18 benchmark combinations of liquid and solid treatment trains, which cover over 70% of the total treatment capacity of publicly owned treatment works (POTWs) in the Contiguous United States. These configurations were based on the Water Environment Research Foundation (WERF, now a part of the Water Research Foundation, WRF), report on net-zero energy solutions for WRRFs [3]_.

These simulation models have been implemented in the `werf EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/werf>`_. More details can be found in :ref:`the interactive page <wrrf_interactive>`.

.. figure:: ../images/wrrf_configs_light.png
   :class: only-light

.. figure:: ../images/wrrf_configs_dark.png
   :class: only-dark

   Distinguishing features of benchmark WRRF configurations. A WRRF configuration is referred to as a unique combination of a liquid code (a) and a solid code (b) as defined in [2]_. Configurations followed by a star (\*) are implemented in ``EXPOsan``.


Other WRRFs
***********

#. Validation of QSDsan implementation of the Anaerobic Digestion Model No. 1 (ADM1) [4]_

    * `adm EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm>`_


#. Validation of QSDsan implementation of the Activated Sludge Models (ASMs) [5]_

    * `asm EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/asm>`_


#. Conventional activated sludge process

    * Publication: Shoener et al., 2016 [6]_
    * `cas EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/cas>`_
    * `cas archived codes <https://github.com/QSD-Group/AnMBR>`_


#. Modular encapsulated two-stage anaerobic biological (METAB) system

    * Publication: Zhang et al., 2024 [7]_
    * `metab EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/metab>`_


#. EcoRecover system: microalgae-based tertiary P recovery process

    * Publication: Kim et al., 2025 [8]_
    * `pm2_batch EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_batch>`_
    * `pm2_ecorecover EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pm2_ecorecover>`_


----------

Non-sewered sanitation systems (NSSSs)
--------------------------------------

Biogenic Refinery
*****************
    - Publication: Rowles et al., 2022 [9]_
    - `biogenic_refinery EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biogenic_refinery>`_
    - `biogenic_refinery archived codes <https://github.com/QSD-Group/EXPOsan/releases/tag/archive%2FBR_OmniProcessor>`_

Bwaise
******
    - Publication: Trimmer et al., 2020 [10]_
    - `bwaise EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bwaise>`_
    - `bwaise Trimmer et al archived codes <https://github.com/QSD-Group/Bwaise-sanitation-alternatives>`_; `bwaise Li and Zhang et al archived codes <https://pypi.org/project/exposan/1.1.4>`_

Eco-San
*******
    - `eco_san EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/eco_san>`_

NEWgenerator
************
    - Publication: Watabe et al., 2023 [11]_
    - `new_generator EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/new_generator>`_
    - `new_generator archived codes <https://github.com/QSD-Group/EXPOsan-private/tree/newgen/exposan/newgen>`_

Reclaimer
*********
    - `reclaimer EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/reclaimer>`_

SCG Zyclonic
************

    * `scg_zyclonic EXPOsan module <https://github.com/QSD-Group/EXPOsan-private/tree/main/exposan/scg_zyclonic>`_


----------

Other Systems
-------------

#. Hydrothermal systems for biobinder and biofuels from food waste
    
    * Publication: Ahmand et al., 2026 [12]_
    * `biobinder EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/biobinder>`_


#. Hydroxyapatite (HAp) Synthesis from Urine
    
    * Publication:  Müller et al., 2025 [13]_
    * `hap EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/hap>`_


#. Hydrothermal systems for fuel and fertilizer production from wet organic wastes
    
    * Publication: Feng et al., 2024 [14]_
    * `htl EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/htl>`_


#. Point-of-use disinfection technologies
    
    * Publication: Elijah et al., 2024 [15]_
    * `pou_disinfection EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/pou_disinfection>`_


#. Hydrothermal systems for sustainable aviation fuel from food waste
    
    * Publication: Si et al., 2026 [16]_
    * `saf EXPOsan module <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/saf>`_


----------

Other Publications
------------------

The following publications also used QSDsan for process modeling and/or sustainability assessment:

#. `DMsan: A Multi-Criteria Decision Analysis Framework and Package to Characterize Contextualized Sustainability of Sanitation and Resource Recovery Technologies <https://doi.org/10.1021/acsenvironau.2c00067>`_

#. `Defining Economic and Environmental Typologies across 77 Countries to Prioritize Opportunities for Nonsewered Sanitation <https://doi.org/10.1021/acs.est.5c02064>`_

#. `Advancing the Economic and Environmental Sustainability of Rare Earth Element Recovery from Phosphogypsum <https://doi.org/10.1021/acs.est.5c04952>`_

#. `N2O as reactant rather than pollutant at wastewater treatment plants: Life Cycle Assessment and Techno-Economic Analysis of N2O-to-phenol <https://chemrxiv.org/doi/10.26434/chemrxiv-2026-m4vwz>`_


----------

.. References
.. [1] Li, Y.; Zhang, X.; Morgan, V. L.; Lohman, H. A. C.; Rowles, L. S.; Mittal, S.; Kogler, A.; Cusick, R. D.; Tarpeh, W. A.; Guest, J. S. QSDsan: An Integrated Platform for Quantitative Sustainable Design of Sanitation and Resource Recovery Systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289–2303. https://doi.org/10.1039/D2EW00455K.

.. [2] Zhang, X.; Rai, S.; Wang, Z.; Li, Y.; Guest, J. S. An Agile Benchmarking Framework for Wastewater Resource Recovery Technologies. npj Clean Water 2025, 9 (1), 4. https://doi.org/10.1038/s41545-025-00537-4.

.. [3] Tarallo, S., Shaw, A., Kohl, P. & Eschborn, R. A Guide to Net-Zero Energy Solutions for Water Resource Recovery Facilities. https://iwaponline.com/ebooks/book/293/ (2015).

.. [4] IWA Task Group for Mathematical Modelling of Anaerobic Digestion Processes. Anaerobic Digestion Model No.1 (ADM1); IWA Publishing, 2005. https://doi.org/10.2166/9781780403052

.. [5] Henze, M.; Gujer, W.; Mino, T.; van Loosedrecht, M. Activated Sludge Models ASM1, ASM2, ASM2d and ASM3. Water Intelligence Online 2006, 5 (0), 9781780402369–9781780402369. https://doi.org/10.2166/9781780402369.

.. [6] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S. Design of Anaerobic Membrane Bioreactors for the Valorization of Dilute Organic Carbon Waste Streams. Energy Environ. Sci. 2016, 9 (3), 1102–1112. https://doi.org/10.1039/C5EE03715H.

.. [7] Zhang, X.; Arnold, W. A.; Wright, N.; Novak, P. J.; Guest, J. S. Prioritization of Early-Stage Research and Development of a Hydrogel-Encapsulated Anaerobic Technology for Distributed Treatment of High Strength Organic Wastewater. Environ. Sci. Technol. 2024, 58 (44), 19651–19665. https://doi.org/10.1021/acs.est.4c05389.

.. [8] Kim, G.-Y.; Molitor, H. R.; Zhang, X.; Li, Y.; Shoener, B. D.; Schramm, S. M.; Morgenroth, E.; Snowling, S. D.; Hartnett, E.; Bradley, I. M.; Pinto, A. J.; Guest, J. S. Development of an Open-Source Process Simulator for Microalgae-Based Tertiary Phosphorus Recovery. npj Clean Water 2025, 9 (1), 13. https://doi.org/10.1038/s41545-025-00545-4.

.. [9] Rowles, L. S.; Morgan, V. L.; Li, Y.; Zhang, X.; Watabe, S.; Stephen, T.; Lohman, H. A. C.; DeSouza, D.; Hallowell, J.; Cusick, R. D.; Guest, J. S. Financial Viability and Environmental Sustainability of Fecal Sludge Treatment with Pyrolysis Omni Processors. ACS Environ. Au 2022, 2 (5), 455–466. https://doi.org/10.1021/acsenvironau.2c00022.

.. [10] Trimmer, J. T.; Lohman, H. A. C.; Byrne, D. M.; Houser, S. A.; Jjuuko, F.; Katende, D.; Banadda, N.; Zerai, A.; Miller, D. C.; Guest, J. S. Navigating Multidimensional Social–Ecological System Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement. Environ. Sci. Technol. 2020, 54 (19), 12641–12653. https://doi.org/10.1021/acs.est.0c03296.

.. [11] Watabe, S.; Lohman, H. A. C.; Li, Y.; Morgan, V. L.; Rowles, L. S.; Stephen, T.; Shyu, H.-Y.; Bair, R. A.; Castro, C. J.; Cusick, R. D.; Yeh, D. H.; Guest, J. S. Advancing the Economic and Environmental Sustainability of the NEWgenerator Nonsewered Sanitation System. ACS Environ. Au 2023, 3 (4), 209–222. https://doi.org/10.1021/acsenvironau.3c00001.

.. [12] Ahmad, A.; Kawale, H. D.; Summers, S.; Bogarin Cantero, B. C.; Allen, C. M.; Hajj, R. M.; Davidson, P. C.; Zhang, Y.; Li, Y. Financial Viability and Carbon Intensity of Hydrothermal Waste Valorization Systems for Bio-Based Asphalt Binder. Chemical Engineering Journal 2026, 528, 172283. https://doi.org/10.1016/j.cej.2025.172283.

.. [13] Müller, I. E.; Lin, A. Y. W.; Otani, Y.; Zhang, X.; Wu, Z.-Y.; Kisailus, D.; Mouncey, N. J.; Guest, J. S.; Rad, B.; Ercius, P.; Yoshikuni, Y. Cost-Effective Urine Recycling Enabled by a Synthetic Osteoyeast Platform for Production of Hydroxyapatite. Nat Commun 2025, 16 (1), 4216. https://doi.org/10.1038/s41467-025-59416-8.

.. [14] Feng, J.; Strathmann, T. J.; Guest, J. S. Hydrothermal-Based Wastewater Solids Management for Targeted Resource Recovery and Decarbonization in the Contiguous U.S. Environ. Sci. Technol. 2025. https://doi.org/10.1021/acs.est.5c07190.

.. [15] Elijah, B. C.; Ahmad, A.; Li, Y.; Plazas-Tuttle, J.; Rowles, L. S. Assessing the Relative Sustainability of Point-of-Use Water Disinfection Technologies for Off-Grid Communities. ACS Environ. Au 2024, 4 (5), 248–259. https://doi.org/10.1021/acsenvironau.4c00017.

.. [16] Si, B.; Wang, Z.; Watson, J.; Summers, S.; Li, Y.; Yu, S.; Yang, H.; Yang, Z.; Heyne, J. S.; Jiang, J.; Ren, Z. J.; Ma, H.; Wang, C.; Wang, P.; Zhang, Y. A Circular Hydrothermal Refinery for Sustainable Aviation Fuel from Food Waste. Nat Sustain 2026, 1–11. https://doi.org/10.1038/s41893-026-01848-1.

.. toctree::
   :hidden:

   wrrf_interactive
