Styling and Citation
====================

``QSDsan`` vs. ``qsdsan``
-------------------------
``QSDsan`` vs. ``qsdsan``? We prefer to use the capitalized version when not in coding settings (e.g., ``QSDsan`` instead of ``qsdsan``) because:

- It refers to the platform, not just the core package (i.e., it includes the entire ecosystem that supports the core package).
- We style the name to convey the name's meaning (e.g., the "QSD" part stands for "quantitative sustainable design").

But names of the actual packages are all in lower cases per `PEP-8 <https://www.python.org/dev/peps/pep-0008/#package-and-module-names>`_:

   *Modules should have short, all-lowercase names. Underscores can be used in the module name if it improves readability. Python packages should also have short, all-lowercase names, although the use of underscores is discouraged.*


Citing QSDsan
-------------
If you use ``QSDsan`` in published work, please cite (Li and Zhang et al., 2022):

  Li, Y.; Zhang, X.; Morgan, V. L.; Lohman, H. A. C.; Rowles, L. S.; Mittal, S.; Kogler, A.; Cusick, R. D.; Tarpeh, W. A.; Guest, J. S. QSDsan: An Integrated Platform for Quantitative Sustainable Design of Sanitation and Resource Recovery Systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289–2303. https://doi.org/10.1039/D2EW00455K.


BibTeX:

.. code:: bibtex

    @article{li2022qsdsan,
      author  = {Li, Yalin and Zhang, Xinyi and Morgan, Victoria L. and Lohman, Hannah A. C. and Rowles, Lewis S. and Mittal, Smiti and Kogler, Anna and Cusick, Roland D. and Tarpeh, William A. and Guest, Jeremy S.},
      title   = {{QSDsan}: An Integrated Platform for Quantitative Sustainable Design of Sanitation and Resource Recovery Systems},
      journal = {Environmental Science: Water Research \& Technology},
      year    = {2022},
      volume  = {8},
      number  = {10},
      pages   = {2289--2303},
      doi     = {10.1039/D2EW00455K}
    }

That is the journal publication. If you would rather cite the software package itself, cite the archived release on Zenodo:

- ``QSDsan``: https://doi.org/10.5281/zenodo.20256569
- ``EXPOsan``: https://doi.org/10.5281/zenodo.20256578

For specific subpackages or unit models, also cite the relevant primary references listed in their docstrings.
