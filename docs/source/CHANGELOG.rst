==========
Change Log
==========

This document records notable changes to `QSDsan <https://github.com/QSD-Group/QSDsan>`_. We aim to follow `Semantic Versioning <https://semver.org/>`_.


Ongoing
-------
- ``LCA.get_normalized_impacts`` was replaced by ``LCA.get_allocated_impacts`` for flexible allocation options.
- Reformatted all documents, added instructions on documentation.
- Added brief instructions on contributing and code of conduct.
- Updated UML diagram.


`0.0.3`_ (2021-01-10)
---------------------
- More flexible setting of ``ImpactItem`` for ``WasteStream``.
- Add status badge to README.rst
- Add CHANGELOG.rst
- Tutorial updates:

	- New:
		- ``TEA`` and ``LCA``
	- Updated:
		-  ``Component`` and ``WasteStream``
		-  ``SanUnit`` and ``System``


`0.0.2`_ (2021-01-07)
---------------------
- Added the all three sanitation scenarios as described in `Trimmer et al.`_, including uncertainty/sensitivity analyses with tutorial.
- Inclusion of GPX models for estimation of ``WasteStream`` properties.
- Live documentation for the `stable package`_ and `beta version`_.
- New classes:

    - All units in `Trimmer et al.`_
    - ``Descriptor`` and ``Checker`` decorators to check user-input values.
    - ``AttrSetter``, ``DictAttrSetter``, and ``FuncGetter`` for batch-setting of uncertainty analysis parameters.

- Added ``save_report`` function to ``LCA`` class for report exporting.


`0.0.1`_ (2020-12-23)
---------------------
- First public release.


.. Other links
.. _stable package: https://qsdsan.readthedocs.io/en/latest/
.. _beta version: https://qsdsan-beta.readthedocs.io/en/latest/
.. _Trimmer et al.: https://doi.org/10.1021/acs.est.0c03296

.. Commit links
.. _0.0.3: https://github.com/QSD-Group/QSDsan/commit/e20222caccc58d9ee414ca08d8ec55f3a44ffca7
.. _0.0.2: https://github.com/QSD-Group/QSDsan/commit/84653f5979fbcd76a80ffb6b22ffec1c5ca2a084
.. _0.0.1: https://github.com/yalinli2/QSDsan/commit/f95e6172780cfe24ab68cd27ba19837e010b3d99

