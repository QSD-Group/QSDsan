Contributing to QSDsan
======================

Below are some brief instructions on how to contribute to QSDsan. If you have any questions regarding the process, feel free to `submit an issue on GitHub <https://github.com/QSD-Group/QSDsan/issues>`_. Thank you in advance for your contribution!


Forking and Cloning
-------------------
#. Fork QSDsan by going to its `GitHub homepage <https://github.com/QSD-Group/QSDsan>`_ and click the "Fork" button at the top right corner.

#. Copy the link to **your** fork of QSDsan (you can find it in the green "Code" button on your forked GitHub QSDsan page), it should be something like:

	.. code:: bash

	    https://github.com/<YOUR_USERNAME>/QSDsan.git

	.. image:: ./images/code.png
		:height: 200
		:align: center


#. In your command prompt, navigate to your preferred location by using ``cd``, e.g.,

	.. code:: bash

	    cd research/coding


#. Clone QSDsan to your local by (use your own link copied from step 2):

	.. code:: bash

	    git clone https://github.com/<YOUR_USERNAME>/QSDsan.git --depth 1

	- If you don't have ``git``, follow the `instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ to install it.
	- The ``--depth-1`` flag is to tell ``git`` just clone the latest commit, you can change the depth number or just remove this flag completely, but then ``git`` will download more historical commits, which takes longer to clone and takes more space.

#. Add the root QSDsan as the upstream:

	.. code:: bash

	    git remote add upstream https://github.com/QSD-Group/QSDsan.git

#. Check your remote settings:

	.. code:: bash

	    git remote -v

	- This should show something like (origin is your fork and upstream is the root repository):

	.. code:: bash

		origin	https://github.com/<YOUR_USERNAME>/QSDsan.git (fetch)
		origin	https://github.com/<YOUR_USERNAME>/QSDsan.git (push)
		upstream	https://github.com/QSD-Group/QSDsan.git (fetch)
		upstream	https://github.com/QSD-Group/QSDsan.git (push)

Note
^^^^
GitHub has really detailed documentation on `forking <https://docs.github.com/en/github/getting-started-with-github/fork-a-repo>`_ (and almost everything else).


Developing Modules
------------------
#. Adding/modifying modules locally.

#. `Commit <https://git-scm.com/docs/git-commit>`_ your changes and concisely summarize your changes in the commit message.

	- You can have multiple `branches <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_ for different features.

#. Push your local changes to your remote fork:

	.. code:: bash

	    git push origin master

	- As your develop your contributions, the root repository may update, you should merge these changes and resolve any conflicts before your final push.

	.. code:: bash

	    git pull upstream master


Submitting Pull Request
-----------------------
#. Once you are satisfied with your changes and push all commits to your fork, go to you GitHub fork of ``QSDsan``, and submit a `pull request <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.

	- You can confirm that you have pulled all updates from the root repository if there's a message showing that your branch is X commits ahead of QSD-Group:master as shown (not X commits, Y commits behind).

	.. image:: ./images/commit.png
		:align: center

#. One of the Quantitative Sustainable Design Group members will review your changes and accept or discuss with you if edits are needed.


Documentation
-------------
Whenever new modules or functions are added, concise and thorough documents should be added with examples for `doctest <https://docs.python.org/3/library/doctest.html>`_. Please also include yourself (contact method is optional) to the list of contributors on the top of the module.

``QSDsan`` uses `numpydoc docstring style <https://numpydoc.readthedocs.io/en/latest/format.html>`_ with some modifications for better rendering. Some important notes:

- Both quotes ('') and double quotes ("") work.
- If you want some notes in your docstring, use `directives <https://docutils.sourceforge.io/docs/ref/rst/directives.html>`_ so that it can be rendered in sphinx.
	
	.. code::

		# This can be rendered by sphinx
		.. note::

			Something to notes

		# This won't be rendered by sphinx
		Notes
		-----

		# This can be rendered by sphinx but won't be recognized as numpydoc docstring
		Note
		----

- Use directives like ``:class:`package.class``` and ``:func:`class.function``` to indicate classes and functions, this will automatically add links to the corresponding documents.

	- Use single back ticks (``) in error messages and warnings since directives won't be rendered.

- If you want to refer to documents of other internal modules or external packages, please include it in the "See Also" section (refer to :class:`sanunits.AnaerobicDigestion` and :class:`Component` as examples).


Most of the documentations will be automatically generated through `sphinx's autodoc extension <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_. If your contribution involves new classes or modules, please add a new .rst file in docs/source/. and add it to the appropriate section in the ``index.rst`` file. You can refer to any of the existing files for examples.

Tutorials are prepared in `Jupyter Notebook <https://jupyter.org/>`_ and potential contributors are encouraged to use the `templates <https://github.com/QSD-Group/QSDsan/tree/master/docs/source/for_developers>`_.


Testing
-------
``QSDsan`` uses `Travis CI <https://travis-ci.com/>`_ for testing. A pull request will only be accepted when the branch has not conflicts with the root repository and all tests have been passed. More instructions on testing will be added.


