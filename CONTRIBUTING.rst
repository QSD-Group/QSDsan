Contributing to ``QSDsan``
==========================

Below are some brief instructions on how to contribute to ``QSDsan``. If you have any questions regarding the process, feel free to `submit an issue on GitHub <https://github.com/QSD-Group/QSDsan/issues>`_. Thank you in advance for your contribution!

Authorship
----------
The following guideline is adapted from `BioSTEAM <https://biosteam.readthedocs.io/en/latest/CONTRIBUTING.html#authorship>`_, we welcome inputs from the community for enhancement. If you feel that your contributions are not acknowledged or adequately acknowledged, please do contact us.

#. Contributions must be acknowledged at the module-level with a short description for:

	- Code development. The primary author is encouraged (but not required) to include contact info in the module.
	- Module development (i.e., math algorithms, codes in other languages).
	- Instrumental comments and suggestions through discussion.

#. All contributors will be added to the `author list <https://qsdsan.readthedocs.io/en/latest/AUTHORS.html>`_.

#. If any code or implementation was copied from a third party, it should be noted in the module-level documentation.

#. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). If license of the third-part package is different from ``QSDsan``, the module should add the third-party license as an option (i.e., dual licensing).


Forking and Cloning
-------------------

Via command-line interface
^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Fork ``QSDsan`` by going to its `GitHub homepage <https://github.com/QSD-Group/QSDsan>`_ and click the "Fork" button at the top right corner.

#. GitHub will open a new page showing your fork, click the green "Code" button on the top and copy the HTTPS address (there's a handy copy button next to the address), it should be something like:

	.. code:: bash

	    https://github.com/<YOUR_USERNAME>/QSDsan.git


#. In your command-line interface (e.g., Anaconda prompt, terminal), navigate to your preferred location by using ``cd``, e.g.,

	.. code:: bash

	    cd research/coding


#. Clone ``QSDsan`` to your local by (use your own link copied from step 2):

	.. code:: bash

	    git clone https://github.com/<YOUR_USERNAME>/QSDsan.git --depth=1

	- If you don't have ``git``, follow the `instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ to install it.
	- The ``--depth-1`` flag is to tell ``git`` just clone the latest commit, you can change the depth number or just remove this flag completely, but then ``git`` will download more historical commits, which takes longer time to clone and needs more space.

	.. note::
	 	
	 	This will only clone the main branch, if you want other branches, then use the ``no-single-branch`` flag, i.e.

		.. code:: bash

		    git clone https://github.com/<YOUR_USERNAME>/QSDsan.git --depth=1 --no-single-branch

		Without the ``no-single-branch`` flag, the reference of the remote branch (when you do ``git fetch``) is set to the main branch only (instead of all of the existing and future new branches), i.e., when you do

		.. code:: bash

		    git config --get remote.origin.fetch

		you will see

		.. code:: bash

		    +refs/heads/main:refs/remotes/origin/main

		Because it only tracks the main branch, so if didn't include the ``no-single-branch`` flag when cloning but later wanted to pull/push other branches, you will need to update the fetch reference to all branches using:


		.. code:: bash

		    git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"

		and you can double-check again to confirm the fetch reference has been updated.

#. Navigate into the cloned QSDsan:

	.. code:: bash

	    cd QSDsan

#. Install required packages:

	.. code:: bash

	    pip install –r requirements.txt


#. Add the root ``QSDsan`` as the upstream:

	.. code:: bash

	    git remote add upstream https://github.com/QSD-Group/QSDsan.git

#. Check your remote settings:

	.. code:: bash

	    git remote -v

	This should show something like (origin is your fork and upstream is the root repository):

	.. code:: bash

		origin	https://github.com/<YOUR_USERNAME>/QSDsan.git (fetch)
		origin	https://github.com/<YOUR_USERNAME>/QSDsan.git (push)
		upstream	https://github.com/QSD-Group/QSDsan.git (fetch)
		upstream	https://github.com/QSD-Group/QSDsan.git (push)

#. Pull in upstream changes:

	.. code:: bash

	    git pull upstream main

#. If you are working on a new feature (rather than some quick work like fixing a small bug), then it is recommended to checkout a new branch (note that branch names are case-sensitive):

	.. code:: bash

	    git checkout -b <REPLACE-ME-WITH-FEATURE-NAME>


Via GitHub Desktop
^^^^^^^^^^^^^^^^^^

If you are new to command-line interface, `GitHub Desktop <https://desktop.github.com/>`_ can be a good way to get started as it has a graphic interface, though less powerful.

To see screenshots of the different interface, visit GitHub's documentations on `Cloning a repository from GitHub to GitHub Desktop <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop>`_

#. Download and install GitHub Desktop.

#. Fork ``QSDsan`` by going to its `GitHub homepage <https://github.com/QSD-Group/QSDsan>`_ and click the "Fork" button at the top right corner.

#. GitHub will open a new page showing your fork, click the green "Code" button on the top and select "Open with GitHub Desktop".

#. GitHub Desktop will automatically open, and it will ask you where you want to clone it, select a place that you like.

#. Next, you will be prompted to select whether you want to contribute to the parent repository or for you own purpose, we would appreciate your contributing back to QSDsan, so please select "To contribute to the parent repository" :). You can read more about this, including how to change this setting, in this post about `fork behavior <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop#managing-fork-behavior>`_.

#. In the opened dialogue, click on the "Fetch origin" button on the top, then if you click the "Current Branch" button (next to the "Fetch origin" button), you should see a list of the branches on your fork (start with "origin", e.g., "origin/main") and those from the root repo managed by us (start with "upstream", e.g., "upstream/main"). All branches on your fork are copied from the corresponding branch from the root repo (i.e., "origin/main" copied from "upstream/main") at this moment. You can choose which one you would like to work on, if unsure, just select main (i.e., "origin/main").

#. You can work on your changes locally, `make commits <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/making-changes-in-a-branch/committing-and-reviewing-changes-to-your-project>`_, then `push <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/making-changes-in-a-branch/pushing-changes-to-github>`_ to your fork remote (i.e., on GitHub's website). Pushing them online would allow you to save/back up the history of your changes, and makes it super easy for us to help you debug.

#. In the future, whenever you want to merge changes from QSDsan (e.g., we just release a new feature), click on the "Current Branch" button, then click the "Choose a branch to merge into main" ("main" would be the name of the branch that you are working on) on the bottom of the drop-down, then select the branch from the root repo (starting with "upstream", e.g., "upstream/main") that you want to pull changes from, and click the "Create a merge commit" button on the bottom. Note that you can control whether Git does the pull ("merge", "rebase", etc.), check Git/GitHub's documentation if you want to know more. Also note that sometimes you need to `resolve conflicts <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/addressing-merge-conflicts/resolving-a-merge-conflict-on-github>`_ prior to merging.


Note
^^^^
#. We use fork as the default way for collaboration (i.e., for all first-time contributors). If you are a constant contributor and have independently made at least one successful and meaningful contribution through forking, you will be given the write access to ``QSDsan`` and you can use branch for easier code syncing. We will also jinvite you to join the ``QSDsan`` team.
#. GitHub has really detailed documentation on `forking <https://docs.github.com/en/github/getting-started-with-github/fork-a-repo>`_ (and almost everything else).
#. As QSDsan is public, all created forks would be public as well. We would appreciate if you make your work public and contribute back, but we understand it if you would like to create a private fork of QSDsan. To do so, please check our tip on creating the `private fork <https://qsdsan.readthedocs.io/en/latest/FAQ.html#private-fork>`_.


Developing Modules
------------------
#. Adding/modifying modules locally.

#. `Commit <https://git-scm.com/docs/git-commit>`_ your changes and concisely summarize your changes in the commit message.

	- You can have multiple `branches <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_ for different features.

#. Push your local changes to your remote fork:

	.. code:: bash

	    git push origin main # or the name of the new branch

	- As your develop your contributions, the root repository may update, you should merge these changes and resolve any conflicts before your final push.

	.. code:: bash

	    git pull upstream main


Submitting Pull Request
-----------------------
#. Once you are satisfied with your changes and push all commits to your fork, go to you GitHub fork of ``QSDsan``, and submit a `pull request <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.

	- You can confirm that you have pulled all updates from the root repository if there's a message showing that your branch is X commits ahead of QSD-Group:main (not X commits ahead, Y commits behind).

#. One of the Quantitative Sustainable Design Group members will review your changes and accept or discuss with you if edits are needed.


Documentation
-------------
Whenever new modules or functions are added, concise and thorough documents should be added with examples for `doctest`_. Please also include yourself (contact method is optional) to the list of contributors on the top of the module.

``QSDsan`` uses `numpydoc docstring style <https://numpydoc.readthedocs.io/en/latest/format.html>`_ with some modifications for better rendering. Some important notes:

- Both quotes ('') and double quotes ("") are good.
- If you want some notes in your docstring, use `directives <https://docutils.sourceforge.io/docs/ref/rst/directives.html>`_ so that it can be rendered by `Sphinx <https://www.sphinx-doc.org/en/master/>`_.
	
	.. code::

		# This can be rendered by Sphinx and as docstring
		.. note::

			Something to notes.

			[1] If you need to have a numbered list, be careful about line-wrapping and indentation.
			The start of the second line should align with the number, not the first character after the number. 

			[2] Second point.

		# This won't be rendered by Sphinx
		Notes
		-----

		# This can be rendered by Sphinx but won't be recognized as docstring
		Note
		----

- Use directives like ``:class:`package.class``` and ``:func:`class.function``` to indicate classes and functions, this will automatically add links to the corresponding documents.

	- Use single back ticks (``) in error messages and warnings since directives won't be rendered.

- If you want to refer to documents of other internal modules or external packages, please include it in the "See Also" section (refer to :class:`qsdsan.sanunits.AnaerobicDigestion` and :class:`qsdsan.Component` as examples).
- Here is a great `memo on reStructuredText and Sphinx <https://rest-sphinx-memo.readthedocs.io/en/latest/>`_.


Most of the documentations will be automatically generated through `Sphinx's autodoc extension <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_. If your contribution involves new classes or modules, please add a new .rst file in docs/source/. and add it to the appropriate section in the ``index.rst`` file. You can refer to any of the existing files for examples.


We recommend generating the documentation locally prior to push to GitHub/send in the pull request to make sure links, formatting, etc. are working properly. This `YouTube video <https://www.youtube.com/watch?v=oJsUvBQyHBs>`_ provides a good walk-through example/demonstration.


Tutorials are prepared in `Jupyter Notebook <https://jupyter.org/>`_ and potential contributors are encouraged to use the `templates <https://github.com/QSD-Group/QSDsan/tree/main/docs/source/for_developers>`_ which includes proper license and contribution information.


Testing
-------
``QSDsan`` uses `AppVeyor <https://www.appveyor.com/>`_ to test all pushes and pull requests. A pull request will only be accepted when:

#. Meaningful contributions have been made.
#. The branch has no conflicts with the root repository.
#. All tests have been passed.

You can run the test locally using `pytest <https://docs.pytest.org/en/6.2.x/>`_:

	.. code:: bash

	    python3 -m pytest

This runs all tests under the QSDsan/tests directory as well as all examples in the documentation through `doctest`_. Test results will be similar to the screenshot below, where a green dot indicates the test has been successfully passed and a red F indicates a failure. The number of dots and Fs indicate how many test functions or doctests are run for each moduel. Detailed error traceback on each failed test will be listed to help you fix the bug.

.. figure:: ../../docs/source/images/pytest.png
   :width: 600
   :align: center


.. Links
.. _doctest: https://docs.python.org/3/library/doctest.html