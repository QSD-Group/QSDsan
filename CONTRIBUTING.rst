Contributors and Guidelines
===========================

Contributors
------------
``QSDsan`` is developed and maintained by the Quantitative Sustainable Design Group and the broader community. Please refer to the `GitHub contributors <https://github.com/QSD-Group/QSDsan/graphs/contributors>`_ for the latest record of contributions.

Systems constructed using ``QSDsan`` are stored in the `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ repository; please refer to its `contributors <https://github.com/QSD-Group/EXPOsan/graphs/contributors>`_ and `commit history <https://github.com/QSD-Group/EXPOsan/commits/main/>`__ for contribution records related to system modules.

If you would like to join the effort, please review the guidelines below. If you have any questions about the process, feel free to `submit an issue on GitHub <https://github.com/QSD-Group/QSDsan/issues>`_. Thank you in advance for your contribution!


Authorship and Licensing
------------------------
The following guideline is adapted from `BioSTEAM <https://biosteam.readthedocs.io/en/latest/CONTRIBUTING.html#authorship>`_; we welcome inputs from the community for enhancement. If you feel that your contributions are not adequately acknowledged, please contact us.

#. Contributions should be acknowledged at the module level with a short description for:

   - Code development. The primary author is encouraged (but not required) to include contact info in the module.
   - Module development (i.e., math algorithms, code in other languages).
   - Instrumental comments and suggestions through discussion.

#. If any code or implementation was copied from a third party, note it in the module-level documentation.

#. Any third-party package copied into ``QSDsan`` must be strictly open-source (not copy-left nor open-access). If the license of the third-party package differs from ``QSDsan``, add the third-party license as an option (i.e., dual licensing).


Setting Up
----------
You can set up using the command line or `GitHub Desktop <https://desktop.github.com/>`_, a graphical alternative that is friendlier if you are new to the command line. ``QSDsan`` requires **Python 3.12 or newer**.

Via the command line
^^^^^^^^^^^^^^^^^^^^^

#. Fork ``QSDsan`` by going to its `GitHub homepage <https://github.com/QSD-Group/QSDsan>`_ and clicking the "Fork" button at the top right.

#. On your fork's page, click the green "Code" button and copy the HTTPS address; it looks like ``https://github.com/<YOUR_USERNAME>/QSDsan.git``.

#. In a terminal, navigate to where you keep your projects (for example, ``cd ~/projects`` on macOS/Linux or ``cd %USERPROFILE%\projects`` on Windows) and clone your fork. The ``git clone`` command creates a new ``QSDsan/`` folder inside the current directory, then ``cd QSDsan`` moves you into it:

   .. code:: bash

       git clone --depth=1 --no-single-branch https://github.com/<YOUR_USERNAME>/QSDsan.git
       cd QSDsan

   If you don't have ``git``, see the `installation instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_. The ``--depth=1`` flag makes a faster, smaller clone by fetching only recent history; ``--no-single-branch`` keeps all branches available. If you later need full history (for example, for ``git blame`` or ``git bisect``), run ``git fetch --unshallow``.

   The remaining setup steps all run from inside this ``QSDsan/`` folder. If you close your terminal and come back later, ``cd`` back into it before continuing.

#. Create and activate a virtual environment so ``QSDsan`` and its dependencies stay isolated from your other projects. Pick **one** of the two options below; both produce a ``.venv/`` folder inside ``QSDsan/``. ``QSDsan`` requires Python 3.12 or newer, so make sure the interpreter you use here meets that.

   **Option A: standard library (pip).**

   .. code:: bash

       # from inside the QSDsan/ folder
       python -m venv .venv

       # then activate it:
       #   Windows (PowerShell):  .venv\Scripts\Activate.ps1
       #   Windows (cmd.exe):     .venv\Scripts\activate.bat
       #   macOS / Linux:         source .venv/bin/activate

   **Option B:** `uv <https://docs.astral.sh/uv/>`_ **(faster, drop-in alternative).**

   ``uv`` manages Python interpreters itself, so you don't need a pre-installed Python 3.12. From inside the ``QSDsan/`` folder:

   .. code:: bash

       # 1. have uv download Python 3.12 (skip if you already have 3.12+ on PATH)
       uv python install 3.12

       # 2. create .venv/ pinned to that interpreter
       uv venv --python 3.12

       # 3. activate it (same paths as Option A):
       #   Windows (PowerShell):  .venv\Scripts\Activate.ps1
       #   Windows (cmd.exe):     .venv\Scripts\activate.bat
       #   macOS / Linux:         source .venv/bin/activate

   Once activated, your shell prompt typically shows ``(.venv)`` at the start. Every command in the rest of this guide assumes you are still inside ``QSDsan/`` with this environment activated.

#. Install ``QSDsan`` in editable mode with the development dependencies. This installs ``QSDsan`` from your local clone along with the packages needed for testing and building the documentation. Use the command matching the option you chose above:

   .. code:: bash

       # from inside the QSDsan/ folder, with .venv activated

       # Option A (pip):
       pip install -e ".[dev]"

       # Option B (uv):
       uv pip install -e ".[dev]"

   The ``.`` in ``-e "."`` refers to the current directory (``QSDsan/``), which is why this step must be run from inside the cloned folder.

#. Add the root ``QSDsan`` repository as the ``upstream`` remote, then check your remotes (still from inside ``QSDsan/``):

   .. code:: bash

       git remote add upstream https://github.com/QSD-Group/QSDsan.git
       git remote -v

   ``origin`` should point to your fork and ``upstream`` to ``QSD-Group/QSDsan``.

#. Pull the latest changes from upstream:

   .. code:: bash

       git pull upstream main

#. Create a branch for your work (branch names are case-sensitive):

   .. code:: bash

       git checkout -b <your-feature-name>

Via GitHub Desktop
^^^^^^^^^^^^^^^^^^

If you are new to the command line, `GitHub Desktop <https://desktop.github.com/>`_ offers a graphical interface. GitHub's guide on `cloning a repository to GitHub Desktop <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop>`_ includes screenshots.

#. Download and install GitHub Desktop.

#. Fork ``QSDsan`` from its `GitHub homepage <https://github.com/QSD-Group/QSDsan>`_.

#. On your fork's page, click the green "Code" button, choose "Open with GitHub Desktop", and pick where to clone it. When prompted whether you are contributing to the parent repository or working on your own, select "To contribute to the parent repository" (see GitHub's note on `fork behavior <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop#managing-fork-behavior>`_).

#. You can now `make commits <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/making-changes-in-a-branch/committing-and-reviewing-changes-to-your-project>`_ and `push <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/making-changes-in-a-branch/pushing-changes-to-github>`_ to your fork. To pull in later updates from ``QSDsan``, use the "Current Branch" menu to merge the corresponding ``upstream`` branch; you may occasionally need to `resolve conflicts <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/addressing-merge-conflicts/resolving-a-merge-conflict-on-github>`_.

#. Open a terminal in the cloned ``QSDsan/`` folder (in GitHub Desktop: **Repository → Open in Terminal**) and complete steps 4–5 of the command-line instructions above to create a virtual environment and install ``QSDsan`` in editable mode.

.. _also-contributing-to-exposan:

Also contributing to EXPOsan?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Systems built with ``QSDsan`` live in the separate `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ repository. If you only plan to change ``QSDsan`` itself, you can stop after the steps above. If you also plan to contribute to ``EXPOsan`` (for example, by adding a new system or modifying an existing one), use the side-by-side layout below instead so a change you make in ``QSDsan`` is immediately visible to ``EXPOsan``.

Put both clones inside a single workspace folder, with **one shared virtual environment at the workspace root**:

.. code-block:: text

   qsdsan-dev/          workspace folder (create this; open it in your editor)
   ├── QSDsan/          git clone of your QSDsan fork
   ├── EXPOsan/         git clone of your EXPOsan fork
   └── .venv/           shared virtual environment (one level above each repo)

The setup differs from the single-repo flow in two places:

#. Create the workspace folder first, then clone **both** forks into it (fork ``EXPOsan`` the same way you forked ``QSDsan``):

   .. code:: bash

       mkdir qsdsan-dev
       cd qsdsan-dev
       git clone --depth=1 --no-single-branch https://github.com/<YOUR_USERNAME>/QSDsan.git
       git clone --depth=1 --no-single-branch https://github.com/<YOUR_USERNAME>/EXPOsan.git

#. Create the virtual environment **at the workspace root** (``qsdsan-dev/``), not inside either repo, then install both packages in editable mode from that same folder:

   .. code:: bash

       # still in qsdsan-dev/, with .venv created and activated as in step 4 above

       # Option A (pip):
       pip install -e "./QSDsan[dev]"
       pip install -e "./EXPOsan[dev]"

       # Option B (uv):
       uv pip install -e "./QSDsan[dev]"
       uv pip install -e "./EXPOsan[dev]"

The remaining ``git remote add upstream`` / branching steps are unchanged, but run them **inside each repo separately** (``cd QSDsan`` to configure ``QSDsan``, ``cd ../EXPOsan`` to configure ``EXPOsan``); they have independent upstreams and independent feature branches.

.. note::

   Keep this workspace folder out of cloud-synced locations (OneDrive, Dropbox, iCloud). A virtual environment is tens of thousands of files; syncing them is slow and can corrupt the environment's tooling.

.. tip:: Keep the environment relocatable

   If you create the environment with ``uv``, add the ``--relocatable`` flag (``uv venv --relocatable``). A normal environment hard-codes the absolute path of its Python interpreter into every console-script launcher (``pytest``, ``sphinx-build``, and so on); renaming or moving ``qsdsan-dev/`` afterwards breaks those launchers with errors like *"failed to canonicalize script path"*. ``--relocatable`` makes them use relative paths, so the folder can be safely moved. If you hit this with an environment that already exists, recreate it.

Notes
^^^^^

#. Forking is the default workflow for all first-time contributors. Constant contributors who have made at least one successful and meaningful contribution through a fork may be given write access and invited to join the ``QSDsan`` team, after which branches can be used directly for easier syncing.
#. As ``QSDsan`` is public, all forks are public as well. If you would prefer a private fork, see the tip on creating a `private fork <FAQ.html#private-fork>`_.
#. GitHub has detailed documentation on `forking <https://docs.github.com/en/github/getting-started-with-github/fork-a-repo>`_ and related workflows.


Development Workflow
--------------------
Develop your contribution on the branch you created, documenting and testing as you go, then push to your fork and open a pull request.

Developing modules
^^^^^^^^^^^^^^^^^^^

Add or modify modules on your feature branch. Keep `commits <https://git-scm.com/docs/git-commit>`_ focused and write concise commit messages; you can use multiple `branches <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_ for different features. As you work, periodically merge in upstream changes and resolve any conflicts:

.. code:: bash

    git pull upstream main

Documentation
^^^^^^^^^^^^^

Whenever new modules or functions are added, include concise, thorough docstrings with examples for `doctest`_. Please also add yourself (contact method optional) to the list of contributors at the top of the module.

``QSDsan`` uses the `numpydoc docstring style <https://numpydoc.readthedocs.io/en/latest/format.html>`_ with some modifications for better rendering:

- Both single quotes ('') and double quotes ("") are fine.
- For notes in a docstring, use `directives <https://docutils.sourceforge.io/docs/ref/rst/directives.html>`_ so they render with `Sphinx <https://www.sphinx-doc.org/en/master/>`_:

  .. code::

      # Rendered by Sphinx and recognized as a docstring note
      .. note::

          Something to note.

      # NOT rendered by Sphinx
      Notes
      -----

      # Rendered by Sphinx but not recognized as a docstring note
      Note
      ----

- Use directives like ``:class:`package.class``` and ``:func:`class.function``` for classes and functions; these add links to the corresponding documents. Use single back ticks in error messages and warnings, since directives are not rendered there.
- To refer to other internal modules or external packages, include them in a "See Also" section (see :class:`qsdsan.unit_operations.AnaerobicDigestion` and :class:`qsdsan.Component` for examples).
- This `memo on reStructuredText and Sphinx <https://rest-sphinx-memo.readthedocs.io/en/latest/>`_ is a helpful reference.

Most documentation is generated automatically through `Sphinx's autodoc extension <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_. If your contribution adds new classes or modules, add a new ``.rst`` file under ``docs/source/`` and reference it from the appropriate section of ``index.rst`` (refer to existing files for examples). Tutorials are written as `Jupyter notebooks <https://jupyter.org/>`_; when adding or updating them, follow the structure of the existing notebooks in ``docs/source/tutorials``.

Build the documentation locally before submitting a pull request to confirm links and formatting are correct. This `YouTube walk-through <https://www.youtube.com/watch?v=oJsUvBQyHBs>`_ demonstrates the process.

Testing
^^^^^^^

``QSDsan`` uses `GitHub Actions <https://github.com/QSD-Group/QSDsan/actions>`_ to test every push and pull request. A pull request is accepted only when:

#. Meaningful contributions have been made.
#. The branch has no conflicts with the root repository.
#. All tests pass.

Run the tests locally before submitting. From your activated environment, in the cloned ``QSDsan`` directory:

.. code:: bash

    pytest    # if this doesn't work, try `python -m pytest`

This runs everything under ``tests/`` as well as the documentation examples via `doctest`_. A dot is a passing test and an ``F`` is a failure; tracebacks for any failures are printed to help you debug.

To confirm Python is using your local clone rather than another installed copy:

.. code:: bash

    python -c "import qsdsan; print(qsdsan.__path__)"

The printed path should point inside your cloned ``QSDsan`` folder.

.. figure:: ../../docs/source/images/pytest.png
   :width: 600
   :align: center

Committing and pushing
^^^^^^^^^^^^^^^^^^^^^^^

Push your branch to your fork. Pushing your work online backs up your history and makes it easier for maintainers to help you debug:

.. code:: bash

    git push origin <your-feature-name>


Submitting a Pull Request
-------------------------
#. Once you are satisfied with your changes and have pushed all commits to your fork, go to your GitHub fork of ``QSDsan`` and submit a `pull request <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_. Before submitting, pull the latest changes from the upstream repository (``git pull upstream main``) and resolve any conflicts, so your branch is ahead of the upstream ``main`` and not behind it.

#. A member of the Quantitative Sustainable Design Group will review your changes and accept them or discuss any needed edits with you.


For Maintainers
---------------
Maintainers reviewing pull requests that touch ``BioSTEAM``/``Thermosteam`` imports, stream/unit APIs, unit registries, or the ``Stream``/``SanStream``/``WasteStream`` taxonomy should follow the PR audit checklist at ``docs/maintainer/pr_audit_checklist.md``, which also explains how to install it as a Codex or Claude Code skill.


.. Links
.. _doctest: https://docs.python.org/3/library/doctest.html
