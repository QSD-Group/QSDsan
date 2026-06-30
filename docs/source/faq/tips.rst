Tips
====

Setting up a development environment
------------------------------------
If you cloned ``QSDsan`` (and/or ``EXPOsan``) and want to work on the code, install in editable mode so Python imports your local copy:

.. code:: bash

    pip install -e ".[dev]"

Use a dedicated virtual environment (``.venv`` or a ``conda`` env) to keep the editable install isolated. To use that environment in Jupyter:

.. code:: bash

    python -m ipykernel install --user --name qsdsan-dev

Then pick ``qsdsan-dev`` as the kernel when you open a tutorial notebook. If you also work on ``EXPOsan``, repeat the editable install in its clone too; sibling editable installs let cross-repo changes show up immediately.

On Windows, paths with spaces (e.g., ``C:\Users\Your Name\Documents``) occasionally break tools that don't quote them. If a build or test fails with a "file not found" error involving such a path, try a path without spaces.


Archive Branch
--------------
If you want to archive a branch but don't want to let it clutter your branch list, you can `archive it <https://stackoverflow.com/questions/1307114/how-can-i-archive-git-branches>`_. Essentially, you would need to

.. code::

    git checkout <BRANCH_TO_BE_ARCHIVED>
    git tag archive/<BRANCH_TO_BE_ARCHIVED> # "archive/<BRANCH_TO_BE_ARCHIVED>" will be the tag name, you can change it however you like
    git push origin archive/<BRANCH_TO_BE_ARCHIVED>

Then you'll see the new tag appears on GitHub and you can safely remove the archived branch from local and remote.


Pickle Protocol
---------------
``QSDsan`` saves some of the default components and processes as `pickle <https://docs.python.org/3/library/pickle.html>`_ files to reduce loading time. If loading cached data fails after upgrading ``QSDsan`` or Python, reinstall ``QSDsan`` in a fresh environment or clear the relevant cached files and try again.


Private Fork
------------
While ``QSDsan`` (and other supporting packages such as ``EXPOsan``) will stay open-source, it is totally understandable that you may want to create private forks of these packages (e.g., because of non-disclosure agreement).

However, GitHub does not allow you to directly create a private fork (or more accurately, this is a separate repo mirror the public repo ``QSDsan``). You can follow these steps for a work-around (modified from an original post `here <https://gist.github.com/0xjac/85097472043b697ab57ba1b1c7530274>`_, you need to do all following in your command-line interface):

#. Create a bare clone of the repository (this is temporary and will be removed):

    .. code::

        git clone --bare https://github.com/QSD-Group/QSDsan.git

    .. note::

        You should firstly navigate (i.e., ``cd``) to wherever you want the repository to be saved.

#. `Create a new private repository on Github <https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository>`_ and name it ``QSDsan`` (this name actually doesn't matter too much and you can use alternatives that you like, but you'll need to update the clone address below).
#. Mirror-push your bare clone to your new ``QSDsan`` repository (replace ``<YOUR_USERNAME>`` with your actual Github username in the url below, without the ``<>``):

    .. code::

        cd QSDsan.git
        git push --mirror https://github.com/<YOUR_USERNAME>/QSDsan.git

#. Remove the temporary local repository you created in step 1 (since we already pushed it to remote).

    .. code::

        cd ..
        rm -rf QSDsan.git

#. You can now clone your ``QSDsan`` repository to your local.

    .. code::

        git clone https://github.com/<YOUR_USERNAME>/QSDsan.git

#. It's also recommend to add the root ``QSDsan`` repo as remote to fetch future changes. Make sure you also disable push on the remote:

    .. code::

        git remote add upstream https://github.com/QSD-Group/QSDsan.git
        git remote set-url --push upstream DISABLED

    .. note::

        Don't forget to firstly navigate to the ``QSDsan`` directory by ``cd QSDsan``

#. To double-check things have been set up correctly, you can check the remote url using ``git remove -v``, and you should see something like:

    .. code::

        origin  https://github.com/<YOUR_USERNAME>/QSDsan.git (fetch)
        origin  https://github.com/<YOUR_USERNAME>/QSDsan.git (push)
        upstream    https://github.com/QSD-Group/QSDsan.git (fetch)
        upstream    DISABLE (push)

#. In the future, you'll want to push to ``origin`` to update your remote fork. To pull updates from the root ``QSDsan`` (i.e., ``upstream``):

    .. code::

        git fetch upstream
        git rebase upstream/main

**Other notes**

#. If you have never used ``git`` in your CLI, GitHub would ask for authentication and requires you create to a personal access token (instead of using your username and password), follow the instructions from `GitHub <https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token>`_ to create the token.
#. For Mac users, you'll probably run into an error related to ``/Library/Developer/CommandLineTools`` if you don't have Xcode Command Line (i.e., ``xcode-select``), follow these `instructions <https://www.freecodecamp.org/news/install-xcode-command-line-tools/>`_ to install it. Note that as you can see in the linked post, even the ``xcode-select``, which is much smaller than the full Xcode app, requires 1GB+ space.
#. After you clone ``QSDsan``, install it from the cloned repository with ``pip install -e ".[dev]"`` so Python imports your local copy.


Upgrade Python
--------------
``QSDsan`` currently requires Python 3.12 or newer.

If you need to upgrade Python but have a lot of existing packages, creating a virtual environment may be the best way to avoid conflicts. If you are using ``conda``, it has related documentation on `Python upgrading <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html>`_.


EXPOsan orientation
-------------------
``EXPOsan`` (the EXPOsition of sanitation and resource recovery systems) is a companion package that catalogs published sanitation and resource-recovery systems built on ``QSDsan``. The packages are released together but live in `separate <https://github.com/QSD-Group/QSDsan>`_ `repos <https://github.com/QSD-Group/EXPOsan>`_.

To run an EXPOsan system:

.. code:: python

    from exposan import bsm1
    bsm1.load()
    bsm1.sys.simulate()
    bsm1.sys.diagram()

See the `EXPOsan documentation <https://exposan.readthedocs.io/>`_ for the full system catalog.
