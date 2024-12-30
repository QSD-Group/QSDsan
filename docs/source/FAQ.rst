.. _faq:

FAQ
===

Common Errors
-------------

``graphviz``
************
When using :func:`diagram`, if you run into a ``graphviz`` error similar to:

   .. code:: bash

       FileNotFoundError: [Errno 2] No such file or directory: 'dot'


or

   .. code::

       ExecutableNotFound: failed to execute ['dot', '-Kdot', '-Tpng'], make sure the Graphviz executables are on your systems' PATH


or if you cannot get any diagram at all. It is likely that your ``graphviz`` is not configured correctly.

.. note::
    For the case where you don't have any diagram, the ``biosteam.RAISE_GRAPHVIZ_EXCEPTION`` is set to ``False``, you can see the error by changing that to ``True``.

This `post <https://stackoverflow.com/questions/35064304/runtimeerror-make-sure-the-graphviz-executables-are-on-your-systems-path-aft>`_ provides a lot of useful information, and this normally can be solved by:

    .. code::

       conda install graphviz # if you are using conda


or

    .. code::

       brew install graphviz # if you are using brew


.. note::

    If you have already installed graphviz (both the actual software and the Python interface) but still getting the same error, your probably need to add the path of the graphviz software to your system path. To do that, you need to firstly locate where the graphviz software is, add the graphviz path to your system path (for Windows, the post above has instruction on how to add to your path; for macOS, you add ``export PATH="<REPLACE_WITH_GRAPHVIZ_PATH>:$PATH"`` to your shell profile).


``ModuleNotFoundError``
***********************
Sometimes (even though you have downloaded/cloned/installed ``qsdsan``), you still cannot see:

   .. code::

       ModuleNotFoundError: No module named 'qsdsan'


There are multiple possible reasons:

- If you have multiple conda environments, make sure you firstly do ``conda activate <ENV NAME>`` (replace ``<ENV NAME>`` with the actual name of your environment) to activate the environment.
- If you are using the downloaded/cloned version of qsdsan, make sure you have added the path to the cloned version to your system path (more details on the tutorial `Helpful Basics <https://qsdsan.readthedocs.io/en/latest/tutorials/1_Helpful_Basics.html>`_).
- If you are using Jupyter Notebook
    
    - If you are using the downloaded/cloned version of ``qsdsan``, note that Jupyter Notebook does not know about the path you configured in other editors (e.g., Spyder), so you may need to change directory (with ``os.chdir``) or set ``sys.path``.
    - If you are using pip-installed ``qsdsan``, try to do this in your command-line interface (CLI, e.g., Anaconda prompt, terminal; firstly do ``conda activate <ENV NAME>`` if you are using an virtual environment, and replace ``<KERNEL NAME>`` with the name you like):

        .. code::

            python -m ipykernel install --user --name <KERNEL_NAME>


        .. note::

            If you do not have ``ipykernel``, firstly do ``conda install ipykernel``


        Then when you open the Jupyter Notebook, select the ``<KERNEL NAME>`` kernel when you create a new notebook you can find more details in this post about `enabling multiple kernels in Jupyter Notebook <https://medium.com/@ace139/enable-multiple-kernels-in-jupyter-notebooks-6098c738fe72>`_.


``underlying object has vanished``
**********************************
This error is related to ``numba`` caching, we haven't figured out the exact mechanism, but clearing cache will help resolve it. One/both of the following approaches should work:

1. Clear cache. Remove all ``.pyc``, ``.nbc``, and ``.nbi`` files, you can do this in your CLI using (replace <DIR> with the directory to your ``thermosteam``, ``biosteam``, ``qsdsan``, and ``exposan`` directory):

   .. code::

       get-childitem . -recurse -include *.pyc, *.nbc, *.nbi | remove-item

2. Uninstalling and reinstalling a different version of ``numba``. Suppose you now have 0.58.1 and the newest version is 0.60.0, you can do:

   .. code::

       pip uninstall numba
       pip install --no-cache-dir numba==0.60.0

The ``--no-cache-dir`` option is to do a fresh installation rather than using previously downloaded packages. Note that you need to exit out your editor/any other programs that are currently using numba. Otherwise the uninstallation is incomplete, you might be prompted to do a manual removal, or this won't work.


``UnicodeDecodeError``
**********************
When using non-English operating systems, you may run into errors similar to (cp949 is the case of Korean Windows):

   .. code::

       UnicodeDecodeError: 'cp949' codec can't decode byte Oxe2 in position 3426: multibyte sequence


To fix this, Windows users can look at this `thread <https://stackoverflow.com/questions/57131654/using-utf-8-encoding-chcp-65001-in-command-prompt-windows-powershell-window>`_ on updating the character encoding in the Windows console to UTF-8. We are not sure if this error will appear for Mac users, but let us know if you run into this and we will be happy to help with troubleshooting.


Tips
----

Archive Branch
**************
If you want to archive a branch but don't want to let it clutter your branch list, you can `archive it <https://stackoverflow.com/questions/1307114/how-can-i-archive-git-branches>`_. Essentially, you would need to

.. code::

    git checkout <BRANCH_TO_BE_ARCHIVED>
    git tag archive/<BRANCH_TO_BE_ARCHIVED> # "archive/<BRANCH_TO_BE_ARCHIVED>" will be the tag name, you can change it however you like
    git push origin archive/<BRANCH_TO_BE_ARCHIVED>

Then you'll see the new tag appears on GitHub and you can safely remove the archived branch from local and remote.


Pickle Protocol
***************
``QSDsan`` saves some of the default components and processes as `pickle <https://docs.python.org/3/library/pickle.html>`_ files to reduce the loading time, Python pickle has different protocols, and Protocol 5 is used in ``QSDsan``. The default ``pickle`` module in Python 3.5-3.7 uses Protocol 4 thus not compatible. For Python 3.5-3.7 users, ``QSDsan`` will prompt a warning to install the `package <https://pypi.org/project/pickle5/>`_ ``pickle5`` for compatibility. For Python 3.4 and below, longer loading time is expected as no pre-saved data files are used.


Private Fork
************
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

        Don't forget to firstly navigate to the ``QSDsan`` folder by ``cd QSDsan``

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
#. After you cloned ``QSDsan``, you'll need to configure your system path to make sure that you are importing the cloned ``QSDsan``, which means you might need to uninstalled any ``pip``-installed version and add the cloned path to your IDE (e.g., Spyder).


Upgrade Python
**************
``QSDsan`` is currently compatible with and tested for Python 3.7 and 3.8. However, with ``BioSTEAM`` moving to Python 3.8 (see this `issue <https://github.com/BioSTEAMDevelopmentGroup/biosteam/issues/56>`_), qsdsan may be only compatible with Python 3.8 and higher in the future. 

If you need to upgrade Python but having a lot of existing packages, creating a virtual environment may be the best way to avoid conflicts. If you are using ``conda``, its has related documentations on `Python upgrading <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html>`_.


Styling
-------
``QSDsan`` vs. ``qsdsan``? We prefer to use the capitalized version when not in coding settings (e.g., ``QSDsan`` instead of ``qsdsan``) because:

- It refers to the platform, not just the core package (i.e., it includes the entire ecosystem that supports the core package).
- We style the name to convey the name's meaning (e.g., the "QSD" part stands for "quantitative sustainable design").

But names of the actual packages are all in lower cases per `PEP-8 <https://www.python.org/dev/peps/pep-0008/#package-and-module-names>`_:

   *Modules should have short, all-lowercase names. Underscores can be used in the module name if it improves readability. Python packages should also have short, all-lowercase names, although the use of underscores is discouraged.*