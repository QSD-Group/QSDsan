Common Errors
=============

.. note::
   For modeling-time surprises (a number that looks wrong, a unit that sized weirdly, a sim that won't converge), see :doc:`Modeling Notes & Pitfalls <../tutorials/14_Modeling_Notes_and_Pitfalls>`.

.. _graphviz-installation:

``Graphviz`` Installation
-------------------------
Install `Graphviz <https://graphviz.org/download/>`_ if you want to use QSDsan's ``diagram`` methods.

Graphviz installation has two parts:

1. The Graphviz software, which provides the ``dot`` executable used to render diagrams.
2. The Python ``graphviz`` package, which lets Python call Graphviz.

For the Graphviz software, follow the `official Graphviz download instructions <https://graphviz.org/download/>`_ for your operating system. Then install the Python interface in your active Python environment:

.. code::

   pip install graphviz

To check that Graphviz is installed correctly, run:

.. code::

   dot -V

If ``dot`` is not found, restart your terminal/editor and try again. If it still fails, follow the Graphviz instructions for adding the Graphviz ``bin`` directory to your system ``PATH``.

When using :func:`diagram`, if you run into a ``graphviz`` error similar to:

   .. code:: bash

       FileNotFoundError: [Errno 2] No such file or directory: 'dot'


or

   .. code::

       ExecutableNotFound: failed to execute ['dot', '-Kdot', '-Tpng'], make sure the Graphviz executables are on your systems' PATH


or if you cannot get any diagram at all. It is likely that your ``graphviz`` is not configured correctly.

.. note::
    For the case where you don't have any diagram, the ``biosteam.preferences.raise_exception`` is set to ``False``, you can see the error by changing that to ``True``.

In most cases, this means either the Graphviz software is not installed or the ``dot`` executable is not available on your system ``PATH``.


``ModuleNotFoundError``
-----------------------
Sometimes (even though you have downloaded/cloned/installed ``qsdsan``), you still cannot see:

   .. code::

       ModuleNotFoundError: No module named 'qsdsan'


There are multiple possible reasons:

- If you use a virtual environment, make sure it is activated first (run the activation script for a ``.venv``, or ``conda activate <ENV NAME>`` for a ``conda`` env, replacing ``<ENV NAME>`` with the actual name of your environment).
- If you are using a cloned version of ``QSDsan`` for development, install it from the cloned repository with ``pip install -e ".[dev]"`` (or ``uv pip install -e ".[dev]"``).
- If you are using Jupyter Notebook

    - Make sure Jupyter is using the same environment where ``QSDsan`` is installed.
    - To add your active environment as a Jupyter kernel, run this in your command-line interface (CLI, e.g., a terminal or Anaconda Prompt; activate your environment first if you use one, and replace ``<KERNEL NAME>`` with the name you like):

        .. code::

            python -m ipykernel install --user --name <KERNEL_NAME>


        .. note::

            If you do not have ``ipykernel``, install it first with ``pip install ipykernel`` (or ``conda install ipykernel`` in a ``conda`` env).


        Then when you open the Jupyter Notebook, select the ``<KERNEL NAME>`` kernel when you create a new notebook you can find more details in this post about `enabling multiple kernels in Jupyter Notebook <https://medium.com/@ace139/enable-multiple-kernels-in-jupyter-notebooks-6098c738fe72>`_.


``underlying object has vanished``
----------------------------------
This error is related to ``numba`` caching, we haven't figured out the exact mechanism, but clearing cache will help resolve it. One/both of the following approaches should work:

1. Clear cache. Remove all ``.pyc``, ``.nbc``, and ``.nbi`` files from the relevant package or project directory:

   .. code::

       get-childitem . -recurse -include *.pyc, *.nbc, *.nbi | remove-item

2. Uninstalling and reinstalling a different version of ``numba``. Suppose you now have 0.58.1 and the newest version is 0.60.0, you can do:

   .. code::

       pip uninstall numba
       pip install --no-cache-dir numba==0.60.0

The ``--no-cache-dir`` option is to do a fresh installation rather than using previously downloaded packages. Note that you need to exit out your editor/any other programs that are currently using numba. Otherwise the uninstallation is incomplete, you might be prompted to do a manual removal, or this won't work.


``UnicodeDecodeError``
----------------------
When using non-English operating systems, you may run into errors similar to (cp949 is the case of Korean Windows):

   .. code::

       UnicodeDecodeError: 'cp949' codec can't decode byte Oxe2 in position 3426: multibyte sequence


To fix this, Windows users can look at this `thread <https://stackoverflow.com/questions/57131654/using-utf-8-encoding-chcp-65001-in-command-prompt-windows-powershell-window>`_ on updating the character encoding in the Windows console to UTF-8. We are not sure if this error will appear for Mac users, but let us know if you run into this and we will be happy to help with troubleshooting.


Version compatibility
---------------------
``QSDsan`` requires Python 3.12 or newer and pins specific minimum versions of ``biosteam``, ``thermosteam``, and other dependencies. The authoritative version table lives in `pyproject.toml <https://github.com/QSD-Group/QSDsan/blob/main/pyproject.toml>`_; the table here is a quick reference and may lag.

If you see import errors after upgrading one of those dependencies independently of ``QSDsan``, the most reliable fix is to reinstall ``QSDsan`` in a fresh environment so the resolver picks compatible versions of everything at once.
