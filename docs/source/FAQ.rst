
Common Errors
=============

``graphviz``
------------
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


``UnicodeDecodeError``
----------------------
When using non-English operating systems, you may run into errors similar to (cp949 is the case of Korean Windows):

   .. code::

       UnicodeDecodeError: 'cp949' codec can't decode byte Oxe2 in position 3426: multibyte sequence


To fix this, Windows users can look at this `thread <https://stackoverflow.com/questions/57131654/using-utf-8-encoding-chcp-65001-in-command-prompt-windows-powershell-window>`_ on updating the character encoding in the Windows console to UTF-8. We are not sure if this error will appear for Mac users, but let us know if you run into this and we will be happy to help with troubleshooting.


Tips
====

Upgrade Python
--------------
``QSDsan`` is currently compatible with and tested for Python 3.7 and 3.8. However, with ``BioSTEAM`` moving to Python 3.8 (see this `issue <https://github.com/BioSTEAMDevelopmentGroup/biosteam/issues/56>`_), qsdsan may be only compatible with Python 3.8 and higher in the future. 

If you need to upgrade Python but having a lot of existing packages, creating a virtual environment may be the best way to avoid conflicts. If you are using ``conda``, its has related documentations on `Python upgrading <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html>`_.


Pickle Protocol
---------------
``QSDsan`` saves some of the default components and processes as `pickle <https://docs.python.org/3/library/pickle.html>`_ files to reduce the loading time, Python pickle has different protocols, and Protocol 5 is used in ``QSDsan``. The default ``pickle`` module in Python 3.5-3.7 uses Protocol 4 thus not compatible. For Python 3.5-3.7 users, ``QSDsan`` will prompt a warning to install the `package <https://pypi.org/project/pickle5/>`_ ``pickle5`` for compatibility. For Python 3.4 and below, longer loading time is expected as no pre-saved data files are used.