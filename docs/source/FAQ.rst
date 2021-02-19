Common Errors
=============

``graphviz``
------------
When using :func:`System.diagram`, If you run into a ``graphviz`` error similar to:

   .. code:: bash

       FileNotFoundError: [Errno 2] No such file or directory: 'dot'


or

   .. code:: bash

       ExecutableNotFound: failed to execute ['dot', '-Kdot', '-Tpng'], make sure the Graphviz executables are on your systems' PATH


It is likely that your ``graphviz`` is not configured correctly. This `post <https://stackoverflow.com/questions/35064304/runtimeerror-make-sure-the-graphviz-executables-are-on-your-systems-path-aft>`_ provides a lot of useful information, and this normally can be solved by:

    .. code:: bash

       conda install graphviz # if you are using conda


or

    .. code:: bash

       brew install graphviz # if you are using brew

