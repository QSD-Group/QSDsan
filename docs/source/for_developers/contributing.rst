Contributing to QSDsan
======================

Forking and cloning
-------------------
#. Fork ``QSDsan`` by going to its `GitHub homepage <https://github.com/QSD-Group/QSDsan>` and click the "Fork" button at the top right corner.

#. Copy the link to **your** fork of ``QSDsan`` (you can find it in the green "Code" button on your forked GitHub ``QSDsan`` page), it should be something like:

	.. code:: bash

	    https://github.com/<your_user_name>/QSDsan.git

	.. image:: ./images/code.png
		:height: 200
		:align: center


#. In your command prompt, navigate to your preferred location by using ``cd``, e.g.,

	.. code:: bash

	    cd research/coding


#. Clone ``QSDsan`` to your local by (use your own link copied from the previous step):

	.. code:: bash

	    git clone https://github.com/<your_user_name>/QSDsan.git

	- If you don't have ``git``, follow the `instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ to install it.


Documentation and Testing
-------------------------
Whenever new modules or functions are added, concise and thorough documents should be added with examples for `doctest <https://docs.python.org/3/library/doctest.html>`_. A pull request will only be accepted when all tests have been passed. More instructions on documentation and testing will be added.



Templates
---------
Templates for code and tutorials are available `here <https://github.com/QSD-Group/QSDsan/tree/master/docs/source/for_developers>`_.

