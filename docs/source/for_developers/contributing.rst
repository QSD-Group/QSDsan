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
#. Once you are satisfied with your changes and push all commits to your fork, go to you GitHub fork of QSDsan, and submit a `pull request <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.

	- You can confirm that you have pulled all updates from the root repository if there's a message showing that your branch is X commits ahead of QSD-Group:master as shown (not X commits, Y commits behind).

	.. image:: ./images/commit.png
		:align: center

#. One of the Quantitative Sustainable Design Group members will review your changes and accept or discuss with you if edits are needed.


Documentation and Testing
-------------------------
Whenever new modules or functions are added, concise and thorough documents should be added with examples for `doctest <https://docs.python.org/3/library/doctest.html>`_. A pull request will only be accepted when the branch has not conflicts with the root repository and all tests have been passed. More instructions on documentation and testing will be added.


Templates
---------
Templates for code and tutorials are available `here <https://github.com/QSD-Group/QSDsan/tree/master/docs/source/for_developers>`_.

