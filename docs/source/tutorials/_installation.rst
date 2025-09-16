Extended Installation Instructions
==================================

If you are new to Python and not even sure how to set up your Python environment, please also refer to our `beginner tutorials <https://uofi.box.com/s/49wf5usk5hz3gdjmcswo3voiokmbzekb>`_ (work in progress).

Note that this is just ONE way to set up your environment and cloning the latest ``QSDsan``-related packages (not the pip-installed version which are only released periodically), not THE way.

#. Download and install `Anaconda <https://www.anaconda.com/products/distribution>`_

   - You will see the Python version that comes with the Anaconda when you download it.``QSDsan`` is mainly developed on Python 3.12 at the time that this tutorial was written.

#. Make a new conda environment with the same version of Python as the base environment in the Anaconda you just installed, e.g., for Python 3.12, in your command line interface (CLI, Anaconda Prompt for Windows and terminal for Mac):

   .. code::

      conda create --name qsdsan python=3.12
      
      conda activate qsdsan

   - The qsdsan is simply the name for the environment, you can change it to other ones that make more sense (e.g., dev), all of the following CLI codes should be executed on this environment

#. Install `Spyder <https://www.spyder-ide.org/>`_, in your CLI

   .. code::
      
      pip install spyder

#. Download GitHub Desktop if you are not comfortable with git in CLI, but it does not have all the capacities of CLI
#. Clone repositories, you'll want to clone the qsdsan branches of ``BioSTEAM``/``Thermosteam`` (dependencies of QSDsan that are being actively developed), and the main branches of ``QSDsan``/``EXPOsan``. In your CLI (``QSDsan`` as the example):

   .. code::
      
      git clone https://github.com/QSD-Group/QSDsan.git --depth=1 --no-single-branch

   - The ``depth`` flag allows you to clone only the most recent commits so the repo size will be smaller; the ``no-single-branch`` flag allows you to switch between the different branches
   - If you need to clone a branch that is not the default branch (e.g., the qsdsan branch of ``BioSTEAM``, add ``@<BRANCH_NAME>`` after ``.git``, i.e., ``git clone https://github.com/BioSTEAMDevelopmentGroup/biosteam.git@qsdsan --depth=1 --no-single-branch``)
   - `GitHub Desktop <https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-a-repository-from-github-to-github-desktop>`_ instructions

#. Install graphviz (you may run into all sorts of issues, search engines are very helpful…)

   - `Windows <https://forum.graphviz.org/t/new-simplified-installation-procedure-on-windows/224#format-svg-not-recognized-use-one-of>`_ instructions 
      
      - Read the instructions… or it might not work
      - You might need to add Graphviz/bin to Spyder’s PYTHONPATH manager

   - Mac instructions

      .. code::
         
         pip install graphviz

         conda install graphviz

   - If this still won’t work (i.e., you cannot get diagrams using ``QSDsan``), check out this `FAQ <https://qsdsan.readthedocs.io/en/latest/FAQ.html#graphviz>`_

#. Install dependency packages

   .. code::

      pip install exposan
      pip uninstall biosteam thermosteam qsdsan exposan

#. Add the path to ALL of the repositories you cloned in Spyder
#. Open Spyder

   - In your CLI, just do ``spyder``
   - Read all of the instructions upon opening Spyder
   - Add the path in Spyder’s PYTHONPATH manager