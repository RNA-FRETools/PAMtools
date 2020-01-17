
|logo| PAMtools   
===============

.. |logo| image:: _static/PAMtools_logo.png
   :width: 50


What is PAMtools
----------------

PAMtools is a set of helper functions written in Matlab and Python to export and plot single-molecule data that has been processed using PAM. The Matlab functions (**MatPAM**) are used to extract data from PAM's BurstBrowser and FCSfit, while the Python functionality (**PyPAM**) is dedicated to visualizing the extracted data for example in a jupyter notebook.


Download and install
--------------------

Clone or :download:`download <https://github.com/fdsteffen/PAMtools/archive/master.zip>` PAMtools into a directory of your choice. ::

    git clone https://github.com/fdsteffen/PAMtools.git

- Add the **MatPAM** folder to the Matlab search path
- Install the **PyPAM** module via pip ::

    pip install --user git+https://github.com/fdsteffen/PAMtools.git
    or
    pip install --user -e <path/to/PAMtools>


Dependencies
------------

- **MatPAM** uses the function `jsonencode()` and therefore requires Matlab 2016b (9.1) or newer

- **PyPAM** depends on the following Python packages:

    - numpy
    - scipy
    - pandas
    - matplotlib
    - seaborn
    - uncertainties
    - naturalcolors (get it from here_)

.. _here : https://github.com/fdsteffen/naturalcolors.git


Getting started
---------------

For an introduction into the functionality of PyPAM visit the :doc:`tutorial <tutorial/pamtools_tutorial>`. The Jupyter Notebook can be downloaded :download:`here <tutorial/pamtools_tutorial.ipynb>`.