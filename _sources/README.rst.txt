|logo| PAMtools   
===============

.. |logo| image:: _static/PAMtools_logo.png
   :width: 50

What is PAMtools?
-----------------

PAMtools is a set of helper functions written in Matlab and Python to export and plot single-molecule data that has been processed using PIE Analysis with Matlab (PAM) software for single-molecule fluorescence spectroscopy. The Matlab functions (**MatPAM**) are used to extract data from PAM's BurstBrowser and FCSfit, while the Python functionality (**PyPAM**) is dedicated to visualizing the extracted data for example in a jupyter notebook.


Download and install
--------------------

1. Make sure that you have Matlab (including the curvefit and Statistics and machine Learning toolbox installed).

2. Install [PAM](https://gitlab.com/PAM-PIE/PAM)

3. Clone or :download:`download <https://github.com/fdsteffen/PAMtools/archive/master.zip>` PAMtools into a directory of your choice. ::

    git clone https://github.com/fdsteffen/PAMtools.git

- Add the **MatPAM** folder to the Matlab search path
- Install the **PyPAM** module via pip ::

    pip install <path/to/PAMtools>

4. Copy the file in `matpam/profiles/profile PIE.mat`` into the profiles directory in PAM.

5. Copy the FCS models in `matpam/models` to the Models directory in PAM


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
