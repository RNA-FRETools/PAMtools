# <img src="docs/source/_static/PAMtools_logo.png" width="200">   

[![Docs Status](https://github.com/fdsteffen/PAMtools/workflows/PAMtools%20docs/badge.svg)](https://github.com/fdsteffen/PAMtools/actions)

What is PAMtools?
-----------------

PAMtools is a set of helper functions written in Matlab and Python to export and plot single-molecule data that has been processed using the single-molecule software package "PIE Analysis with MATLAB" (PAM) (https://gitlab.com/PAM-PIE).<sup>[1](#pymol)</sup> The Matlab functions (**MatPAM**) are used to extract data from PAM's BurstBrowser and FCSfit, while the Python API (**PyPAM**) is used to visualize and plot the extracted data for example in a jupyter notebook.

Documentation
-------------

PAMtools is documented [here](https://fdsteffen.github.io/PAMtools/).


Download and install
--------------------

Clone or [download](https://github.com/fdsteffen/PAMtools/archive/master.zip) PAMtools into a directory of your choice.
```
git clone https://github.com/fdsteffen/PAMtools.git
```

- Add the **MatPAM** folder to the Matlab search path
- Install the **PyPAM** module via pip
```
pip install --user git+https://github.com/fdsteffen/PAMtools.git
or
pip install --user -e <path/to/PAMtools>
```

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
    - naturalcolors (get it from [here](https://github.com/fdsteffen/naturalcolors.git))

---
<sup><a name="pymol">1</a></sup> W. Schrimpf, A. Barth, J. Hendrix, D.C. Lamb, *Biophys. J.* **2018**, *114*, 1518–1528. [![](https://img.shields.io/badge/DOI-10.1016/j.bpj.2018.02.035-blue.svg)](https://doi.org/10.1016/j.bpj.2018.02.035)
