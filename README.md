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

1. Make sure that you have Matlab (including the curvefit and Statistics and machine Learning toolbox installed). **MatPAM** uses the function `jsonencode()` and therefore requires Matlab 2016b (9.1) or newer

2. Install [PAM](https://gitlab.com/PAM-PIE/PAM)

3. Clone PAMtools into a directory of your choice.
```
git clone https://github.com/fdsteffen/PAMtools.git
```

- Add the **MatPAM** folder to the Matlab search path
- Install the **PyPAM** module via pip
```
pip install .
```

4. Copy the file in `src/pamtools/matpam/profiles/profile PIE.mat` into the profiles directory in PAM.

5. Copy the FCS models in `src/pamtools/matpam/models` to the Models directory in PAM

---
<sup><a name="pymol">1</a></sup> W. Schrimpf, A. Barth, J. Hendrix, D.C. Lamb, *Biophys. J.* **2018**, *114*, 1518–1528. [![](https://img.shields.io/badge/DOI-10.1016/j.bpj.2018.02.035-blue.svg)](https://doi.org/10.1016/j.bpj.2018.02.035)
