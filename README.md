## Introduction
---

This repository is intended to provide a foundation for visualizing simulation data produced by the [BlackHoles@Home](https://blackholesathome.net/) project using python and VisIt.

## Setup
---

Getting Python to cooperate with VisIt can be a hastle depending on the operating system. Below are the basic steps for setting things up.

### Windows 11 installation
1. In the directory where you would like the repository, open terminal and run this codeblock. Alternatively run each line individually as needed.
```py
git clone https://github.com/tyndalestutz/bh_vis.git
cd bh_vis
Python -m venv .venv
.venv/scripts/Activate.ps1
pip install -U numpy "vtk<9.4" pillow pytest traitsui PyQt5 
pip install -U git+https://github.com/enthought/mayavi.git

```
 (note that you must have a version of python installed and added to PATH):