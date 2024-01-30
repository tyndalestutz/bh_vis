## Introduction
---

This repository is intended to provide a foundation for visualizing simulation data produced by the [BlackHoles@Home](https://blackholesathome.net/) project using python and VisIt.

## Setup
---

Getting Python to cooperate with VisIt can be a hastle depending on the operating system. Below are the basic steps for setting things up.

1. **Installing VisIt:** Installation files and instructions can be found [here](https://visit-dav.github.io/visit-website/releases-as-tables/), according to your operating system. 
1. **Appending Directories:**
    1. Locate site-packages:
        ```
        import sys
        sys.path.append("/path/to/visit/<version>/<architecture>/lib/site-packages")
        ```
    1. Import modules and launch:
        ```
        import visit
        visit.Launch()
        import visit
        ```