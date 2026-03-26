# UFS_domain_select

A Cartopy script to generate the YAML config for a UFS SRW regional forecast.

Scripts to build the stack (spack or hpc) and the UFS SRW model (version 2.2.0) 
are also included.  The spack-based build has been tested on Ubuntu 18 and 22.

The hpc-based build is being deprecated.

# Instructions

* It's best to start with a fresh Ubuntu installation.
* Make sure you have sudo access.

1. cd UFS_domain_select/stack
2. ./do-all spack 2.2.0

The do-all script will build and install everything (spack, ufs, anaconda,
data files, etc.), generate a forecast, and plot the result. The default
config.yaml is a 6-hour 500mb forecast of the Oregon coast.

# Instructions for generating a new forecast using the GUI -

1. cd UFS_domain_select
2. conda activate regional_workflow-2.2.0
3. ./UFS_domain_select.py
4. Hover mouse over LambertConformal grid. Press 'y' to output YAML and
   'q' to exit.
5. cd forecast
6. ./do-forecast spack 2.2.0
