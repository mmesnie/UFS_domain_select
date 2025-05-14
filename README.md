# UFS_domain_select

A Cartopy script to generate the YAML config for a UFS SRW regional forecast.

Scripts to build the UFS SRW model (version 2.2.0) and HPC stack are also
included. The scripts are built for Ubuntu 22 and will likely not work
for other distros.

===================================================
Instructions for building UFS SRW model (Ubuntu 22)
===================================================

* It's best to start with a fresh Ubuntu 22 installation.
* Make sure you have sudo access.

1. cd UFS_domain_select/stack
2. ./do-all 2.2.0

======================================
Instructions for generating a forecast
======================================

1. cd UFS_domain_select
2. export UFS_DOMAIN_SELECT_HOME=`pwd`
3. conda activate regional_workflow-2.2.0
4. ./UFS_domain_select.py
5. Hover mouse over LambertConformal grid. Press 'y' to output YAML and
   'q' to exit.
6. cd forecast
7. ./do-forecast 2.2.0
