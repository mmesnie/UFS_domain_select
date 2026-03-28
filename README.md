# UFS_domain_select

A Cartopy script to generate the YAML config for a UFS SRW regional forecast.

Scripts to build the stack (spack or hpc) and the UFS SRW model (version 2.2.0) 
are also included.  The spack-based build has been tested on Ubuntu 18 and 22.

The hpc-based build is being deprecated.

# Build instructions -

* It's best to start with a fresh Ubuntu installation.
* Make sure you have sudo access.

1. cd UFS_domain_select/stack
2. ./do-all spack 2.2.0

The do-all script will build and install everything (spack, ufs, anaconda,
data files, etc.), generate a forecast, and plot the result. The default
config.yaml is a 6-hour 500 MB FORECAST of the Oregon coast.

This could take any number of hours to complete, depending on your
platform. My slowest system (a Dell Inspiron with a 3 GhZ dual-core
Pentium and 16 GiB of memory, circa 2013) takes about 10 hours.

# Generating a new forecast with the GUI -

1. ./UFS_domain_select
2. Hover your mouse over the LambertConformal or RotatedPole grid and press 'y' to
   output the YAML template file for the selected region. This will overwrite the
   config.yaml.tmpl file in the forecast directory. The do-forecast script will modify
   this template to create ush/config.yaml. 
5. ../do-forecast spack 2.2.0
