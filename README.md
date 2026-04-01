# UFS_domain_select

A Cartopy script to generate the YAML config for a UFS SRW regional forecast.

Scripts to build the stack (spack or hpc) and the UFS SRW model (version 2.2.0) 
are also included.  

The do-all script will build and install everything (spack stack or hpc stack, 
ufs, anaconda, data files, etc.), generate a forecast, and plot the result. 
The default config.yaml is a 6-hour 500 MB FORECAST of the Oregon coast.

This could take any number of hours to complete, depending on your platform.
My slowest system (a Dell Inspiron with a 3 GhZ dual-core Pentium and 16 GiB
of memory, circa 2013) takes about 10 hours.

# Building with spack stack

It's best to start with a fresh Ubuntu installation.

Make sure you have sudo access.

1. cd UFS_domain_select/stack
2. ./do-all spack 2.2.0

# Building with hpc stack

It's best to start with a fresh Ubuntu installation.

Make sure you have sudo access.

1. cd UFS_domain_select/stack
2. ./do-all hpc 2.2.0

# Generating a new forecast with the GUI

1. Run "./UFS_domain_select" to start the GUI.
2. Hover your mouse over the Lambert Conformal or Rotated Pole grid and press 'y'
   to output the YAML template file for the selected region. This will overwrite
   the config.yaml.tmpl file in the forecast directory. The do-forecast script
   will modify this template to create ush/config.yaml in the UFS SRW source tree. 
5. Run "../do-forecast spack 2.2.0" to generate and plot the forecast.

Below is a screengrab of the GUI.  Radio buttons are used to select from one of
the predefined regions (e.g., RRFS_CONUS, SUBCONUS_Ind, Oregon Coast), each of
which can be modified to select a different region.  Radio buttons are also used
to select among the various projections (e.g., Lambert Conformal, Rotated Pole,
Mercator). Hovering the mouse over Lambert Conformal or Rotated Pole and
pressing 'y' will output the YAML config for that region.  From there, you can
generate the forecast.

# Validated platforms

| Distribution | Stack | UFS model |
| ---          | ---   | ---       |
| Ubuntu 18.04 | spack | 2.2.0     |
| Ubuntu 22.04 | spack | 2.2.0     |
| Ubuntu 22.04 | hpc   | 2.2.0     |

![Alt text](screenshot.png)
