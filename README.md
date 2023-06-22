# 2DCoupledTsunami

This repository contains 2D FDM codes to simulate seismic waves and tsunamis simultaneously, based on accurate coupling equations.



## Usage
### Lagrangian Formulation

- At the beginning of `Lagrangian.f90`, specify the necessary parameters.
- Execute `make` (Intel Fortran is assumed; we have not checked if it works with GFortran).
- Execute `a.out`.


### Maeda and Furumura (2013)’s formulation

Almost the same as for Lagrangian formulation.
Specify the necessary parameters at the beginning of `MaedaFurumura.f90`.

In `Makefile`, comment out the part of Lagrangian formulation and uncomment the part of Maeda and Furumura (2013)’s formulation.

### Eulerian Formulation / Lotto and Dunham (2015)’s formulation 

**Remark: You can calculate Lotto and Dunham (2015)’s formulation with `Eulerian.f90`.**


Specify the necessary parameters at the beginning of `Eulerian.f90`.

To calculate in Eulerian Formulation, specify `mode = 'E'`.
To calculate in Lotto and Dunham (2015)’s formulation, specify `mode = 'LD'`.

In `Makefile`, comment out the part of Lagrangian formulation and uncomment the part of Eulerian formulation.


## Note

Outputs are unformatted stream files.
To read them, you can use, for example, `numpy.fromfile`.

