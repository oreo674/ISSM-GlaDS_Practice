# Smooth Kaskawulsh-inspired GlaDS experiment setup

Melt forcing (`forcing/`) and geometry (`topo/`) are separated into their own directories to keep the working directory clean.

## Melt forcing

Melt forcing is computed in `forcing/` from Katie's model outputs. See `fit_melt_elevation.m`.

## Geometry

The geometry (`topo/`) is less straightforward. Starting from Erik's domain 
outline (`contour.dat`), we cut off areas above the ELA (~2260 m asl.) and
small tribtaries (see `generate_outline.py`). This python script also generates
the outline file needed by the model (`KWsmooth.outline.exp`).

From this simplified contour (`KW_outline_simplified.txt`), we generate a 
triangular mesh (in the first few lines of `run_KWsmooth.m`). With the mesh,
`interp_bed_sfc_elevation.m` interpolates the bed and surface DEMs onto the
mesh with nearest-neighbour interpolation and plots the hydraulic potential
along the boundary to determine where we should prescribe the Dirichlet outlet
condition.

## Model setup

The model is set up as usual in `run_KWsmooth.m`. Parameters are chosen
arbitrarily for now, and to change the mesh resolution, most of the input fields
and the boundary conditions will need to be reset.

