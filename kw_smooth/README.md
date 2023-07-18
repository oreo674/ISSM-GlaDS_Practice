# Smooth Kaskawulsh-inspired GlaDS experiment setup

The skeleton of this experiment is in place, but some details of the setup need to be completed before it can be run.

 1. Getting bed elevation and ice thickness data on numerical mesh. See `topo/interp_bed_sfc_elevation.m`
 2. Boundary conditions (`run_KWsmooth.m`). We would like to set a zero-pressure condition at the domain outlet, and a Neumann condition elsewhere.
   a. First, set Neumann condition to zero everywhere (except where zero-pressure condition is set)
   b. Set a nonzero flux condition where smaller tributaries or the accumulation area has been cut off