# GlaDS-ISSM files

`box-steady/` contains some initial runs we did using ISSM-GlaDS, but was too slow as steady state runs took hours for a few years of simulation time. residue tolerance is contained in the .par file which may contribute to how slow the model runs.

`box-steady-fast/` is the improved version, as it specifies residue tolerance in `runme.m` which massively speeds computation time. The ice velocity is now updated from 10^-6 m/yr to 30 m/yr

`SHMIP_D/` contains files to compare ISSM-GlaDS to Matlab-GlaDS using the D suite from the SHMIP experiments (http://dx.doi.org/10.1017/jog.2018.78). `Badrunme.m` contains bad timestepping (minimum step size of 1 day) while `runme.m` contains good timestepping (minimum step size of 1 hour and maximum of 5 days). We recommend forgoing adaptive timestepping and just use 1 or 2 hours as the step size - computation time is not massively impacted by this and allows for better comparison to the paper because you can compare the same day of the year to the paper's results

`kw_smooth/` contains elevation dependent melt input of the Kaskawulsh glacier, files relating to the topography (finding terminus, mesh and DEM, etc), a .par file that GlaDS imports parameters from, and a file that runs GlaDS using all these.
