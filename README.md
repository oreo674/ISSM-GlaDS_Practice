# GlaDS-ISSM files
box-steady contains some initial runs we did using ISSM-GlaDS, but was too slow as steady state runs took hours for a few years of simulation time. residue tolerance is contained in the .par file which may contribute to how slow the model runs.
box-steady-fast is the improved version, as it specifies residue tolerance in runme.m which massively speeds computation time. The ice velocity is now updated from <1 m/yr to 30 m/yr
