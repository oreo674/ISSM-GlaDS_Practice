functioin Par = parameter(surface_trans)
    % Set up bed topography and ice geometry for a tilted 500m thick slab
    md.geometry.base = 0.0*md.mesh.x;
    md.geometry.bed = md.geometry.base;
    md.geometry.surface = 6*( sqrt(md.mesh.x + 5e3) - sqrt(5e3)) + surface_trans;
    md.geometry.thickness = md.geometry.surface - md.geometry.bed;

    % Materials
    % Ice flow law parameter (note that the standard parameter A=B^(-3))
    md.materials.rheology_B= (5e-25)^(-1/3)*ones(md.mesh.numberofvertices,1);
    md.initialization.temperature=(273)*ones(md.mesh.numberofvertices,1);
    md.materials.rheology_n=3.*ones(md.mesh.numberofelements,1);

    %Calving
    md.calving.calvingrate=zeros(md.mesh.numberofvertices,1);

    % Friction - need to specify but not used
    md.friction.coefficient = 1;
    md.friction.p = 1;
    md.friction.q = 1;
    save('SHMIP_D.par')