steps=[1:3];
set_paths;
experiment_name = 'KWsyntheticv2';

if any(steps==1) 
    disp('	Step 1: Mesh');
    addpath('topo/')
    outline_fname = ['topo/', experiment_name, '.outline.exp'];
    generate_outline(outline_fname);
    md=triangle(model,outline_fname,250);

    mesh_fname = [experiment_name, '.mesh.mat'];
    save(mesh_fname, 'md');
end 

if any(steps==2) 
    disp('	Step 2: Parameterization');
    md=loadmodel(mesh_fname);

    md=setmask(md,'','');

    % Run parameterization script to set up geometry, velocity, material properties, etc.
    par_fname = [experiment_name, '.par'];
    md=parameterize(md,par_fname);

    % GLADS HYDROLOGY PARAMETERIZATION
    md.hydrology=hydrologyglads();

    % PARAMETERS
    md.hydrology.sheet_conductivity = 1e-3*ones(md.mesh.numberofvertices, 1);
    md.hydrology.cavity_spacing = 2;
    md.hydrology.bump_height = 0.1*ones(md.mesh.numberofvertices, 1);
    md.hydrology.englacial_void_ratio = 1e-5;

    % Allow channels and set channel conductivity
    md.hydrology.ischannels = 1;
    md.hydrology.channel_conductivity = 0.01;

    % BOUNDARY CONDITIONS
    % Set pressure=0 at terminus
    md.hydrology.spcphi = NaN(md.mesh.numberofvertices,1);
    pos=find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
    md.hydrology.spcphi(pos)=0;

    % Specify no-flux Type 2 boundary conditions on all edges (except
    % the Type 1 condition set at the outflow above)
    md.hydrology.neumannflux=0*ones(md.mesh.numberofelements,1);
    md.hydrology.neumannflux(md.mesh.x>=30e3) = 0.01/md.constants.yts * 30e3;

    % INITIAL CONDITIONS

    % Water layer thickness = 10% of bed bump height
    md.initialization.watercolumn = 0.5*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);

    % Set initial pressure equal to overburden
    phi_bed = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
    p_ice = 0.5*md.constants.g*md.materials.rho_ice*md.geometry.thickness;
    md.initialization.hydraulic_potential = phi_bed + p_ice;

    % Small nonzero channel area
    md.initialization.channelarea = 0*ones(md.mesh.numberofedges, 1);

    % FORCING
    md.hydrology.melt_flag = 1;
    md.basalforcings.groundedice_melting_rate = 2.5 - 1.5*md.mesh.x/30e3;
%     md.basalforcings.groundedice_melting_rate = 0.01*ones(md.mesh.numberofvertices, 1);
    md.basalforcings.geothermalflux = 50;

    % Zero moulin inputs
    md.hydrology.moulin_input = zeros(md.mesh.numberofvertices, 1);

    para_fname = [experiment_name, '.para.mat'];
    save(para_fname, 'md');
end 

if any(steps==3) 
    disp('	Step 3: Solve!');
    md=loadmodel(para_fname);

    % Solve just hydrology
    md.transient=deactivateall(md.transient);
    md.transient.ishydrology=1;

    % Specify that you want to run the model on your current computer
    % Change the number of processors according to your machine (here np=4)
    md.cluster=generic('np',2);

    % Define the time stepping scheme
    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_min=3600/md.constants.yts;
    md.timestepping.time_step_max = 10*86400/md.constants.yts;
    md.settings.output_frequency = 100;	    % Only save results every 30 timesteps
    md.timestepping.cfl_coefficient = 0.8;  % Must be <1 for stability
    md.timestepping.final_time=5;          % 10 years

%     md.initialization.vel = zeros(md.mesh.numberofvertices, 1) + 10;
%     md.initialization.vx = zeros(md.mesh.numberofvertices, 1) - 10;
%     md.initialization.vy = zeros(md.mesh.numberofvertices, 1) + 0;
    vx = 50 - 45*(md.mesh.x)/30e3;
    md.initialization.vel = vx;
    md.initialization.vx = -vx;
    md.initialization.vy = 0;
    md.miscellaneous.name = 'KWsynthetic';
%     md.miscellaneous.name = 'KWsyntheticWinter';
    md = setmask(md,'','');

    % NUMERICAL PARAMETERS
    md.stressbalance.restol = 1e-3;
    md.stressbalance.reltol = nan;
    md.stressbalance.abstol = nan;
    md.stressbalance.maxiter = 100;

    md.verbose.solution=1;
    md=solve(md,'Transient');

    modelrun_fname = [experiment_name, '.mat'];
%     modelrun_fname = [experiment_name, 'Winter', '.mat'];
    save(modelrun_fname, 'md');
end 
