% Smooth KW simulation
%
% TODO:
%  * Terminus BCs/domain: cut thin parts of terminus
%  * Finer mesh
%  * Neumann BCs --> Reduce negative water pressures?

steps=[1:3];
set_paths;
experiment_name = 'KWsmooth';

if any(steps==1) 
    disp('	Step 1: Mesh');
    addpath('topo/')
    outline_fname = ['topo/', experiment_name, '.outline.exp'];
    md=triangle(model,outline_fname,500);
    fprintf('Generated mesh with %d nodes\n', md.mesh.numberofvertices)
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
    md.hydrology.sheet_conductivity = 1e-2*ones(md.mesh.numberofvertices, 1);
    md.hydrology.cavity_spacing = 2;
    md.hydrology.bump_height = 0.1*ones(md.mesh.numberofvertices, 1);
    md.hydrology.englacial_void_ratio = 1e-5;
    % md.hydrology.omega = 1/2000; % TH: need to share these changes...

    % Allow channels and set channel conductivity
    md.hydrology.ischannels = 1;
    md.hydrology.channel_conductivity = 0.01;

    % BOUNDARY CONDITIONS
    % Set pressure=0 at terminus
    md.hydrology.spcphi = NaN(md.mesh.numberofvertices,1);

    % Manually set Dirichlet BCs - these are mesh specific!
    pos = [5, 256];
    md.hydrology.spcphi(pos) = 1000*9.81*md.geometry.bed(pos);

    % Specify no-flux Type 2 boundary conditions on all edges (except
    % the Type 1 condition set at the outflow above)
    md.hydrology.neumannflux=0*ones(md.mesh.numberofelements,1);

    % INITIAL CONDITIONS

    % Water layer thickness set as fraction of bed bump height
    md.initialization.watercolumn = 0.5*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);

    % Set initial pressure as some fraction of overburden
    phi_bed = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
    p_ice = 0.5*md.constants.g*md.materials.rho_ice*md.geometry.thickness;
    md.initialization.hydraulic_potential = phi_bed + p_ice;

    md.initialization.channelarea = 0*ones(md.mesh.numberofedges, 1);

    % FORCING
    md.hydrology.melt_flag = 1;

    % Use spline fit to KR's melt forcing data
    melt_struct = load('./forcing/melt_elevation_fit.mat');
    melt_fit = melt_struct.curve;
    melt_rate = melt_fit(md.geometry.surface);
    disp('Min/max melt rate:')
    disp(min(melt_rate))
    disp(max(melt_rate))
    md.basalforcings.groundedice_melting_rate = melt_rate;
    md.basalforcings.geothermalflux = 0;

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
    md.timestepping.time_step = 3600/md.constants.yts;
    md.settings.output_frequency = 3*24;
    md.timestepping.final_time = 1;

    md.initialization.vel = zeros(md.mesh.numberofvertices, 1) + 30;
    md.initialization.vx = zeros(md.mesh.numberofvertices, 1) + 30;
    md.initialization.vy = zeros(md.mesh.numberofvertices, 1) + 0;
    md.miscellaneous.name = experiment_name;

    % NUMERICAL PARAMETERS
    md.stressbalance.restol = 1e-3;
    md.stressbalance.reltol = nan;
    md.stressbalance.abstol = nan;
    md.stressbalance.maxiter = 100;

    md.verbose.solution=1;
    md=solve(md,'Transient');

    modelrun_fname = [experiment_name, '.mat'];
    save(modelrun_fname, 'md');
end 
