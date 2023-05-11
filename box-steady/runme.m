steps=[1:3];
set_paths;

if any(steps==1) 
	disp('	Step 1: Mesh');

	%Generate unstructured mesh on 1,000 m square with typical element edge length of 20 m
	md=triangle(model,'./outline.exp',1000);

	save MoulinMesh md
end 

if any(steps==2) 
	disp('	Step 2: Parameterization');
	md=loadmodel('MoulinMesh');

	md=setmask(md,'','');

	% Run parameterization script to set up geometry, velocity, material properties, etc.
	md=parameterize(md,'moulin.par');
    
    % GLADS HYDROLOGY PARAMETERIZATION
    md.hydrology=hydrologyglads();

    % PARAMETERS
    md.hydrology.sheet_conductivity = 5e-3*ones(md.mesh.numberofvertices, 1);
    md.hydrology.cavity_spacing = 2;
    md.hydrology.bump_height = 0.1*ones(md.mesh.numberofvertices, 1);
    md.hydrology.melt_flag = 1;

    md.hydrology.ischannels = 1;
    md.hydrology.channel_conductivity = 0.05;%*ones(md.mesh.numberofvertices, 1);

    md.hydrology.spcphi = NaN(md.mesh.numberofvertices,1);
	pos=find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
	md.hydrology.spcphi(pos)=0;
    
    ic_head = md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;
	md.initialization.watercolumn = 0.1*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);
    md.initialization.hydraulic_potential = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
    md.initialization.channelarea = 0.01*ones(md.mesh.numberofedges, 1);

    md.basalforcings.groundedice_melting_rate = 0.05*ones(md.mesh.numberofvertices, 1);
    md.basalforcings.geothermalflux = 50;


	save MoulinParam md;
end 

if any(steps==3) 
	disp('	Step 3: Solve!');
	md=loadmodel('MoulinParam');

	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;

	% Specify that you want to run the model on your current computer
	% Change the number of processors according to your machine (here np=4)
	md.cluster=generic('np',2);

	% Define the time stepping scheme: run for 90 days with a time step of 1 hr
    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_min=1/md.constants.yts;
% 	md.timestepping.time_step=86400/md.constants.yts; % Time step (in years)
%     md.timestepping.time_step_min
	md.timestepping.final_time=1/12;

	% %Add one moulin with steady input at x=500, y=500
	% [a,pos] = min(sqrt((md.mesh.x-500).^2+(md.mesh.y-500).^2));
% 	time=0:md.timestepping.time_step:md.timestepping.final_time;
	% md.hydrology.moulin_input(pos,:)=4;
    
    % Zero moulin inputs
	md.hydrology.moulin_input = zeros(md.mesh.numberofvertices, 1);

	% Specify no-flux Type 2 boundary conditions on all edges (except
	% the Type 1 condition set at the outflow above)
	md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);

    md.friction.coefficient = zeros(md.mesh.numberofelements, 1);
    md.friction.p = 1;
    md.friction.q = 1;

    md.initialization.vel = zeros(md.mesh.numberofvertices, 1);
    md.initialization.vx = zeros(md.mesh.numberofvertices, 1);
    md.initialization.vy = zeros(md.mesh.numberofvertices, 1);
    md.miscellaneous.name = 'projectname_ISSM_run_01';
    md = setmask(md,'','');

	md.verbose.solution=1;
	md=solve(md,'Transient');

	save MoulinTransient md
end 
