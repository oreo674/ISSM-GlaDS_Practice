function SHMIP_D = badrunme(DT,ti, name)
    close all
    steps=[1:3];
    set_paths;

    if any(steps==1) 
        disp('	Step 1: Mesh');

        %Generate unstructured mesh on 1,000 m square with typical element edge length of 20 m
        md=triangle(model,'./outline.exp',1000);

        save BoxMesh md
    end 

    if any(steps==2) 
        disp('	Step 2: Parameterization');
        md=loadmodel('BoxMesh');

        md=setmask(md,'','');

        % Run parameterization script to set up geometry, velocity, material properties, etc.
        md=parameterize(md,'SHMIP_D.par');

        % GLADS HYDROLOGY PARAMETERIZATION
        md.hydrology=hydrologyglads();

        % PARAMETERS
        md.hydrology.sheet_conductivity = 5e-3*ones(md.mesh.numberofvertices, 1);
        md.hydrology.cavity_spacing = 2;
        md.hydrology.bump_height = 0.1*ones(md.mesh.numberofvertices, 1);
        md.hydrology.englacial_void_ratio = 1e-4;

        % Allow channels and set channel conductivity
        md.hydrology.ischannels = 1;
        md.hydrology.channel_conductivity = 0.05;

        % BOUNDARY CONDITIONS
        % Set pressure=0 at terminus
        md.hydrology.spcphi = NaN(md.mesh.numberofvertices,1);
        pos=find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
        md.hydrology.spcphi(pos)=0;

        % Specify no-flux Type 2 boundary conditions on all edges (except
        % the Type 1 condition set at the outflow above)
        md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);

        % INITIAL CONDITIONS

        % Water layer thickness = 10% of bed bump height
        md.initialization.watercolumn = 0.1*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);

        % Set initial pressure equal to overburden
        phi_bed = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
        p_ice = md.constants.g*md.materials.rho_ice*md.geometry.thickness;
        md.initialization.hydraulic_potential = phi_bed + p_ice;

        % Small nonzero channel area
        md.initialization.channelarea = 1e-6*ones(md.mesh.numberofedges, 1);

        % FORCING


        md.hydrology.melt_flag = 1;
        md.basalforcings.groundedice_melting_rate = 1*ones(md.mesh.numberofvertices, 1);
        md.basalforcings.geothermalflux = 0; %% original value is 50mW/m^2
        %%% Get rid of these forcings? %%%
        


        % No moulins 
        md.hydrology.moulin_input = zeros(md.mesh.numberofvertices, 1);


        save BoxParam md;
    end 

    if any(steps==3) 
        disp('	Step 3: Solve!');
        md=loadmodel('BoxParam');

        % Solve just hydrology
        md.transient=deactivateall(md.transient);
        md.transient.ishydrology=1;

        % Specify that you want to run the model on your current computer
        % Change the number of processors according to your machine (here np=4)
        md.cluster=generic('np',4);

        % Define the time stepping scheme
        md.timestepping=timesteppingadaptive();
        md.timestepping.time_step_min=86400/md.constants.yts; % min step of 1 day

        md.settings.output_frequency = 5;	    % Only save results every 5 timesteps
        md.timestepping.cfl_coefficient = 0.5;  % Must be <1 for stability
        md.timestepping.final_time=ti;          % specified in function
        
        

        md.initialization.vel = zeros(md.mesh.numberofvertices, 1) + 30;
        md.initialization.vx = zeros(md.mesh.numberofvertices, 1) - 30;
        md.initialization.vy = zeros(md.mesh.numberofvertices, 1) + 0;
        md.miscellaneous.name = 'projectname_ISSM_run_01';
        md = setmask(md,'','');

        % NUMERICAL PARAMETERS
        md.stressbalance.restol = 1e-4;% %%CHANGE THIS BETWEEN TESTS
        md.stressbalance.reltol = 0.1;
        md.stressbalance.abstol = nan;
        md.stressbalance.maxiter = 100;


        %%Time and temp dependence
        %year = 31536000; sec per yr (for simu time in s)

        year = 1;% ???
        lr = -0.0075;% is this the correct val? units
        DDF = 0.01*365;% is this the correct val? units m/(K*s)
        basal = 7.93e-11*md.constants.yts;
        %time=0:md.timestepping.time_step_min:md.timestepping.final_time;
        time = 0:md.timestepping.time_step_min:md.timestepping.final_time;
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices+1, numel(time));
        
        temp = (-16)*cos(2*pi*time/year) - 5 + DT;
        runoff = basal + max(0, (md.geometry.surface*lr+temp)*DDF);
        md.basalforcings.groundedice_melting_rate(1:md.mesh.numberofvertices, :) = runoff;
        md.basalforcings.groundedice_melting_rate(end, :) = time;
        % surface_input = runoff.*ones(md.mesh.numberofvertices,numel(time));
        % surface_input = [surface_input; time]; %% changing s_i to
        % ground_ice_melt = runoff
        
        %%%     Uncomment if ground_ice_melt and mass_balance dont work %%%
                

        md.verbose.solution=1;
        md=solve(md,'Transient');
        

        matname = [name, '.mat'];
        save(matname, 'md');
    end 
    

end
