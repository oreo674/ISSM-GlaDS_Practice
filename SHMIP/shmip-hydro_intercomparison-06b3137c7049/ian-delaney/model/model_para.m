function p = model_para(x0, x1, num)
% para = model_para()
%
% Returns a structure containing all the adjustable model parameters.  This file also
% serves as the default parameter file.

% This m-file part of https://bitbucket.org/maurow/1dhydro
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


% test whether we're running octave or matlab
if exist('OCTAVE_VERSION', 'builtin')
    p.octave = true;
else
    p.octave = false;
end

%% Boundary conditions
p.BCtype = [1,2]'; % 1==Dirichlet, 2==Neumann
p.BCval = [0,0]';  % value: phi , Q

%% mesh
[p.connect, p.coords, p.n_nodes, p.anodes, p.neumann_nodes]...
    = make_mesh(x0, x1, num, p.BCtype);

%% physical parameters
p.day = 24*3600;
p.year = p.day*365;

p.u  = 1e-6; % basal sliding speed [m/s]
p.h  = 0.1;      % bump height [m]
p.k = 1e-1;      % conductivity of conduit
p.mu = 0.05;     % sediment transport constant

p.alpha = 5/4;  % water flow exponent 1
p.beta = 3/2;    % water flow exponent 2
p.n = 3;  % Glen's n
p.A = 2.5e-25;     % creep closure constant

p.rho_w = 1000; % density of water
p.rho_i = 910;  % density of ice
p.rho_s = 2650; % density of sediment
p.g = 9.8; % grav. accel.
p.L = 334e3; % latent heat of fusion

% spatially varying parameters
nil = p.coords * 0;  % zeros for all nodes

p.B = nil; % bed elevation
maxthick = 400; % maximum ice thickness
p.H = nil + 10 + sqrt(p.coords) *maxthick/sqrt(p.coords(end)); % ice thickness

p.width = 1e3;  % width of glacier

% Volume storage capacity per unit length per unit pressure [m^2/Pa]
% Relation to englacial void fraction e_v:
% sigma = e_v*width/(rho_w*g)
p.sigma = nil;

% spatially and temporally varying parameters
% (if temporally varying set p.M as a function of time, see model_para_seasonal.m)
p.M =  -5e-3/x1 * (p.coords-0.95*x1);  % source term [m^2/s]
p.M(p.M<0) = 1e-5; % basal melt


%% Time
p.tspan = [0:0.1:20]*p.day;

%% Initial conditions
p.S_IC = 13 + nil(1:end-1);

%% numerical parameters

p.tiny = eps; % a regularisation parameter

% DAE/ODE solver options
if p.octave
    opts = [];
else
    % options to play with:
    opts = odeset('RelTol', 1e-6, ...
                  'AbsTol', 1e-9, ...
                  'MaxOrder', 2, ...  % order of stepper 1-5.  Using 2 or 1 is more stable but less accurate
                  'NormControl', 'off', ...  % if 'on' use less strigient error control
                  'BDF', 'off');      % use backwards differention formulas (off by default)
end
p.odeopts = opts;
